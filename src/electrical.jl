"""Lattice dielectric-breakdown growth with quasi-static field screening.

Conductive occupied cells are held at potential one, and the external box at zero.
Every `event_interval` units of process time, one face-connected frontier cell is
selected with probability proportional to its local voltage drop raised to `eta`.
Each event transfers `deposit_per_cell` nutrient into deposit. This is a CPU-only,
nondimensional dielectric-breakdown surrogate, not a calibrated electrical model.
The random variate is a deterministic hash of `(seed, event index)`.
"""
mutable struct DielectricGrowth
    potential::Array{Float32,3}
    occupied::BitArray{3}
    frontier::BitArray{3}
    seed::UInt64
    events::Int
    elapsed::Float64
    event_interval::Float64
    eta::Float64
    deposit_per_cell::Float32
    tolerance::Float64
    max_iterations::Int
    relaxation::Float64
    last_residual::Float64
    last_iterations::Int
    reached_ground::Bool
end

function DielectricGrowth(w::World; seed=w.seed, root=ntuple(d->cld(w.grid.size[d],2),3),
    event_interval=0.05, eta=2.5, deposit_per_cell=0.8, tolerance=2e-5,
    max_iterations=1200, relaxation=1.8)
    w.backend isa CPU || throw(ArgumentError("DielectricGrowth currently supports CPU worlds"))
    all(>=(5),w.grid.size) || throw(ArgumentError("electrical growth needs a 3D grid at least 5 cells per axis"))
    all(d->1<root[d]<w.grid.size[d],1:3) || throw(ArgumentError("root must lie inside the grounded box"))
    event_interval>0 && isfinite(event_interval) || throw(ArgumentError("positive finite event_interval required"))
    eta>0 && isfinite(eta) || throw(ArgumentError("positive finite eta required"))
    0<deposit_per_cell<=minimum(w.fields.nutrient) || throw(ArgumentError("deposit_per_cell exceeds available nutrient"))
    tolerance>0 && isfinite(tolerance) && max_iterations>0 && 0<relaxation<2 || throw(ArgumentError("invalid Laplace solver controls"))
    p=DielectricGrowth(zeros(Float32,w.grid.size),falses(w.grid.size),falses(w.grid.size),
        reinterpret(UInt64,Int64(seed)),0,0.,Float64(event_interval),Float64(eta),Float32(deposit_per_cell),
        Float64(tolerance),Int(max_iterations),Float64(relaxation),Inf,0,false)
    electrical_mark!(p,CartesianIndex(root))
    solve_potential!(p,w.grid.spacing)
    w.fields.nutrient[CartesianIndex(root)]-=p.deposit_per_cell
    w.fields.deposit[CartesianIndex(root)]+=p.deposit_per_cell
    p
end

@inline electrical_boundary(I,dims)=any(d->I[d]==1 || I[d]==dims[d],1:3)
@inline function electrical_neighbors(I)
    (I+CartesianIndex(1,0,0),I-CartesianIndex(1,0,0),
     I+CartesianIndex(0,1,0),I-CartesianIndex(0,1,0),
     I+CartesianIndex(0,0,1),I-CartesianIndex(0,0,1))
end
function electrical_occupy!(w,p,I)
    p.occupied[I] && throw(ArgumentError("site is already occupied"))
    food=w.fields.nutrient[I]
    food>=p.deposit_per_cell || throw(ArgumentError("insufficient nutrient at growth site"))
    w.fields.nutrient[I]-=p.deposit_per_cell
    w.fields.deposit[I]+=p.deposit_per_cell
    electrical_mark!(p,I)
end
function electrical_mark!(p,I)
    p.occupied[I]=true;p.frontier[I]=false;p.potential[I]=1f0
    dims=size(p.potential)
    for J in electrical_neighbors(I)
        if electrical_boundary(J,dims)
            p.reached_ground=true
        elseif !p.occupied[J]
            p.frontier[J]=true
        end
    end
end

"""Maximum residual of the normalized discrete Laplace equation at free interior cells."""
function potential_residual(p::DielectricGrowth,spacing=(1.,1.,1.))
    v=p.potential;weights=ntuple(d->1/Float64(spacing[d])^2,3);denom=2sum(weights);r=0.
    nx,ny,nz=size(v)
    @inbounds for k in 2:nz-1,j in 2:ny-1,i in 2:nx-1
        p.occupied[i,j,k] && continue
        average=(weights[1]*(v[i-1,j,k]+v[i+1,j,k])+weights[2]*(v[i,j-1,k]+v[i,j+1,k])+weights[3]*(v[i,j,k-1]+v[i,j,k+1]))/denom
        r=max(r,abs(average-v[i,j,k]))
    end
    r
end

"""SOR solve at fixed conductor geometry; warm starts retain the previous field."""
function solve_potential!(p::DielectricGrowth,spacing=(1.,1.,1.))
    v=p.potential;weights=ntuple(d->1/Float64(spacing[d])^2,3);denom=2sum(weights)
    nx,ny,nz=size(v)
    # Boundary values are explicit Dirichlet data; never periodic wraps.
    v[1,:,:].=0;v[end,:,:].=0;v[:,1,:].=0;v[:,end,:].=0;v[:,:,1].=0;v[:,:,end].=0
    v[p.occupied].=1
    for iteration in 1:p.max_iterations
        @inbounds for k in 2:nz-1,j in 2:ny-1,i in 2:nx-1
            p.occupied[i,j,k] && continue
            average=(weights[1]*(v[i-1,j,k]+v[i+1,j,k])+weights[2]*(v[i,j-1,k]+v[i,j+1,k])+weights[3]*(v[i,j,k-1]+v[i,j,k+1]))/denom
            v[i,j,k]=v[i,j,k]+p.relaxation*(average-v[i,j,k])
        end
        # Checking every eight sweeps amortizes the residual pass.
        if iteration%8==0 || iteration==p.max_iterations
            p.last_residual=potential_residual(p,spacing);p.last_iterations=iteration
            p.last_residual<=p.tolerance && return p
        end
    end
    error("electrical potential failed to converge: residual $(p.last_residual) > $(p.tolerance)")
end

@inline function electrical_uniform(seed::UInt64,event::Int)
    # SplitMix64 finalizer: explicit UInt64 arithmetic wraps by definition.
    z=seed+UInt64(event)*0x9e3779b97f4a7c15
    z=(z ⊻ (z>>30))*0xbf58476d1ce4e5b9
    z=(z ⊻ (z>>27))*0x94d049bb133111eb
    z=z ⊻ (z>>31)
    Float64(z>>11)*0x1.0p-53
end

function electrical_event!(w,p::DielectricGrowth)
    p.reached_ground && return false
    candidates=CartesianIndex{3}[];weights=Float64[]
    for I in CartesianIndices(p.frontier)
        p.frontier[I] || continue
        w.fields.nutrient[I]>=p.deposit_per_cell || continue
        # Field magnitude at a conductor face, summed over touching occupied faces.
        weight=0.
        for (axis,h) in enumerate(w.grid.spacing)
            offset=CartesianIndex(ntuple(d->d==axis ? 1 : 0,3))
            for J in (I-offset,I+offset)
                p.occupied[J] && (weight+=(max(0.,1-Float64(p.potential[I]))/h)^p.eta)
            end
        end
        push!(candidates,I);push!(weights,weight)
    end
    total=sum(weights)
    total>0 || return false
    target=electrical_uniform(p.seed,p.events+1)*total;acc=0.;selected=last(candidates)
    for (I,weight) in zip(candidates,weights)
        acc+=weight
        if acc>target;selected=I;break;end
    end
    electrical_occupy!(w,p,selected);p.events+=1
    solve_potential!(p,w.grid.spacing)
    w.fields.activity.=p.frontier .* max.(0f0,1f0.-p.potential)
    true
end

function validate(w::World,p::DielectricGrowth)
    w.backend isa CPU || throw(ArgumentError("DielectricGrowth requires CPU"))
    all(a->size(a)==w.grid.size,(p.potential,p.occupied,p.frontier)) || throw(DimensionMismatch("electrical state and world grid differ"))
    all(isfinite,p.potential) || throw(ArgumentError("nonfinite electrical potential"))
    isfinite(p.elapsed) && p.elapsed>=0 && p.events>=0 || throw(ArgumentError("invalid electrical event clock"))
    isfinite(p.event_interval) && p.event_interval>0 || throw(ArgumentError("invalid electrical event interval"))
    isfinite(p.eta) && p.eta>0 || throw(ArgumentError("invalid field exponent"))
    isfinite(p.deposit_per_cell) && p.deposit_per_cell>0 || throw(ArgumentError("invalid electrical deposit budget"))
    isfinite(p.tolerance) && p.tolerance>0 && p.max_iterations>0 && isfinite(p.relaxation) && 0<p.relaxation<2 || throw(ArgumentError("invalid electrical solver controls"))
end
function advance!(w::World,p::DielectricGrowth,dt::Real)
    isfinite(dt)&&dt>0 || throw(ArgumentError("positive finite dt required"))
    validate(w,p)
    p.elapsed+=Float64(dt)
    desired=floor(Int,(p.elapsed+4eps(Float32)*max(1.,p.elapsed))/p.event_interval)
    while p.events<desired && !p.reached_ground
        electrical_event!(w,p) || break
    end
    w
end

"""Optical map: real deposit is the channel; a one-cell six-neighbor halo is optical only."""
function electrical_optical_fields(w::World,p::DielectricGrowth; halo=0.18f0,carrier=0.08f0)
    deposit=copy(Array(w.fields.deposit));glow=copy(deposit);nx,ny,nz=size(glow)
    @inbounds for k in 2:nz-1,j in 2:ny-1,i in 2:nx-1
        glow[i,j,k]+=halo*(deposit[i-1,j,k]+deposit[i+1,j,k]+deposit[i,j-1,k]+deposit[i,j+1,k]+deposit[i,j,k-1]+deposit[i,j,k+1])/6
    end
    mask=ones(Float32,size(deposit))
    for I in CartesianIndices(mask)
        xyz=ntuple(d->2f0*(I[d]-.5f0)/size(mask,d)-1f0,3)
        r=sqrt(sum(x->x*x,xyz));mask[I]=clamp((.94f0-r)*30f0,0f0,1f0)
    end
    (density=carrier.+0.35f0.*glow,activity=2.5f0.*glow,structure=deposit,mask=mask)
end
