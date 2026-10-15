"""CPU central-force graph carrier. Columns are persistent node IDs; bonds never reorder.
Positions use world physical coordinates. Springs use force/length stiffness,
node mass, time and dimensionless tensile failure strain. No periodic wrapping.
"""
mutable struct BondNetwork
    position::Matrix{Float64}
    velocity::Matrix{Float64}
    reference::Matrix{Float64}
    node_ids::Vector{Int}
    bonds::Matrix{Int}
    bond_ids::Vector{Int}
    rest_length::Vector{Float64}
    stiffness::Vector{Float64}
    threshold::Vector{Float64}
    damage::Vector{Float64}
    peak_strain::Vector{Float64}
    mass::Vector{Float64}
    grip::Vector{Int}
    load::Matrix{Float64}
    damping::Float64
    stretch_rate::Float64
    shear_rate::Float64
    elapsed::Float64
    damage_energy::Float64
    damping_energy::Float64
    external_work::Float64
end

function BondNetwork(w::World; size=(10,10,w.grid.size[3]==1 ? 1 : 10),
    extent=Tuple(.65 .* w.grid.extent), stiffness=8., mass=.02, damping=2.,
    failure_strain=.18, weak_seam=.48, stretch_rate=0., shear_rate=0., jitter=.08, seed=w.seed)
    w.backend isa CPU || throw(ArgumentError("BondNetwork requires CPU"))
    length(size)==3 && all(>(0),size) && size[1]>=2 && size[2]>=2 || throw(ArgumentError("invalid network size"))
    all(x->isfinite(x)&&x>0,extent) && stiffness>0 && mass>0 && damping>=0 && failure_strain>0 && 0<weak_seam<=1 && 0<=jitter<.3 || throw(ArgumentError("invalid mechanical parameters"))
    all(isfinite,(stiffness,mass,damping,failure_strain,stretch_rate,shear_rate,jitter)) || throw(ArgumentError("finite parameters required"))
    dims=Tuple(Int.(size)); inds=LinearIndices(dims); rng=MersenneTwister(seed)
    x=zeros(3,prod(dims));grip=zeros(Int,prod(dims))
    h=ntuple(d->dims[d]>1 ? Float64(extent[d])/(dims[d]-1) : 0.,3)
    for I in CartesianIndices(dims)
        i=inds[I]
        for d in 1:3
            x[d,i]=dims[d]>1 ? (I[d]-1)*h[d]-extent[d]/2 : 0.
            # Grip planes stay flat; interior is quenched geometric disorder.
            if 1<I[1]<dims[1];x[d,i]+=jitter*h[d]*(rand(rng)-.5);end
        end
        grip[i]=I[1]==1 ? -1 : I[1]==dims[1] ? 1 : 0
    end
    edges=Tuple{Int,Int}[]
    for I in CartesianIndices(dims), dz in -1:1,dy in -1:1,dx in -1:1
        J=I+CartesianIndex(dx,dy,dz)
        all(d->1<=J[d]<=dims[d],1:3) || continue
        inds[J]>inds[I] || continue
        push!(edges,(inds[I],inds[J]))
    end
    b=hcat(([i,j] for (i,j) in edges)...); lengths=[norm(x[:,j]-x[:,i]) for (i,j) in edges]
    thresholds=fill(Float64(failure_strain),length(edges))
    # A weakened half-plane is a notch seed, not a painted crack: only spring
    # strength changes. Subsequent topology follows measured tensile strain.
    for (e,(i,j)) in enumerate(edges)
        if x[1,i]*x[1,j]<=0 && (x[2,i]+x[2,j])/2<0
            thresholds[e]*=weak_seam
        end
    end
    BondNetwork(copy(x),zeros(3,prod(dims)),x,collect(1:prod(dims)),b,
      collect(1:length(edges)),lengths,fill(Float64(stiffness),length(edges)),thresholds,
      zeros(length(edges)),zeros(length(edges)),fill(Float64(mass),prod(dims)),grip,
      zeros(3,prod(dims)),Float64(damping),Float64(stretch_rate),Float64(shear_rate),0.,0.,0.,0.)
end

@inline _bond_delta(p::BondNetwork,i,j) =
    (p.position[1,j]-p.position[1,i],p.position[2,j]-p.position[2,i],p.position[3,j]-p.position[3,i])

function bond_forces(p::BondNetwork)
    f=zeros(size(p.position))
    for e in eachindex(p.bond_ids)
        p.damage[e]==1 && continue
        i=p.bonds[1,e];j=p.bonds[2,e];delta=_bond_delta(p,i,j);l=norm(delta)
        l>1e-14 || continue
        scale=(1-p.damage[e])*p.stiffness[e]*(l-p.rest_length[e])/l
        for d in 1:3
            force=scale*delta[d];f[d,i]+=force;f[d,j]-=force
        end
    end
    f
end
function mechanical_energy(p::BondNetwork)
    elastic=0.
    for e in eachindex(p.bond_ids)
        i=p.bonds[1,e];j=p.bonds[2,e];extension=norm(_bond_delta(p,i,j))-p.rest_length[e]
        elastic+=.5*(1-p.damage[e])*p.stiffness[e]*extension^2
    end
    kinetic=sum(.5*p.mass[i]*sum(abs2,(p.velocity[1,i],p.velocity[2,i],p.velocity[3,i])) for i in eachindex(p.mass))
    (elastic=elastic,kinetic=kinetic,damage=p.damage_energy,damping=p.damping_energy,
     external_work=p.external_work,total=elastic+kinetic+p.damage_energy+p.damping_energy)
end
function break_bonds!(p::BondNetwork)
    for e in eachindex(p.bond_ids)
        i=p.bonds[1,e];j=p.bonds[2,e];extension=norm(_bond_delta(p,i,j))-p.rest_length[e]
        strain=extension/p.rest_length[e];p.peak_strain[e]=max(p.peak_strain[e],strain)
        if p.damage[e]<1 && strain>=p.threshold[e]
            p.damage_energy+=.5*(1-p.damage[e])*p.stiffness[e]*extension^2
            p.damage[e]=1.
        end
    end
    p
end
function stability_limit(w::World,p::BondNetwork)
    # Gershgorin bound on the current central-force tangent, including compressed
    # transverse stiffness; conservative safety margin for velocity Verlet.
    rows=zeros(length(p.mass))
    for e in eachindex(p.bond_ids)
        p.damage[e]==1 && continue
        i=p.bonds[1,e];j=p.bonds[2,e];l=norm(_bond_delta(p,i,j))
        tangent=p.stiffness[e]*(1-p.damage[e])*max(1.,abs(1-p.rest_length[e]/max(l,1e-14)))
        rows[i]+=tangent;rows[j]+=tangent
    end
    bound=maximum(2 .* rows ./ p.mass)
    bound==0 ? Inf : .35/sqrt(bound)
end
function validate(w::World,p::BondNetwork)
    w.backend isa CPU || throw(ArgumentError("BondNetwork requires CPU"))
    n=length(p.node_ids);b=length(p.bond_ids)
    size(p.position)==size(p.velocity)==size(p.reference)==size(p.load)==(3,n) || throw(DimensionMismatch("node arrays must have shape (3,n)"))
    size(p.bonds)==(2,b) || throw(DimensionMismatch("bonds must have shape (2,b)"))
    all(a->length(a)==n,(p.mass,p.grip)) && all(a->length(a)==b,(p.rest_length,p.stiffness,p.threshold,p.damage,p.peak_strain)) || throw(DimensionMismatch("mechanical property lengths differ"))
    length(unique(p.node_ids))==n && length(unique(p.bond_ids))==b || throw(ArgumentError("mechanical IDs must be unique"))
    all(i->1<=i<=n,p.bonds) && all(e->p.bonds[1,e]!=p.bonds[2,e],1:b) || throw(ArgumentError("invalid edge endpoint"))
    all(a->all(isfinite,a),(p.position,p.velocity,p.reference,p.load,p.mass,p.rest_length,p.stiffness,p.threshold,p.damage,p.peak_strain)) || throw(ArgumentError("nonfinite mechanical state"))
    all(>(0),p.mass) && all(>(0),p.rest_length) && all(>=(0),p.stiffness) && all(>(0),p.threshold) && all(x->0<=x<=1,p.damage) && all(x->x in (-1,0,1),p.grip) || throw(ArgumentError("invalid mechanical material state"))
    all(isfinite,(p.damping,p.stretch_rate,p.shear_rate,p.elapsed,p.damage_energy,p.damping_energy,p.external_work)) && p.damping>=0 && p.elapsed>=0 && p.damage_energy>=0 && p.damping_energy>=0 || throw(ArgumentError("invalid mechanical controls or ledgers"))
end
function mechanical_grips!(p,t)
    for i in eachindex(p.grip)
        s=p.grip[i];s==0 && continue
        for d in 1:3;p.position[d,i]=p.reference[d,i];p.velocity[d,i]=0.;end
        p.position[1,i]+=s*p.stretch_rate*t/2
        p.position[2,i]+=s*p.shear_rate*t/2
    end
end
function advance!(w::World,p::BondNetwork,dt::Real)
    check_dt(w,p,dt);dt=Float64(dt)
    # Grip displacement is a kinematic external operation; book its exact change
    # in elastic energy separately from integration error and damage loss.
    before=mechanical_energy(p).elastic
    mechanical_grips!(p,p.elapsed+dt)
    p.external_work+=mechanical_energy(p).elastic-before
    break_bonds!(p)
    f=bond_forces(p); decay=exp(-p.damping*dt/2);dx=Vector{Float64}(undef,3)
    for i in eachindex(p.mass)
        p.grip[i]!=0 && continue
        old=sum(abs2,(p.velocity[1,i],p.velocity[2,i],p.velocity[3,i]))
        for d in 1:3;p.velocity[d,i]*=decay;end
        p.damping_energy+=.5*p.mass[i]*(old-sum(abs2,(p.velocity[1,i],p.velocity[2,i],p.velocity[3,i])))
        for d in 1:3
            p.velocity[d,i]+=dt/2*(f[d,i]+p.load[d,i])/p.mass[i]
            dx[d]=dt*p.velocity[d,i];p.position[d,i]+=dx[d]
        end
        # Retain the same strided BLAS dot reduction as the original dense
        # slices; a tuple dot changes rounding on some supported platforms.
        p.external_work+=dot(view(p.load,:,i),dx)
    end
    break_bonds!(p);f=bond_forces(p)
    for i in eachindex(p.mass)
        p.grip[i]!=0 && continue
        for d in 1:3;p.velocity[d,i]+=dt/2*(f[d,i]+p.load[d,i])/p.mass[i];end
        old=sum(abs2,(p.velocity[1,i],p.velocity[2,i],p.velocity[3,i]))
        for d in 1:3;p.velocity[d,i]*=decay;end
        p.damping_energy+=.5*p.mass[i]*(old-sum(abs2,(p.velocity[1,i],p.velocity[2,i],p.velocity[3,i])))
    end
    p.elapsed+=dt;w
end
function record_process_fields(p::BondNetwork)
    names=fieldnames(typeof(p));NamedTuple{names}(Tuple(getfield(p,n) for n in names))
end

"""Rasterize actual current bonds and fracture endpoints; no rest-lattice proxy.
A scar is the broken-bond endpoint neighborhood, not a repaired mechanical bond.
"""
function mechanical_optical_fields(w::World,p::BondNetwork; radius=.012)
    radius>0 || throw(ArgumentError("positive radius required"))
    density=zeros(Float32,w.grid.size);stress=similar(density);fill!(stress,0);scar=copy(stress)
    h=minimum(w.grid.spacing[1:(w.grid.size[3]==1 ? 2 : 3)])
    function splat!(a,x,value)
        index=ntuple(d->floor(Int,(x[d]+w.grid.extent[d]/2)/w.grid.spacing[d])+1,3)
        r=ceil(Int,radius/h)
        for dz in (w.grid.size[3]==1 ? (0:0) : (-r:r)),dy in -r:r,dx in -r:r
            I=CartesianIndex(index .+ (dx,dy,dz))
            all(d->1<=I[d]<=w.grid.size[d],1:3) || continue
            xyz=ntuple(d->(I[d]-.5)*w.grid.spacing[d]-w.grid.extent[d]/2,3)
            distance=sum((xyz[d]-x[d])^2 for d in 1:(w.grid.size[3]==1 ? 2 : 3))
            a[I]=max(a[I],Float32(value*exp(-distance/(2radius^2))))
        end
    end
    for e in eachindex(p.bond_ids)
        i,j=p.bonds[:,e];a=p.position[:,i];b=p.position[:,j]
        if p.damage[e]==1
            splat!(scar,a,1);splat!(scar,b,1);continue
        end
        strain=max(0.,norm(b-a)/p.rest_length[e]-1)/p.threshold[e]
        n=max(1,ceil(Int,norm(b-a)/(h/2)))
        for q in 0:n
            x=a.+(q/n).*(b-a);splat!(density,x,1);splat!(stress,x,strain)
        end
    end
    (density=.12f0.*density,activity=stress.+.5f0.*scar,structure=density.+scar,mask=ones(Float32,w.grid.size))
end
"""OBJ vertices and live line segments with persistent original node order."""
function export_bond_graph(p::BondNetwork,path)
    open(path,"w") do io
        println(io,"# Current deformed bond network; broken segments omitted")
        for i in eachindex(p.node_ids);println(io,"v ",join(p.position[:,i]," "));end
        for e in eachindex(p.bond_ids)
            p.damage[e]==1 && continue
            println(io,"l ",p.bonds[1,e]," ",p.bonds[2,e])
        end
    end
    path
end
