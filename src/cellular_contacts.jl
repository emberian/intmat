"""Persistent cell bonds. Endpoints are stable cell IDs, never column indices."""
mutable struct CellBonds
    endpoints::Vector{Tuple{Int,Int}}
    rest_lengths::Vector{Float64}
end
CellBonds() = CellBonds(Tuple{Int,Int}[], Float64[])

"""Add a bond (including a parent–daughter bond); return whether it was new."""
function add_cell_bond!(b::CellBonds, a::Integer, c::Integer, rest_length::Real)
    a != c || throw(ArgumentError("a cell cannot bond to itself"))
    isfinite(rest_length) && rest_length > 0 || throw(ArgumentError("positive finite bond rest length required"))
    edge = minmax(Int(a), Int(c))
    edge in b.endpoints && return false
    push!(b.endpoints, edge); push!(b.rest_lengths, Float64(rest_length))
    true
end

@inline _cellular_ndim(grid) = grid.size[3] == 1 ? 2 : 3
# Wrapped positions put the shifted displacement within one adjacent period.
# Sterbenz subtraction is exact on [L,2L); the negative branch uses the same
# r+L operation as Base.mod. Preserve the general path for unwrapped callers.
@inline function _cellular_mod_period(u::Float64,L::Float64)
    if 0<=u<L
        u==0 ? 0.0 : u
    elseif -L<u<0
        u+L
    elseif L<=u<2L
        u-L
    else
        mod(u,L)
    end
end

@inline function _cellular_delta(grid, x, i, j)
    @inbounds ntuple(d -> d > _cellular_ndim(grid) ? 0. :
        _cellular_mod_period(Float64(x[d,j])-Float64(x[d,i])+Float64(grid.extent[d])/2,
            Float64(grid.extent[d]))-Float64(grid.extent[d])/2, 3)
end
@inline _cellular_length(delta) = sqrt(sum(abs2, delta))

function _cellular_validate_geometry(grid, positions, radii)
    n = length(radii)
    size(positions) == (3,n) || throw(DimensionMismatch("cell positions must have shape (3,n)"))
    all(isfinite,positions) && all(r -> isfinite(r) && r > 0,radii) ||
        throw(ArgumentError("finite positions and positive finite radii required"))
    nothing
end

"""Return column-index pairs separated by at most `radii[i]+radii[j]+range`.

`range` is an additional physical surface gap. Coordinates are physical periodic
coordinates, with z ignored in 2D. Sparse bins have width at least the maximum
interaction distance; only adjacent bins are visited, including across seams.
"""
function cellular_pairs(grid, positions, radii; range::Real=0.)
    _cellular_validate_geometry(grid,positions,radii)
    isfinite(range) && range >= 0 || throw(ArgumentError("nonnegative finite range required"))
    n = length(radii); pairs = Tuple{Int,Int}[]
    n < 2 && return pairs
    cutoff = 2maximum(radii)+range; nd = _cellular_ndim(grid)
    # Capping bin counts by n prevents integer overflow for tiny physical cells.
    dims = ntuple(d -> d > nd ? 1 : max(1,Int(floor(min(n,grid.extent[d]/cutoff)))),3)
    # Dense linked bin heads avoid a hash lookup per neighboring bin and a
    # separately allocated vector per occupied bin. Bound the dense workspace
    # for tiny interaction radii; the sparse fallback uses the same traversal.
    dense=Float64(dims[1])*dims[2]*dims[3]<=max(4096,16n)
    heads=dense ? zeros(Int,prod(dims)) : Dict{Int,Int}()
    _cellular_pairs!(pairs,grid,positions,radii,range,dims,nd,heads)
end

@inline _cellular_bin_head(heads::Vector{Int},bin)=heads[bin]
@inline _cellular_bin_head(heads::Dict{Int,Int},bin)=get(heads,bin,0)
@inline _cellular_bin_linear(key,dims)=key[1]+dims[1]*((key[2]-1)+dims[2]*(key[3]-1))
@inline function _cellular_bin_neighbor(key,offset,n)
    value=key+offset
    value<1 ? value+n : value>n ? value-n : value
end

function _cellular_pairs!(pairs,grid,positions,radii,range,dims,nd,heads)
    n=length(radii); next=zeros(Int,n)
    keys = Vector{NTuple{3,Int}}(undef,n)
    for i in 1:n
        key = ntuple(d -> d > nd ? 1 : 1+min(dims[d]-1,
            floor(Int,mod(Float64(positions[d,i]),Float64(grid.extent[d]))/grid.extent[d]*dims[d])),3)
        keys[i]=key; bin=_cellular_bin_linear(key,dims)
        next[i]=_cellular_bin_head(heads,bin); heads[bin]=i
    end
    # At most three unique offsets per axis; compute once, not per cell.
    offsets=ntuple(d -> dims[d]==1 ? (0:0) : dims[d]==2 ? (0:1) : (-1:1),3)
    for i in 1:n
        key = keys[i]
        for dz in offsets[3], dy in offsets[2], dx in offsets[1]
            neighbor=(_cellular_bin_neighbor(key[1],dx,dims[1]),
                _cellular_bin_neighbor(key[2],dy,dims[2]),
                _cellular_bin_neighbor(key[3],dz,dims[3]))
            j=_cellular_bin_head(heads,_cellular_bin_linear(neighbor,dims))
            while j!=0
                if j>i
                    _cellular_length(_cellular_delta(grid,positions,i,j)) <= radii[i]+radii[j]+range && push!(pairs,(i,j))
                end
                j=next[j]
            end
        end
    end
    sort!(pairs)
end

function _cellular_validate_bonds(b)
    length(b.endpoints) == length(b.rest_lengths) || throw(DimensionMismatch("bond arrays differ"))
    all(e -> e[1] < e[2],b.endpoints) && length(unique(b.endpoints)) == length(b.endpoints) ||
        throw(ArgumentError("bond endpoints must be unique canonical ID pairs"))
    all(l -> isfinite(l) && l > 0,b.rest_lengths) || throw(ArgumentError("invalid bond rest length"))
end

@inline function _cellular_energy_forces(grid, ids, x, radii, bonds, index, repulsion, stiffness; forces=false,pairs=nothing)
    _cellular_energy_forces(grid,ids,x,radii,bonds,index,repulsion,stiffness,Val(forces),pairs)
end

function _cellular_energy_forces(grid, ids, x, radii, bonds, index, repulsion, stiffness, ::Val{F}, pairs) where F
    f=F ? Matrix{Float64}(undef,size(x)) : nothing
    rows=F ? Vector{Float64}(undef,length(ids)) : Float64[]
    _cellular_energy_forces!(grid,ids,x,radii,bonds,index,repulsion,stiffness,Val(F),pairs,f,rows)
end

# Buffers live for one mechanics call; they are not persistent simulation state.
function _cellular_energy_forces!(grid,ids,x,radii,bonds,index,repulsion,stiffness,::Val{F},pairs,f,rows) where F
    forces=F
    F && (fill!(f,0.);fill!(rows,0.))
    energy=0.;overlaps=0
    candidates=pairs===nothing ? cellular_pairs(grid,x,radii) : pairs
    @inbounds for (i,j) in candidates
        delta = _cellular_delta(grid,x,i,j); distance = _cellular_length(delta)
        overlap = max(0.,radii[i]+radii[j]-distance)
        overlap > 0 || continue
        overlaps += 1; energy += .5repulsion*overlap^2
        if forces
            rows[i] += repulsion; rows[j] += repulsion
            # A reproducible subgradient separates exactly coincident cells.
            direction = distance > 0 ? delta ./ distance : (ids[i] < ids[j] ? (1.,0.,0.) : (-1.,0.,0.))
            for d in 1:3
                force = -repulsion*overlap*direction[d]
                f[d,i] += force; f[d,j] -= force
            end
        end
    end
    @inbounds for e in eachindex(bonds.endpoints)
        a,b = bonds.endpoints[e]; i = get(index,a,0); j = get(index,b,0)
        (i == 0 || j == 0) && continue
        delta = _cellular_delta(grid,x,i,j); distance = _cellular_length(delta)
        extension = max(0.,distance-bonds.rest_lengths[e])
        energy += .5stiffness*extension^2
        if forces && extension > 0
            rows[i] += stiffness; rows[j] += stiffness
            for d in 1:3
                force = stiffness*extension*delta[d]/distance
                f[d,i] += force; f[d,j] -= force
            end
        end
    end
    (;energy,forces=f,rows,overlaps)
end

"""Overlap penalty plus tensile spring energy; no kinetic or chemical energy."""
function cellular_mechanical_energy(grid,ids,positions,radii,bonds::CellBonds;repulsion=10.,stiffness=2.)
    _cellular_validate_geometry(grid,positions,radii); _cellular_validate_bonds(bonds)
    length(ids)==length(radii) && length(unique(ids))==length(ids) || throw(ArgumentError("unique cell IDs required"))
    all(k -> isfinite(k) && k >= 0,(repulsion,stiffness)) || throw(ArgumentError("invalid stiffness"))
    index=Dict(id=>i for (i,id) in enumerate(ids))
    _cellular_energy_forces(grid,ids,positions,radii,bonds,index,repulsion,stiffness).energy
end

function _cellular_refresh_bonds!(grid,ids,x,radii,bonds,index,adhesion,threshold,break_strain,max_degree,bond_range,stiffness; pairs=nothing)
    keep = trues(length(bonds.endpoints)); broken = 0; loss = 0.
    degree = zeros(Int,length(ids))
    for e in eachindex(bonds.endpoints)
        a,b = bonds.endpoints[e]; i=get(index,a,0); j=get(index,b,0)
        alive = i != 0 && j != 0
        extension = alive ? max(0.,_cellular_length(_cellular_delta(grid,x,i,j))-bonds.rest_lengths[e]) : 0.
        if !alive || adhesion[i] < threshold || adhesion[j] < threshold || extension > break_strain*bonds.rest_lengths[e]
            keep[e]=false; broken+=1; loss+=.5stiffness*extension^2
        else
            degree[i]+=1; degree[j]+=1
        end
    end
    bonds.endpoints = bonds.endpoints[keep]; bonds.rest_lengths = bonds.rest_lengths[keep]
    existing = Set(bonds.endpoints); formed = 0
    max_degree==0 && return (;broken,formed,loss)
    candidates = pairs===nothing ? cellular_pairs(grid,x,radii;range=bond_range) :
        [(i,j) for (i,j) in pairs if _cellular_length(_cellular_delta(grid,x,i,j))<=radii[i]+radii[j]+bond_range]
    sort!(candidates;by=ij -> minmax(ids[ij[1]],ids[ij[2]]))
    for (i,j) in candidates
        edge=minmax(ids[i],ids[j])
        (edge in existing || degree[i]>=max_degree || degree[j]>=max_degree || adhesion[i]<threshold || adhesion[j]<threshold) && continue
        rest=max(Float64(radii[i]+radii[j]),_cellular_length(_cellular_delta(grid,x,i,j)))
        add_cell_bond!(bonds,ids[i],ids[j],rest)
        push!(existing,edge); degree[i]+=1;degree[j]+=1;formed+=1
    end
    (;broken,formed,loss)
end

"""Advance physical periodic cell bodies using overdamped, unit-mobility mechanics.

Repulsion is a soft overlap potential; bonds resist tension only. Both forces are
pairwise equal/opposite. Adjacent cells with mutual adhesion at least
`adhesion_threshold` acquire zero-energy bonds (rest length at least their radii
sum). Bonds survive reordering and disappear with an endpoint, withdrawal, or
excess tensile strain. No ancestry gate is imposed.

`motility` is a desired physical velocity. With an `energy` vector, accepted
motor moves debit positive potential increase plus `motor_cost*|dx|²/h` drag.
Potential cost is shared by squared displacement; each cell's motion is reduced
until its share is affordable. Downhill motion never credits chemical energy.
Without `energy`, funding belongs to the caller; `spent_motor_work` reports the
model's work requirement. This discrete budget is not a calibrated motor model.
Adaptive substeps bound displacement and stiffness, and passive relaxation is
backtracked until potential does not increase. Displacements are never clipped.
`relaxation_energy_loss` is the actual potential decrease of this discrete
mechanical model, not a claim of exact chemical/thermal conservation. Soft
contacts can retain a small overlap under sustained forcing or crowding.
"""
# The candidate skin allows each endpoint to move by 2travel_limit. Reanchor
# already at travel_limit, leaving a full factor of two before completeness
# could fail. This test only controls work reuse, never force arithmetic.
function _cellular_inside_pair_skin(grid,x,origin,limit)
    limit2=limit^2;nd=_cellular_ndim(grid)
    @inbounds for i in axes(x,2)
        distance2=0.
        for d in 1:nd
            L=Float64(grid.extent[d])
            delta=_cellular_mod_period(x[d,i]-origin[d,i]+L/2,L)-L/2
            distance2+=delta^2
        end
        distance2<=limit2 || return false
    end
    true
end
function _cellular_reanchor_pairs!(pairs,origin,grid,x,radii,range)
    fresh=cellular_pairs(grid,x,radii;range)
    empty!(pairs);append!(pairs,fresh);copyto!(origin,x)
    nothing
end

function cellular_mechanics!(grid,ids,positions,radii,bonds::CellBonds,dt::Real;
    adhesion=ones(length(ids)),motility=zeros(3,length(ids)),repulsion=10.,stiffness=2.,
    break_strain=.75,max_degree::Integer=6,adhesion_threshold=.5,bond_range=0.,max_substeps::Integer=10000,
    energy=nothing,motor_cost=1.,reuse_candidates::Bool=true)
    _cellular_validate_geometry(grid,positions,radii); _cellular_validate_bonds(bonds)
    n=length(ids); nd=_cellular_ndim(grid)
    length(radii)==length(adhesion)==n && size(motility)==(3,n) || throw(DimensionMismatch("cell action dimensions differ"))
    length(unique(ids))==n || throw(ArgumentError("unique cell IDs required"))
    all(a->isfinite(a)&&0<=a<=1,adhesion) && all(isfinite,motility) || throw(ArgumentError("invalid adhesion or motility"))
    isfinite(dt) && dt>=0 && all(k->isfinite(k)&&k>=0,(repulsion,stiffness,break_strain,bond_range)) &&
        0<adhesion_threshold<=1 && max_degree>=0 && max_substeps>0 || throw(ArgumentError("invalid mechanics parameters"))
    isfinite(motor_cost) && motor_cost>0 || throw(ArgumentError("positive motor drag coefficient required"))
    if energy !== nothing
        length(energy)==n && all(e->isfinite(e)&&e>=0,energy) || throw(ArgumentError("invalid motor energy budget"))
    end
    spent_motor_work=zeros(n);motor_drag_loss=0.;active_distance=zeros(n)
    index=Dict(id=>i for (i,id) in enumerate(ids))
    x=Float64.(positions)
    for i in 1:n,d in 1:3; x[d,i]=d>nd ? 0. : mod(x[d,i],Float64(grid.extent[d])); end
    initial_pairs=cellular_pairs(grid,x,radii;range=bond_range)
    potential_before=_cellular_energy_forces(grid,ids,x,radii,bonds,index,repulsion,stiffness,Val(false),initial_pairs).energy
    broken=0;formed=0;topology_energy_loss=0.;motor_potential_work=0.;relaxation_energy_loss=0.;substeps=0
    refresh=_cellular_refresh_bonds!(grid,ids,x,radii,bonds,index,adhesion,adhesion_threshold,break_strain,max_degree,bond_range,stiffness;pairs=initial_pairs)
    broken+=refresh.broken;formed+=refresh.formed;topology_energy_loss+=refresh.loss
    remaining=Float64(dt); travel_limit=n==0 ? Inf : .2minimum(radii)
    motors=Float64.(motility); nd==2 && (motors[3,:].=0.)
    # Reuse trial state across adaptive substeps and rejected budget attempts.
    budget_motors=similar(motors);velocity=similar(motors)
    motor_x=similar(x);trial=similar(x);motor_spend=zeros(n);weights=zeros(n)
    state_forces=similar(x);state_rows=zeros(n)
    motor_forces=similar(x);motor_rows=zeros(n)
    candidates=Tuple{Int,Int}[];pair_origin=copy(x);pairs_ready=false
    while remaining>0 && n>0
        substeps < max_substeps || throw(ArgumentError("mechanics exceeded max_substeps; reduce dt or forcing"))
        # Every accepted motor and passive displacement is <= travel_limit.
        # A pair can close by at most 4travel_limit across both split moves.
        # This physical skin reuses one broadphase for all trial evaluations.
        if !reuse_candidates || !pairs_ready || !_cellular_inside_pair_skin(grid,x,pair_origin,travel_limit)
            _cellular_reanchor_pairs!(candidates,pair_origin,grid,x,radii,4travel_limit+bond_range)
            pairs_ready=true
        end
        state=_cellular_energy_forces!(grid,ids,x,radii,bonds,index,repulsion,stiffness,Val(true),candidates,state_forces,state_rows)
        copyto!(budget_motors,motors)
        if energy !== nothing
            for i in 1:n
                drag=motor_cost*remaining*sum(abs2,view(budget_motors,:,i))
                drag>energy[i] && (budget_motors[:,i].*=sqrt(energy[i]/drag))
            end
        end
        # Snapshot loop-carried arrays before closure capture. Capturing the
        # reassigned bindings boxes them and makes the numerical loop dynamic.
        speed=let budget_motors=budget_motors, state=state
            maximum(sqrt(sum(abs2,view(budget_motors,:,i)))+sqrt(sum(abs2,view(state.forces,:,i))) for i in 1:n)
        end
        rate=maximum(state.rows)
        h=min(remaining,speed>0 ? travel_limit/speed : Inf,rate>0 ? .45/rate : Inf)
        accepted=false
        for attempt in 1:60
            copyto!(velocity,budget_motors)
            if energy !== nothing
                for i in 1:n
                    drag=motor_cost*h*sum(abs2,view(velocity,:,i))
                    drag>energy[i] && (velocity[:,i].*=sqrt(energy[i]/drag))
                end
            end
            motor_state=state
            funded=false
            for budget_attempt in 1:80
                @. motor_x=x+h*velocity
                for i in 1:n,d in 1:nd;motor_x[d,i]=mod(motor_x[d,i],Float64(grid.extent[d]));end
                if reuse_candidates && !_cellular_inside_pair_skin(grid,motor_x,pair_origin,travel_limit)
                    _cellular_reanchor_pairs!(candidates,pair_origin,grid,x,radii,4travel_limit+bond_range)
                end
                motor_state=_cellular_energy_forces!(grid,ids,motor_x,radii,bonds,index,repulsion,stiffness,Val(true),candidates,motor_forces,motor_rows)
                for i in 1:n
                    weights[i]=sum(abs2,view(velocity,:,i))
                end
                weight_sum=sum(weights)
                potential_cost=max(0.,motor_state.energy-state.energy)
                for i in 1:n
                    motor_spend[i]=motor_cost*h*weights[i]+(weight_sum>0 ? potential_cost*weights[i]/weight_sum : 0.)
                end
                if energy === nothing || all(i->motor_spend[i]<=energy[i],eachindex(motor_spend))
                    funded=true;break
                end
                for i in 1:n
                    motor_spend[i]>energy[i] && (velocity[:,i].*=.9sqrt(energy[i]/motor_spend[i]))
                end
            end
            funded || throw(ArgumentError("motor budget could not resolve displacement"))
            @. trial=motor_x+h*motor_state.forces
            for i in 1:n,d in 1:nd;trial[d,i]=mod(trial[d,i],Float64(grid.extent[d]));end
            passive_travel=let motor_state=motor_state
                h*maximum(sqrt(sum(abs2,view(motor_state.forces,:,i))) for i in 1:n)
            end
            # Reject out-of-skin moves before evaluating the cached pair list.
            trial_energy=Inf
            if passive_travel<=travel_limit
                if reuse_candidates && !_cellular_inside_pair_skin(grid,trial,pair_origin,travel_limit)
                    _cellular_reanchor_pairs!(candidates,pair_origin,grid,x,radii,4travel_limit+bond_range)
                end
                trial_energy=_cellular_energy_forces(grid,ids,trial,radii,bonds,index,repulsion,stiffness,Val(false),candidates).energy
            end
            if trial_energy<=motor_state.energy && passive_travel<=travel_limit
                motor_potential_work+=motor_state.energy-state.energy
                relaxation_energy_loss+=motor_state.energy-trial_energy
                for i in 1:n
                    spent_motor_work[i]+=motor_spend[i]
                    active_distance[i]+=h*sqrt(sum(abs2,view(velocity,:,i)))
                    energy !== nothing && (energy[i]-=motor_spend[i])
                end
                motor_drag_loss+=motor_cost*h*sum(abs2,velocity)
                copyto!(x,trial); accepted=true;break
            end
            h/=2
        end
        accepted && h>0 && remaining-h<remaining || throw(ArgumentError("mechanical relaxation could not resolve timestep"))
        remaining=max(0.,remaining-h);substeps+=1
        refresh=_cellular_refresh_bonds!(grid,ids,x,radii,bonds,index,adhesion,adhesion_threshold,break_strain,max_degree,bond_range,stiffness;pairs=candidates)
        broken+=refresh.broken;formed+=refresh.formed;topology_energy_loss+=refresh.loss
    end
    positions .= x
    potential_after=_cellular_energy_forces(grid,ids,positions,radii,bonds,index,repulsion,stiffness).energy
    # Includes final conversion to the caller's position precision.
    energy_residual=potential_after-potential_before-motor_potential_work+relaxation_energy_loss+topology_energy_loss
    (;potential_before,potential_after,motor_potential_work,relaxation_energy_loss,topology_energy_loss,
      energy_residual,formed,broken,bonds=length(bonds.endpoints),substeps,
      spent_motor_work,motor_drag_loss,active_distance)
end

"""Conservative, positive exchange of intracellular species amounts along bonds.

Each edge applies exact two-compartment concentration equilibration with passive
conductance. Volumes are positive cell volumes (default one); columns are stable
IDs in any order. Canonical edge order makes results independent of array order.
The edge splitting is first-order for networks, not a simultaneous exact solve.
"""
function exchange_cell_metabolites!(amounts,ids,bonds::CellBonds,dt::Real;rate=0.,volumes=ones(length(ids)),openness=ones(length(ids)))
    eltype(amounts)<:AbstractFloat || throw(ArgumentError("metabolite amounts require floating-point storage"))
    size(amounts,2)==length(ids)==length(volumes) || throw(DimensionMismatch("metabolite columns differ"))
    length(unique(ids))==length(ids) || throw(ArgumentError("unique cell IDs required"))
    _cellular_validate_bonds(bonds)
    all(a->isfinite(a)&&a>=0,amounts) && all(v->isfinite(v)&&v>0,volumes) || throw(ArgumentError("invalid metabolite amounts or volumes"))
    length(openness)==length(ids) && all(o->isfinite(o)&&0<=o<=1,openness) || throw(ArgumentError("invalid exchange openness"))
    isfinite(dt)&&dt>=0&&isfinite(rate)&&rate>=0 || throw(ArgumentError("invalid exchange rate or dt"))
    index=Dict(id=>i for (i,id) in enumerate(ids));transferred=0.
    for e in sortperm(bonds.endpoints)
        a,b=bonds.endpoints[e];i=get(index,a,0);j=get(index,b,0)
        (i==0||j==0) && continue
        vi=Float64(volumes[i]);vj=Float64(volumes[j])
        fraction=-expm1(-Float64(rate)*Float64(dt)*openness[i]*openness[j]*(inv(vi)+inv(vj)))
        fraction==0 && continue
        for s in axes(amounts,1)
            ai=Float64(amounts[s,i]);aj=Float64(amounts[s,j]);total=ai+aj
            transfer=fraction*(ai*vj-aj*vi)/(vi+vj)
            transfer==0 && continue
            first=eltype(amounts)(clamp(ai-transfer,0.,total))
            # A low-precision rounded first compartment must not exceed the
            # pair's total, or its complementary compartment becomes negative.
            Float64(first)>total && (first=prevfloat(first))
            second=eltype(amounts)(max(0.,total-Float64(first)))
            amounts[s,i]=first;amounts[s,j]=second
            transferred+=abs(ai-Float64(first))
        end
    end
    (;transferred)
end

exchange_cell_metabolites!(amounts,bonds::CellBonds,ids,volumes,dt::Real;conductance=0.,openness=ones(length(ids))) =
    exchange_cell_metabolites!(amounts,ids,bonds,dt;rate=conductance,volumes,openness)
