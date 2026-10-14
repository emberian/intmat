"""Fixed recurrent particle policy. Rates are nondimensional; no online learning."""
struct ParticleRule{T}
    recurrence::Matrix{T}
    sensing::Matrix{T}
    motion::Matrix{T}
    secretion::Vector{T}
    uptake_rate::T
    secretion_rate::T
    speed::T
    relaxation::T
end

function ParticleRule(; T=Float32, hidden=4, uptake_rate=0.3,
                      secretion_rate=0.12, speed=0.5, relaxation=1.0)
    hidden >= 2 || throw(ArgumentError("particle hidden size must be at least two"))
    all(x -> isfinite(x) && x >= 0, (uptake_rate, secretion_rate, speed, relaxation)) ||
        throw(ArgumentError("particle rates must be finite and nonnegative"))
    w = zeros(T, hidden, hidden)
    for i in 1:hidden
        w[i,i] = T(0.9)
        w[i,mod1(i+1,hidden)] = T(isodd(i) ? -1.4 : 1.4)
    end
    sense = zeros(T, hidden, 3)
    sense[:,1] .= T(0.25)
    sense[:,2] .= T(-0.3)
    sense[:,3] .= T(0.15)
    motion = zeros(T, 3, hidden)
    for d in 1:3
        motion[d,mod1(d,hidden)] = one(T)
    end
    ParticleRule(w,sense,motion,fill(T(1/hidden),hidden),T(uptake_rate),
                 T(secretion_rate),T(speed),T(relaxation))
end

"""CPU moving carrier: columns are particles, positions are periodic grid-cell coordinates.

IDs, hidden state and storage move together under `reorder_particles!`. Exchange
amounts use the same per-cell inventory units as the grid arrays. Hidden states
are private to this rule and are never mixed with grid latent channels.
"""
mutable struct ParticleCarrier{T,R}
    ids::Vector{Int}
    positions::Matrix{T}
    hidden::Matrix{T}
    storage::Vector{T}
    rule::R
end

function ParticleCarrier(positions::AbstractMatrix; rule=ParticleRule(T=eltype(positions)),
                         ids=collect(1:size(positions,2)), hidden=nothing, storage=nothing)
    T = eltype(rule.recurrence)
    size(positions,1) == 3 || throw(DimensionMismatch("positions must have three rows"))
    n = size(positions,2)
    length(ids) == n && length(unique(ids)) == n || throw(ArgumentError("IDs must be unique"))
    h = hidden === nothing ? zeros(T,size(rule.recurrence,1),n) : Matrix{T}(hidden)
    s = storage === nothing ? zeros(T,n) : Vector{T}(storage)
    size(h) == (size(rule.recurrence,1),n) || throw(DimensionMismatch("hidden shape"))
    length(s) == n || throw(DimensionMismatch("storage shape"))
    all(isfinite,positions) && all(isfinite,h) && all(x->isfinite(x)&&x>=0,s) ||
        throw(ArgumentError("particle state must be finite, with nonnegative storage"))
    ParticleCarrier(Int.(ids),Matrix{T}(positions),h,s,rule)
end

function reorder_particles!(p::ParticleCarrier, permutation)
    sort(collect(permutation)) == collect(eachindex(p.ids)) || throw(ArgumentError("invalid permutation"))
    p.ids = p.ids[permutation]
    p.positions = p.positions[:,permutation]
    p.hidden = p.hidden[:,permutation]
    p.storage = p.storage[permutation]
    p
end

total_particle_nutrient(p::ParticleCarrier) = sum(p.storage)

"""Conservative cell-local uptake, recurrence, deposition, then motion.

Particles sharing a cell receive a proportional share of available nutrient.
Buckets are built and reduced in persistent-ID order, so spatial reordering does
not change the trajectory. This method deliberately requires CPU `Array`s:
GPU scatter/exchange is not implemented. Motion speed is physical distance/time;
`spacing` converts it to cell coordinates. Optional grid velocity is sampled in
the occupied cell. Runtime excluding grid validation is O(particles log particles).
"""
function particle_step!(p::ParticleCarrier{T}, nutrient::Array{T,3},
                        deposit::Array{T,3}, dt::Real;
                        spacing=(1,1,1), velocity=nothing) where T
    size(nutrient) == size(deposit) || throw(DimensionMismatch("exchange grids differ"))
    isfinite(dt) && dt >= 0 || throw(ArgumentError("dt must be finite and nonnegative"))
    length(spacing) == 3 && all(x->isfinite(x)&&x>0,spacing) ||
        throw(ArgumentError("spacing must contain three positive finite values"))
    velocity === nothing || (length(velocity) == 3 &&
        all(v->v isa Array{T,3} && size(v)==size(nutrient) && all(isfinite,v),velocity)) ||
        throw(ArgumentError("velocity must contain three matching finite CPU grids"))
    all(x->isfinite(x)&&x>=0,nutrient) && all(x->isfinite(x)&&x>=0,deposit) ||
        throw(ArgumentError("exchange inventories must be finite and nonnegative"))
    dt == 0 && return p
    dims = size(nutrient)
    buckets = Dict{Int,Vector{Int}}()
    particle_cells = zeros(Int,length(p.ids))
    indices = LinearIndices(nutrient)
    for i in sortperm(p.ids)
        for d in 1:3
            p.positions[d,i] = dims[d] == 1 ? zero(T) : mod(p.positions[d,i],T(dims[d]))
        end
        cell = ntuple(d->floor(Int,p.positions[d,i])+1,3)
        index = indices[cell...]
        particle_cells[i] = index
        push!(get!(()->Int[],buckets,index),i)
    end
    r = p.rule
    requests = zeros(T,length(p.ids))
    alpha = T(-expm1(-Float64(r.relaxation)*dt))
    for cell in sort!(collect(keys(buckets)))
        bucket = buckets[cell]
        isempty(bucket) && continue
        requested = zero(T)
        for i in bucket
            observation = T[nutrient[cell],deposit[cell],p.storage[i]]
            target = tanh.(r.recurrence*p.hidden[:,i] + r.sensing*observation)
            p.hidden[:,i] .= (one(T)-alpha).*p.hidden[:,i] .+ alpha.*target
            requests[i] = T(dt)*r.uptake_rate*(T(0.5)+T(0.5)*tanh(p.hidden[1,i]))
            requested += requests[i]
        end
        available = nutrient[cell]
        scale = requested > 0 ? min(one(T),available/requested) : zero(T)
        # Sequential residual cap avoids negative inventory from floating-point sums.
        remaining = available
        deposited = zero(T)
        for i in bucket
            uptake = min(remaining,requests[i]*scale)
            remaining -= uptake
            p.storage[i] += uptake
            gate = T(0.5)+T(0.5)*tanh(sum(r.secretion .* p.hidden[:,i]))
            secretion = min(p.storage[i],T(dt)*r.secretion_rate*gate*p.storage[i])
            p.storage[i] -= secretion
            deposited += secretion
        end
        nutrient[cell] = remaining
        deposit[cell] += deposited
    end
    for i in eachindex(p.ids)
        active_velocity = r.motion*p.hidden[:,i]
        for d in 1:3
            drift = velocity === nothing ? zero(T) : velocity[d][particle_cells[i]]
            p.positions[d,i] = dims[d] == 1 ? zero(T) :
                mod(p.positions[d,i]+T(dt)*(r.speed*tanh(active_velocity[d])+drift)/T(spacing[d]),T(dims[d]))
        end
    end
    p
end

advance!(w::World, p::ParticleCarrier, dt) =
    particle_step!(p, w.fields.nutrient, w.fields.deposit, dt;
        spacing=w.grid.spacing, velocity=(w.fields.vx,w.fields.vy,w.fields.vz))

function validate(w::World,p::ParticleCarrier{T}) where T
    f=w.fields
    all(a->a isa Array{T,3},(f.nutrient,f.deposit,f.vx,f.vy,f.vz)) ||
        throw(ArgumentError("particle exchange requires matching CPU arrays"))
    r=p.rule; n=length(p.ids); h=size(r.recurrence,1)
    h>=2 && size(r.recurrence)==(h,h) && size(r.sensing)==(h,3) &&
        size(r.motion)==(3,h) && length(r.secretion)==h ||
        throw(DimensionMismatch("particle rule dimensions"))
    size(p.positions)==(3,n) && size(p.hidden)==(h,n) && length(p.storage)==n ||
        throw(DimensionMismatch("particle state dimensions"))
    length(unique(p.ids))==n || throw(ArgumentError("particle IDs must be unique"))
    all(a->all(isfinite,a),(p.positions,p.hidden,r.recurrence,r.sensing,r.motion,r.secretion)) &&
        all(x->isfinite(x)&&x>=0,p.storage) &&
        all(x->isfinite(x)&&x>=0,(r.uptake_rate,r.secretion_rate,r.speed,r.relaxation)) ||
        throw(ArgumentError("particle state and rule must be finite, with nonnegative inventories and rates"))
    all(a->size(a)==w.grid.size,(f.nutrient,f.deposit,f.vx,f.vy,f.vz)) ||
        throw(DimensionMismatch("particle world grid dimensions"))
    all(a->all(isfinite,a),(f.vx,f.vy,f.vz)) &&
        all(a->all(x->isfinite(x)&&x>=0,a),(f.nutrient,f.deposit)) ||
        throw(ArgumentError("particle world requires finite velocity and nonnegative finite inventories"))
    nothing
end
