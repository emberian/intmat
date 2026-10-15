# Checkpoint payloads contain only host arrays and immutable metadata, suitable
# for Julia Serialization in the execution environment. Never deserialize untrusted files.
const DEVICE_CELLULAR_SOURCE_PIN=(
    model=bytes2hex(SHA.sha256(read(joinpath(@__DIR__,"device_cellular.jl")))),
    initialization=bytes2hex(SHA.sha256(read(joinpath(@__DIR__,"device_cellular_initialization.jl")))),
    founder_initialization=bytes2hex(SHA.sha256(read(joinpath(@__DIR__,"cellular_initialization.jl")))),
    support=bytes2hex(SHA.sha256(read(joinpath(@__DIR__,"device_cellular_support.jl")))),
    kernels=bytes2hex(SHA.sha256(read(joinpath(@__DIR__,"device_cellular_kernels.jl")))),
    checkpoint=bytes2hex(SHA.sha256(read(@__FILE__))),
    protocol=isfile(joinpath(@__DIR__,"device_cellular_protocol.jl")) ? bytes2hex(SHA.sha256(read(joinpath(@__DIR__,"device_cellular_protocol.jl")))) : nothing)
"""Synchronize and capture a complete source-pinned host checkpoint.

No checkpoint is emitted from a partial/failed step. Scratch arrays, contact
bins, free lists, counters, identities and compensated ledgers are included.
"""
function device_cellular_checkpoint(w::World,p::DeviceCellular)
    device_cellular_validate(w,p)
    snapshot=device_cellular_snapshot(p)
    snapshot.attempted_steps==snapshot.environment.metadata[6] || throw(ArgumentError("checkpoint requires a complete step boundary"))
    (format="Interior.DeviceCellular",version=2,julia_version=VERSION,source=DEVICE_CELLULAR_SOURCE_PIN,
        world_time=w.time,world_seed=w.seed,controller_contract=(sensors=CELLULAR_SENSORS,actions=CELLULAR_ACTIONS),snapshot=snapshot)
end
"""Restore a complete checkpoint onto the selected World backend.

This restores the declared experimental law on the exact source pin. The World
must have the same grid and zero unrelated nutrient. It does not reinterpret a
CellularEcology snapshot or resume a failed partial step.
"""
function device_cellular_restore(w::World,checkpoint)
    checkpoint.format=="Interior.DeviceCellular" && checkpoint.version==2 || throw(ArgumentError("unsupported device checkpoint"))
    checkpoint.julia_version==VERSION || throw(ArgumentError("checkpoint Julia version differs"))
    checkpoint.source==DEVICE_CELLULAR_SOURCE_PIN || throw(ArgumentError("checkpoint source pin differs from this runtime"))
    checkpoint.controller_contract==(sensors=CELLULAR_SENSORS,actions=CELLULAR_ACTIONS) || throw(ArgumentError("checkpoint controller layout differs"))
    x=checkpoint.snapshot;g=x.geometry;q=x.parameters
    w.grid.size==Tuple(Int.(g.dims)) && w.grid.extent==g.extent || throw(ArgumentError("checkpoint grid mismatch"))
    g.spacing==w.grid.spacing && g.nd==(w.grid.size[3]==1 ? 2 : 3) &&
        g.volume==prod(w.grid.spacing[d] for d in 1:Int(g.nd)) || throw(ArgumentError("checkpoint metric differs from grid"))
    all(x->isfinite(x)&&x>=0,Iterators.flatten(v isa Tuple ? v : (v,) for v in values(q))) &&
        all(>(0),(q.exchange_radius,q.density,q.store_capacity,q.energy_capacity,q.half_saturation,q.response_time,q.adaptation_time,q.motor_cost,q.minimum_mass,q.division_mass,q.radius_bound,q.bin_width,q.mechanical_substeps)) &&
        0<=q.efficiency<=1 && 0<=q.mutation_probability<=1 && q.minimum_mass<q.division_mass/2 || throw(ArgumentError("invalid checkpoint parameters"))
    (!hasproperty(w.fields,:nutrient)||all(iszero,Array(w.fields.nutrient))) || throw(ArgumentError("zero unrelated world nutrient before restore"))
    only(x.environment.status)==0 && 0<=x.attempted_steps==x.environment.metadata[6] || throw(ArgumentError("checkpoint is not a complete state"))
    all(a->!(eltype(a)<:AbstractFloat)||all(isfinite,a),values(x.state)) &&
        all(a->!(eltype(a)<:AbstractFloat)||all(isfinite,a),values(x.environment)) || throw(ArgumentError("nonfinite checkpoint"))
    C=Int(g.capacity);h=Int(g.hidden);K=Int(g.degree);S=prod(Int.(g.candidates));V=prod(Int.(g.dims))
    C>0 && h>0 && K>=0 && S>0 && V>0 || throw(ArgumentError("invalid checkpoint capacities"))
    all(a->a==0||a==1,x.state.active) || throw(ArgumentError("invalid checkpoint active flags"))
    g.support_strategy in (0,1) && g.support_radius2==q.exchange_radius^2 || throw(ArgumentError("invalid checkpoint support strategy"))
    expected_candidates=ntuple(d->d>g.nd ? 1 : min(w.grid.size[d],2ceil(Int,q.exchange_radius/w.grid.spacing[d])+2),3)
    expected_bins=ntuple(d->d>g.nd ? 1 : max(1,floor(Int,w.grid.extent[d]/q.bin_width)),3)
    Tuple(Int.(g.candidates))==expected_candidates && Tuple(Int.(g.bins))==expected_bins && K==q.max_degree || throw(ArgumentError("checkpoint support/bin capacities differ from law"))
    all(>=(0),x.environment.metadata[2:6]) || throw(ArgumentError("negative checkpoint counters"))
    all(k->Float64(x.environment.resource_ledger[k,1])+Float64(x.environment.resource_ledger_low[k,1])>=0 &&
        Float64(x.environment.resource_ledger[k,2])+Float64(x.environment.resource_ledger_low[k,2])>=0,1:3) || throw(ArgumentError("negative checkpoint budgets"))
    all(k->all(column->Float32(Float64(x.environment.resource_ledger[k,column])+Float64(x.environment.resource_ledger_low[k,column]))==x.environment.resource_ledger[k,column],1:2),1:3) || throw(ArgumentError("checkpoint remaining budgets are not normalized compensated pairs"))
    active=findall(==(1),x.state.active)
    length(active)==x.environment.metadata[2] && length(unique(x.state.ids[active]))==length(active) &&
        all(id->0<id<x.environment.metadata[1],x.state.ids[active]) || throw(ArgumentError("invalid checkpoint identities"))
    state,env=device_cellular_arrays(w.backend,C,h,K,S,V,Int.(g.bins);
        support_strategy=g.support_strategy==1 ? :axis : :materialized,axis_capacity=maximum(Int.(g.candidates)))
    for (target,source) in ((state,x.state),(env,x.environment))
        keys(target)==keys(source) || throw(ArgumentError("checkpoint array roles differ"))
        for key in keys(target)
            dst=getproperty(target,key);src=getproperty(source,key)
            size(dst)==size(src) && eltype(dst)==eltype(src) || throw(ArgumentError("checkpoint array shape/type differs: $key"))
            copyto!(dst,src)
        end
    end
    p=DeviceCellular(state,env,w.backend,q,g,x.attempted_steps,device_cellular_contract(state,env,q,g))
    device_cellular_launch!(p,dc_validate_complete!,state,env,g;ndrange=max(C,V))
    device_cellular_status!(p)
    w.time=checkpoint.world_time;w.seed=checkpoint.world_seed
    p
end
