# Opt-in ensemble composition. No numerical kernel or biological law is copied.
using Interior, SHA
isdefined(Main,:CellularWorldBatch) || include("cellular_world_batches.jl")
isdefined(Interior,:DeviceCatalyticHabitat) || error("Device habitat implementation required")

struct ConstructiveBatchCells{B,R}
    base::B
    habitat_seeds::R
end
@inline Base.getproperty(p::ConstructiveBatchCells,k::Symbol)=k in (:base,:habitat_seeds) ? getfield(p,k) : getproperty(getfield(p,:base),k)
Base.setproperty!(p::ConstructiveBatchCells,k::Symbol,v)=setproperty!(getfield(p,:base),k,v)
Interior.device_cellular_validate(w::World,p::ConstructiveBatchCells)=Interior.device_cellular_validate(w,p.base)
Interior.device_cellular_status!(p::ConstructiveBatchCells)=(Interior.device_cellular_status!(p.base);p)
Interior.device_cellular_timestep_limit(p::ConstructiveBatchCells)=Interior.device_cellular_timestep_limit(p.base)
Interior.device_support_descriptor(p::ConstructiveBatchCells)=Interior.device_support_descriptor(p.state,p.geometry)
@inline batch_parameters(q,w,seeds::NamedTuple{(:cells,:habitat)})=batch_parameters(q,w,seeds.cells)
@inline batch_habitat_parameters(q,w,seeds)=batch_parameters(q,w,seeds.habitat)
@inline batch_anabolic_argument(x::Interior.DeviceCatalyticArgument,w,seeds)=map(a->batch_world(a,w),x)
function Interior.device_cellular_launch!(p::ConstructiveBatchCells,kernel,args...;ndrange)
    dims=ndrange isa Tuple ? ndrange : (ndrange,)
    seeds=(cells=p.seeds,habitat=p.habitat_seeds)
    CELLULAR_BATCH_KERNELS[kernel](p.backend,128)(args...,seeds;ndrange=(dims...,p.worlds))
end
# Public composition validation, stability and hooks are reused without adapters.
function Interior.device_cellular_step!(w::World,p::ConstructiveBatchCells,dt::Real;
        check_status=true,execution=:projected,anabolism=nothing,lifecycle=:prefix)
    anabolism isa Interior.DeviceCatalyticHabitat && anabolism.cells===p || throw(ArgumentError("constructive batch requires its owned habitat"))
    Interior.device_cellular_validate(w,p)
    execution in (:staged,:projected)||throw(ArgumentError("unknown device execution recipe"))
    lifecycle in (:prefix,:dense)||throw(ArgumentError("unknown batch lifecycle recipe"))
    delta=Float32(dt)
    isfinite(delta)&&0<delta<=Interior.stability_limit(w,Interior.DeviceConstructiveEcology(p,anabolism)) || throw(ArgumentError("invalid constructive batch timestep"))
    p.attempted_steps+=1
    Interior.device_cellular_enqueue!(w,p,delta,execution,anabolism,UInt32(p.attempted_steps),lifecycle)
    check_status&&Interior.device_cellular_status!(p)
    p
end
function install_constructive_batch_kernels!()
    path=joinpath(@__DIR__,"..","src","device_habitat_kernels.jl")
    for ex in Meta.parseall(read(path,String)).args
        ex isa Expr&&ex.head==:macrocall&&ex.args[1]==Symbol("@kernel")||continue
        old,new,generated=batch_kernel_ast(deepcopy(ex);habitat=true)
        Core.eval(Interior,generated)
        CELLULAR_BATCH_KERNELS[getproperty(Interior,old)]=Base.invokelatest(getproperty,Interior,new)
    end
    (base=bytes2hex(sha256(read(joinpath(@__DIR__,"cellular_world_batches.jl")))),
     composition=bytes2hex(sha256(read(@__FILE__))),habitat=Interior.DEVICE_HABITAT_SOURCE_PIN)
end
const CONSTRUCTIVE_BATCH_SOURCE=install_constructive_batch_kernels!()
function batch_upload(backend,a::WorldBatchArray)
    storage=Interior.KernelAbstractions.zeros(backend,eltype(a.storage),size(a.storage)...)
    copyto!(storage,a.storage);WorldBatchArray(storage,a.shape)
end
function batch_rehome(p::CellularWorldBatch,backend)
    backend==p.backend&&return p
    s=map(a->batch_upload(backend,a),p.state);e=map(a->batch_upload(backend,a),p.environment)
    seeds=Interior.KernelAbstractions.zeros(backend,UInt32,p.worlds);copyto!(seeds,Array(p.seeds))
    CellularWorldBatch(s,e,backend,p.parameters,p.geometry,Interior.device_cellular_contract(s,e,p.parameters,p.geometry;resource_workspace=batch_resource_workspace(backend,size(e.fields,1),p.worlds;inventory_type=eltype(e.fields)),lifecycle_workspace=batch_lifecycle_workspace(backend,Int(p.geometry.capacity),vec(Array(Interior.device_lifecycle_workspace(p).highwater.storage)))),p.attempted_steps,p.worlds,seeds)
end
function constructive_batch_carrier(base,arrays,feedarrays,q,fq,habitat_seeds)
    c=ConstructiveBatchCells(base,habitat_seeds)
    feed=(arrays=feedarrays,parameters=fq,argument=merge(feedarrays,(status=c.environment.status,)))
    argument=Interior.device_habitat_argument(c,arrays,q)
    ownership=(cells=c,arrays=arrays,feedstock=feed,parameters=q,argument=argument,
        array_shapes=map(size,arrays),feed_shapes=map(size,feedarrays))
    h=Interior.DeviceCatalyticHabitat(c,arrays,feed,q,argument,ownership)
    Interior.DeviceConstructiveEcology(c,h)
end
function constructive_world_batch(ps::AbstractVector;backend=nothing)
    isempty(ps)&&throw(ArgumentError("empty constructive ensemble"))
    backend===nothing&&(backend=first(ps).cells.backend)
    (backend isa CPU||nameof(typeof(backend))===:ROCBackend)||throw(ArgumentError("constructive batch requires CPU or AMD"))
    reference=first(ps).habitat;q=reference.parameters;fq=reference.feedstock.parameters
    for p in ps
        h=p.habitat;o=h.ownership
        h.cells===p.cells&&o.cells===p.cells&&h.arrays===o.arrays&&h.feedstock===o.feedstock&&h.argument===o.argument || throw(ArgumentError("constructive input ownership differs"))
        merge(h.parameters,(seed=q.seed,))==q&&h.feedstock.parameters==fq || throw(ArgumentError("constructive physical parameters differ"))
    end
    base=batch_rehome(CellularWorldBatch(getproperty.(ps,:cells)),backend)
    arrays=pack_world_arrays([p.habitat.arrays for p in ps]);feed=pack_world_arrays([p.habitat.feedstock.arrays for p in ps])
    if backend!=first(ps).cells.backend
        arrays=map(a->batch_upload(backend,a),arrays);feed=map(a->batch_upload(backend,a),feed)
    end
    seeds=Interior.KernelAbstractions.zeros(backend,UInt32,length(ps));copyto!(seeds,UInt32[p.habitat.parameters.seed for p in ps])
    constructive_batch_carrier(base,arrays,feed,q,fq,seeds)
end
function constructive_world_batch_snapshot(p::Interior.DeviceConstructiveEcology{<:ConstructiveBatchCells},i;check_status=true)
    check_status&&Interior.device_cellular_status!(p.cells)
    c=p.cells;h=p.habitat
    1<=i<=c.worlds||throw(BoundsError())
    owned(group)=map(a->copy(selectdim(Array(a.storage),ndims(a.storage),i)),group)
    (cells=cellular_world_batch_snapshot(c.base,i),arrays=owned(h.arrays),feedstock=owned(h.feedstock.arrays),
        parameters=merge(h.parameters,(seed=Array(c.habitat_seeds)[i],)),feed_parameters=h.feedstock.parameters)
end
# Rebuild only ordinary typed owners around downloaded arrays. The strict
# production checkpoint writer/reader remains the schema and validation oracle.
function constructive_batch_host_carrier(x)
    s=x.cells
    c=Interior.DeviceCellular(s.state,s.environment,CPU(),s.parameters,s.geometry,s.attempted_steps,
        Interior.device_cellular_contract(s.state,s.environment,s.parameters,s.geometry))
    a=x.arrays;q=x.parameters
    feed=(arrays=x.feedstock,parameters=x.feed_parameters,argument=merge(x.feedstock,(status=c.environment.status,)))
    arg=Interior.device_habitat_argument(c,a,q)
    own=(cells=c,arrays=a,feedstock=feed,parameters=q,argument=arg,array_shapes=map(size,a),feed_shapes=map(size,x.feedstock))
    Interior.DeviceConstructiveEcology(c,Interior.DeviceCatalyticHabitat(c,a,feed,q,arg,own))
end
function constructive_batch_resource_nonalias(c,arrays,feedarrays)
    work=Interior.device_resource_workspace(c)
    for scratch in (work.grants,(a for level in work.levels for a in (level.high,level.low))...,Interior.device_lifecycle_workspace(c).highwater)
        for physical in (values(arrays)...,values(feedarrays)...)
            Base.mightalias(scratch.storage,physical.storage)&&throw(ArgumentError("batch resource scratch aliases habitat or raw food"))
        end
    end
    true
end
function constructive_batch_validate_ownership(p)
    c=p.cells;h=p.habitat;o=h.ownership
    h.cells===c&&o.cells===c&&h.arrays===o.arrays&&h.feedstock===o.feedstock&&h.parameters===o.parameters&&h.argument===o.argument || throw(ArgumentError("constructive batch ownership changed"))
    map(size,h.arrays)==o.array_shapes&&map(size,h.feedstock.arrays)==o.feed_shapes || throw(ArgumentError("constructive batch shapes changed"))
    h.argument.status===c.environment.status&&h.feedstock.argument.status===c.environment.status || throw(ArgumentError("constructive batch failure latch alias changed"))
    constructive_batch_resource_nonalias(c,h.arrays,h.feedstock.arrays)
    nothing
end
function constructive_batch_validate_worlds(worlds,p)
    length(worlds)==p.cells.worlds||throw(DimensionMismatch("one World per batch member required"))
    length(unique(objectid.(worlds)))==length(worlds)||throw(ArgumentError("batch members require distinct World clocks"))
    for w in worlds
        Interior.device_cellular_validate(w,p.cells)
        !hasproperty(w.fields,:nutrient)||all(iszero,Array(w.fields.nutrient))||throw(ArgumentError("zero unrelated world nutrient before batch execution"))
    end
    nothing
end
function constructive_world_batch_checkpoint(worlds,p::Interior.DeviceConstructiveEcology{<:ConstructiveBatchCells};members=collect(eachindex(worlds)))
    constructive_batch_validate_worlds(worlds,p)
    indices=Int.(collect(members))
    !isempty(indices)&&length(unique(indices))==length(indices)&&all(i->1<=i<=length(worlds),indices) || throw(ArgumentError("invalid selected batch members"))
    constructive_batch_validate_ownership(p)
    snapshots=map(indices) do i
        Interior.device_cellular_validate(worlds[i],p.cells)
        host=constructive_batch_host_carrier(constructive_world_batch_snapshot(p,i;check_status=false))
        w=World(worlds[i].grid,CPU(),NamedTuple(),NamedTuple(),worlds[i].time,worlds[i].seed)
        # The ordinary writer rejects this member if latched or incomplete;
        # failed neighbours do not prevent exporting explicitly selected peers.
        Interior.device_habitat_checkpoint(w,host)
    end
    (format="Interior.ExperimentalConstructiveWorldBatch",version=1,julia_version=VERSION,
        adapter_source=CONSTRUCTIVE_BATCH_SOURCE,source_member_count=p.cells.worlds,
        member_indices=indices,members=snapshots)
end
function constructive_world_batch_restore(worlds,x)
    x.format=="Interior.ExperimentalConstructiveWorldBatch"&&x.version==1&&x.julia_version==VERSION || throw(ArgumentError("unsupported constructive batch checkpoint"))
    x.adapter_source==CONSTRUCTIVE_BATCH_SOURCE||throw(ArgumentError("constructive batch adapter source differs"))
    length(worlds)==length(x.members)>0||throw(DimensionMismatch("batch member/world count differs"))
    length(x.member_indices)==length(x.members)&&length(unique(x.member_indices))==length(x.members)&&all(i->1<=i<=x.source_member_count,x.member_indices) || throw(ArgumentError("invalid checkpoint member indices"))
    length(unique(objectid.(worlds)))==length(worlds)||throw(ArgumentError("batch restore requires distinct World clocks"))
    backend=first(worlds).backend
    (backend isa CPU||nameof(typeof(backend))===:ROCBackend)||throw(ArgumentError("constructive batch restore requires CPU or AMD"))
    all(w->w.backend==backend,worlds)||throw(ArgumentError("mixed restore backends"))
    all(w->!hasproperty(w.fields,:nutrient)||all(iszero,Array(w.fields.nutrient)),worlds)||throw(ArgumentError("zero unrelated world nutrient before restore"))
    # Validate and restore ordinary payloads on CPU first. No standalone GPU
    # populations are allocated in addition to the packed GPU carrier.
    members=map(eachindex(worlds)) do i
        w=World(worlds[i].grid,CPU(),NamedTuple(),NamedTuple(),0.,worlds[i].seed)
        Interior.device_habitat_restore(w,x.members[i])
    end
    result=constructive_world_batch(members;backend)
    Interior.validate(first(worlds),result)
    for i in eachindex(worlds)
        worlds[i].time=x.members[i].cells.world_time;worlds[i].seed=x.members[i].cells.world_seed
    end
    result
end

"""Advance the shared composition once and each successfully completed World clock.

All members share a submitted step index. Failed members remain latched; use
check_status=false for diagnostic continuation of healthy peers, then explicitly
select only complete peers when checkpointing. No failed state is resumed.
"""
function constructive_world_batch_step!(worlds,p::Interior.DeviceConstructiveEcology{<:ConstructiveBatchCells},dt;
        check_status=true,execution=:projected,lifecycle=:prefix)
    constructive_batch_validate_worlds(worlds,p)
    constructive_batch_validate_ownership(p)
    before=copy(Array(p.cells.environment.metadata.storage)[6,:])
    # Do not globally reject an already-latched member before its healthy peers
    # can run. Every original numerical kernel honors its own world's latch.
    Interior.device_cellular_step!(first(worlds),p.cells,dt;check_status=false,execution,lifecycle,anabolism=p.habitat)
    Interior.KernelAbstractions.synchronize(p.cells.backend)
    after=Array(p.cells.environment.metadata.storage)[6,:]
    for i in eachindex(worlds)
        after[i]==before[i]+1&&(worlds[i].time+=Float64(Float32(dt)))
    end
    check_status&&Interior.device_cellular_status!(p.cells)
    p
end
