# Opt-in study: independent worlds share launches, never state or material.
# Math is sourced from the loaded implementation; only indexing is transformed.
using Interior, SHA
isdefined(Interior,:DeviceCellular) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_cellular.jl"))

struct WorldBatchArray{T,N,A<:AbstractArray{T}} <: AbstractArray{T,N}
    storage::A
    shape::NTuple{N,Int}
end
Base.size(a::WorldBatchArray)=a.shape
Interior.KernelAbstractions.get_backend(a::WorldBatchArray)=Interior.KernelAbstractions.get_backend(a.storage)
const WorldBatchAdapt=Interior.KernelAbstractions.Adapt
function WorldBatchAdapt.adapt_structure(to,a::WorldBatchArray)
    WorldBatchArray(WorldBatchAdapt.adapt(to,a.storage),a.shape)
end
# A world is a contiguous slice. Encode its offset once instead of routing
# every scalar access through SubArray's general multidimensional indexing.
struct WorldBatchSlice{T,N,A} <: DenseArray{T,N}
    storage::A
    shape::NTuple{N,Int}
    offset::Int
end
Base.size(a::WorldBatchSlice)=a.shape
Base.IndexStyle(::Type{<:WorldBatchSlice})=IndexLinear()
@inline function Base.getindex(a::WorldBatchSlice,i::Int)
    @boundscheck checkbounds(a,i)
    @inbounds a.storage[a.offset+i]
end
@inline function Base.setindex!(a::WorldBatchSlice,v,i::Int)
    @boundscheck checkbounds(a,i)
    @inbounds a.storage[a.offset+i]=v
end
Base.pointer(a::WorldBatchSlice)=pointer(a.storage,a.offset+1)
Base.pointer(a::WorldBatchSlice,i::Integer)=pointer(a.storage,a.offset+i)
Base.unsafe_convert(::Type{Ptr{T}},a::WorldBatchSlice{T}) where T=pointer(a)
@inline batch_world(a::WorldBatchArray{T,N},w) where {T,N}=WorldBatchSlice{T,N,typeof(a.storage)}(a.storage,a.shape,(Int(w)-1)*length(a))
@inline batch_world(x,w)=x
# Project only fields a kernel actually touches; eagerly mapping every field
# produces dozens of unused SubArray objects at every logical thread.
struct WorldBatchTuple{T}
    packed::T
    world::Int
end
@inline batch_world(x::NamedTuple,w)=WorldBatchTuple(x,Int(w))
@inline Base.getproperty(x::WorldBatchTuple,k::Symbol)=batch_world(getproperty(getfield(x,:packed),k),getfield(x,:world))
Base.propertynames(x::WorldBatchTuple)=keys(getfield(x,:packed))
mutable struct CellularWorldBatch{S,E,B,Q,G,O,R}
    state::S
    environment::E
    backend::B
    parameters::Q
    geometry::G
    ownership::O
    attempted_steps::Int
    worlds::Int
    seeds::R
end
function pack_world_arrays(xs)
    map(keys(first(xs))) do k
        arrays=map(x->getproperty(x,k),xs);a=first(arrays)
        all(x->size(x)==size(a)&&eltype(x)==eltype(a),arrays) || error("incompatible world arrays")
        storage=cat(arrays...;dims=ndims(a)+1)
        WorldBatchArray(storage,size(a))
    end |> values -> NamedTuple{keys(first(xs))}(values)
end
# Each world's resource reduction has disjoint transient storage. The shared
# production phase sees ordinary per-world shapes through the batch projection.
function batch_resource_workspace(backend,n,worlds;inventory_type=Float32)
    array(rows)=WorldBatchArray(Interior.KernelAbstractions.zeros(backend,inventory_type,rows,6,worlds),(rows,6))
    grants=array(n);levels=NamedTuple[]
    while true
        n=cld(n,128);push!(levels,(high=array(n),low=array(n)))
        n==1&&break
    end
    typed=Tuple(levels)
    contract=(grants=grants,levels=typed,grant_shape=size(grants),level_shapes=map(x->(size(x.high),size(x.low)),typed))
    Interior.DeviceResourceWorkspace(grants,typed,contract)
end
function batch_resource_validate(p)
    work=Interior.device_resource_workspace(p)
    Interior.device_resource_workspace_validate(work)
    arrays=(work.grants,(a for l in work.levels for a in (l.high,l.low))...)
    for (i,a) in enumerate(arrays)
        a isa WorldBatchArray&&size(a.storage)==(size(a)...,p.worlds)&&eltype(a)===eltype(p.environment.fields) ||
            throw(ArgumentError("invalid world-isolated resource storage"))
        Interior.KernelAbstractions.get_backend(a.storage)==p.backend || throw(ArgumentError("resource batch backend differs"))
        for other in arrays[1:i-1]
            Base.mightalias(a.storage,other.storage)&&throw(ArgumentError("resource batch scratch parents alias"))
        end
        for other in (values(p.state)...,values(p.environment)...)
            Base.mightalias(a.storage,other.storage)&&throw(ArgumentError("resource batch scratch aliases process state"))
        end
    end
    true
end
function batch_lifecycle_workspace(backend,capacity,markers)
    values=Int32.(collect(markers))
    all(h->0<=h<=capacity,values)||throw(ArgumentError("batch lifecycle marker exceeds capacity"))
    parent=Interior.KernelAbstractions.zeros(backend,Int32,1,length(values));copyto!(parent,reshape(values,1,:))
    highwater=WorldBatchArray(parent,(1,))
    Interior.DeviceLifecycleWorkspace(highwater,backend,Int(capacity),(highwater=highwater,shape=size(highwater)))
end
function batch_lifecycle_validate(p)
    work=Interior.device_lifecycle_workspace(p);h=work.highwater
    Interior.device_lifecycle_workspace_validate(work,p.backend,Int(p.geometry.capacity))
    h isa WorldBatchArray&&size(h.storage)==(1,p.worlds)||throw(ArgumentError("invalid world-local lifecycle storage"))
    resource=Interior.device_resource_workspace(p)
    others=(values(p.state)...,values(p.environment)...,resource.grants,(a for level in resource.levels for a in (level.high,level.low))...)
    any(a->Base.mightalias(h.storage,a.storage),others)&&throw(ArgumentError("batch lifecycle marker aliases process or resource storage"))
    true
end
function CellularWorldBatch(ps::AbstractVector)
    isempty(ps)&&error("empty ensemble");a=first(ps)
    all(p->p.geometry==a.geometry&&merge(p.parameters,(seed=a.parameters.seed,))==a.parameters&&p.backend==a.backend&&p.attempted_steps==a.attempted_steps,ps)||error("batch requires identical geometry, physical parameters, and step clock")
    s=pack_world_arrays(getproperty.(ps,:state));e=pack_world_arrays(getproperty.(ps,:environment))
    resource=batch_resource_workspace(a.backend,size(e.fields,1),length(ps);inventory_type=eltype(e.fields))
    markers=Int32[only(Array(Interior.device_lifecycle_workspace(p).highwater)) for p in ps]
    lifecycle=batch_lifecycle_workspace(a.backend,Int(a.geometry.capacity),markers)
    o=Interior.device_cellular_contract(s,e,a.parameters,a.geometry;resource_workspace=resource,lifecycle_workspace=lifecycle)
    seeds=Interior.KernelAbstractions.zeros(a.backend,UInt32,length(ps))
    copyto!(seeds,UInt32[p.parameters.seed for p in ps])
    CellularWorldBatch(s,e,a.backend,a.parameters,a.geometry,o,a.attempted_steps,length(ps),seeds)
end
# These study methods do not replace any production method.
function Interior.device_cellular_validate(w::Interior.World,p::CellularWorldBatch)
    w.backend==p.backend && w.grid.size==Tuple(Int.(p.geometry.dims)) && w.grid.extent==p.geometry.extent || error("batch world mismatch")
    p.state===p.ownership.state && p.environment===p.ownership.environment || error("batch arrays replaced")
    batch_resource_validate(p)
    batch_lifecycle_validate(p)
    true
end
Interior.device_cellular_timestep_limit(p::CellularWorldBatch)=begin
    q=p.parameters;g=p.geometry;D=maximum(q.diffusion);S=sum(inv(Float64(g.spacing[d])^2) for d in 1:Int(g.nd))
    min(.25,Float64(q.response_time/2),D>0 ? .45/(D*S) : Inf)
end
function Interior.device_cellular_status!(p::CellularWorldBatch)
    Interior.KernelAbstractions.synchronize(p.backend)
    codes=Array(p.environment.status.storage)
    all(iszero,codes)||throw(ArgumentError("independent-world status bits: $codes"))
    p
end
@inline batch_parameters(q,w,seeds)=hasproperty(q,:seed) ? merge(q,(seed=seeds[w],)) : q
@inline batch_anabolic_argument(x,w,seeds)=x
const CELLULAR_BATCH_KERNELS=IdDict()
function Interior.device_cellular_launch!(p::CellularWorldBatch,kernel,args...;ndrange)
    dims=ndrange isa Tuple ? ndrange : (ndrange,)
    CELLULAR_BATCH_KERNELS[kernel](p.backend,128)(args...,p.seeds;ndrange=(dims...,p.worlds))
end
# Rewrite parsed syntax, not textual bodies. Refuse unfamiliar index forms.
function batch_kernel_ast(ex;habitat=false)
    fn=ex.args[end];sig=fn.args[1];body=fn.args[2]
    name=sig.args[1];sig.args[1]=Symbol("worldbatch_",name)
    indexpos=findfirst(x->x isa Expr&&x.head==:(=)&&x.args[2] isa Expr&&x.args[2].head==:macrocall&&x.args[2].args[1]==Symbol("@index"),body.args)
    indexpos===nothing&&error("no top-level kernel index in $name")
    index=body.args[indexpos];mode=index.args[2].args[end]
    mode in (:Linear,:NTuple)||error("unsupported index mode")
    lhs=index.args[1];indices=lhs isa Symbol ? [lhs] : lhs.args
    index.args[1]=Expr(:tuple,indices...,:batch_world_id)
    index.args[2].args[end]=:NTuple
    # Original kernel argument symbols keep their names after world slicing.
    insertions=Expr[]
    for j in 2:length(sig.args)
        arg=sig.args[j]
        if arg in (:s,:e,:q,:a,:feed,:descriptor,:cq,:anabolism,:grants,:high,:low,:input_high,:input_low,:status,:highwater)
            packed=Symbol("packed_",arg)
            sig.args[j]=packed
            if arg===:anabolism
                push!(insertions,:($arg=$(GlobalRef(@__MODULE__,:batch_anabolic_argument))($packed,batch_world_id,batch_seeds)))
            elseif arg===:q && habitat
                push!(insertions,:($arg=$(GlobalRef(@__MODULE__,:batch_habitat_parameters))($packed,batch_world_id,batch_seeds)))
            elseif arg in (:q,:cq)
                push!(insertions,:($arg=$(GlobalRef(@__MODULE__,:batch_parameters))($packed,batch_world_id,batch_seeds)))
            else
                push!(insertions,:($arg=$(GlobalRef(@__MODULE__,:batch_world))($packed,batch_world_id)))
            end
        end
    end
    push!(sig.args,:batch_seeds)
    splice!(body.args,indexpos+1:indexpos,insertions)
    name,sig.args[1],ex
end
function install_cellular_world_batch!()
    normpath(dirname(pathof(Interior)))==normpath(joinpath(@__DIR__,"..","src")) || error("batch study must use Interior from the same source snapshot")
    source=joinpath(@__DIR__,"..","src","device_cellular_kernels.jl")
    resource_source=joinpath(@__DIR__,"..","src","device_resource_exchange.jl")
    for ex in (Meta.parseall(read(source,String)).args...,Meta.parseall(read(resource_source,String)).args...)
        ex isa Expr&&ex.head==:macrocall&&ex.args[1]==Symbol("@kernel")||continue
        old,new,generated=batch_kernel_ast(deepcopy(ex))
        Core.eval(Interior,generated)
        CELLULAR_BATCH_KERNELS[getproperty(Interior,old)]=Base.invokelatest(getproperty,Interior,new)
    end
    host=joinpath(@__DIR__,"..","src","device_cellular.jl")
    isdefined(Interior,:device_cellular_enqueue!) || error("batch study requires the shared generic enqueue interface")
    (kernels=bytes2hex(sha256(read(source))),resource=bytes2hex(sha256(read(resource_source))),schedule=bytes2hex(sha256(read(host))))
end
const CELLULAR_WORLD_BATCH_SOURCE=install_cellular_world_batch!()
# Only checked entry bookkeeping is local; every physical stage is the common
# enqueue schedule, including its mechanical host loop and extension hooks.
function cellular_world_batch_step!(w::World,p::CellularWorldBatch,dt::Real;
        check_status=true,execution=:projected,anabolism=nothing,lifecycle=:prefix)
    anabolism===nothing || throw(ArgumentError("batched anabolic extension is not implemented"))
    lifecycle in (:prefix,:dense)||throw(ArgumentError("unknown batch lifecycle recipe"))
    Interior.device_cellular_validate(w,p)
    execution in (:staged,:projected) || throw(ArgumentError("unknown device execution recipe"))
    delta=Float32(dt)
    isfinite(delta)&&0<delta<=Interior.device_cellular_timestep_limit(p)||throw(ArgumentError("invalid device cellular timestep"))
    p.attempted_steps+=1
    Interior.device_cellular_enqueue!(w,p,delta,execution,nothing,UInt32(p.attempted_steps),lifecycle)
    check_status&&Interior.device_cellular_status!(p)
    p
end

function cellular_world_batch_snapshot(p::CellularWorldBatch,i::Integer)
    1<=i<=p.worlds || throw(BoundsError())
    # Audit copies are outside timing. Download the packed parent in bulk before
    # host slicing: our custom CPU slice must never scalar-index ROCArray data.
    (state=map(a->copy(selectdim(Array(a.storage),ndims(a.storage),i)),p.state),
        environment=map(a->copy(selectdim(Array(a.storage),ndims(a.storage),i)),p.environment),
        parameters=merge(p.parameters,(seed=Array(p.seeds)[i],)),geometry=p.geometry,attempted_steps=p.attempted_steps)
end
