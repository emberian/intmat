# Dense GPU-capable continuation of the frozen CPU storage tranche.
module DeviceMemoryDense
using Interior
isdefined(Main,:DeviceMemorySpecialization) || Base.include(Main,joinpath(@__DIR__,"device_memory_specialization.jl"))
const M=Main.DeviceMemorySpecialization
const KA=Interior.KernelAbstractions
const Adapt=KA.Adapt
const SOURCE_PIN=(adapter=bytes2hex(Interior.SHA.sha256(read(@__FILE__))),storage=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"device_memory_specialization.jl")))),cellular=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat=Interior.DEVICE_HABITAT_SOURCE_PIN)

"""Logical zero-tail columns with a device-safe latched invariant failure."""
struct DeviceColumns{T,A,S} <: AbstractMatrix{T}
    data::A
    columns::Int
    status::S
end
DeviceColumns(data::AbstractMatrix{T},columns,status) where T=DeviceColumns{T,typeof(data),typeof(status)}(data,columns,status)
Base.size(a::DeviceColumns)=(size(a.data,1),a.columns)
Base.IndexStyle(::Type{<:DeviceColumns})=IndexCartesian()
@inline Base.getindex(a::DeviceColumns{T},i::Integer,k::Integer) where T=k<=size(a.data,2) ? a.data[i,k] : zero(T)
@inline function Base.setindex!(a::DeviceColumns{T},value,i::Integer,k::Integer) where T
    if k<=size(a.data,2)
        a.data[i,k]=value
    elseif !isequal(value,zero(T))
        KA.@atomic a.status[1]|=Int32(4096)
    end
    value
end
Adapt.adapt_structure(to,a::DeviceColumns)=DeviceColumns(Adapt.adapt(to,a.data),a.columns,Adapt.adapt(to,a.status))
function Adapt.adapt_structure(to,a::M.MaskColumns{T}) where T
    columns=Adapt.adapt(to,a.columns)
    M.MaskColumns{T,typeof(columns)}(a.rows,columns)
end
KA.get_backend(a::DeviceColumns)=KA.get_backend(a.data)
Base.dataids(a::DeviceColumns)=(Base.dataids(a.data)...,Base.dataids(a.status)...)
Base.dataids(a::M.MaskColumns)=Tuple(Iterators.flatten(c isa M.UniformColumn ? () : Base.dataids(c) for c in a.columns))
Base.dataids(::M.UniformVector)=()
# Explicit bulk readout; never perform forbidden host scalar indexing of GPU storage.
function Base.Array(a::DeviceColumns{T}) where T
    data=Array(a.data);out=zeros(T,size(a));out[:,1:size(data,2)].=data;out
end
function Base.Array(a::M.MaskColumns{T}) where T
    out=Matrix{T}(undef,size(a))
    for k in eachindex(a.columns)
        c=a.columns[k];out[:,k].=c isa M.UniformColumn ? c.value : Array(c)
    end
    out
end
Base.Array(a::M.UniformVector)=fill(a.value,a.rows)
Base.all(f,a::DeviceColumns)=all(f,a.data)&&(size(a.data,2)==a.columns||f(zero(eltype(a))))
Base.all(f,a::M.MaskColumns)=all(c->c isa M.UniformColumn ? f(c.value) : all(f,c),a.columns)
Base.all(f,a::M.UniformVector)=a.rows==0||f(a.value)

mutable struct Transfer{B,S}
    backend::B
    status::S
    seen::IdDict{Any,Any}
end
function transfer_array!(t,a::Array)
    get!(t.seen,a) do
        out=KA.zeros(t.backend,eltype(a),size(a)...);copyto!(out,a);out
    end
end
function (t::Transfer)(a)
    if a isa Array
        transfer_array!(t,a)
    elseif a isa SubArray
        view(t(parent(a)),Base.parentindices(a)...)
    elseif a isa M.ActiveColumns
        DeviceColumns(t(a.data),a.columns,t.status)
    elseif a isa M.MaskColumns
        columns=map(c->c isa M.UniformColumn ? c : t(c),a.columns)
        M.MaskColumns{eltype(a),typeof(columns)}(a.rows,columns)
    elseif a isa M.UniformVector || a isa M.UniformColumn
        a
    elseif a isa NamedTuple || a isa Tuple
        map(t,a)
    else
        throw(ArgumentError("unsupported dense transfer role $(typeof(a))"))
    end
end

"""Direct target allocation from canonical host checkpoint. No dense target carrier.
Host validation currently uses the canonical CPU restore and may have a dense
host peak. Large fresh initialization needs a later allocator hook; no hidden
claim that checkpoint construction or CPU validation has become sparse.
"""
function restore(w::World,checkpoint;shared_resource=true)
    cpu=World(w.grid,CPU(),NamedTuple(),NamedTuple(),0.,w.seed)
    canonical=Interior.device_habitat_restore(cpu,checkpoint)
    compact=M.specialize(cpu,canonical;pages=false)
    transfer(w,cpu,compact;shared_resource)
end
function transfer(w::World,cpu,p;shared_resource=true)
    c=p.cells;h=p.habitat
    status=KA.zeros(w.backend,Int32,size(c.environment.status)...);copyto!(status,c.environment.status)
    t=Transfer(w.backend,status,IdDict{Any,Any}(c.environment.status=>status))
    arena=shared_resource ? KA.zeros(w.backend,Float32,size(c.environment.fields,1),6) : nothing
    if shared_resource
        scratch=view(arena,:,1:3);copyto!(scratch,c.environment.scratch)
        t.seen[c.environment.scratch]=scratch
    end
    state=t(c.state);env=t(c.environment)
    # The earned core workspace explicitly leases the arena only through finish.
    contract=if shared_resource
        work=Interior.DeviceResourceWorkspace(w.backend,size(env.fields,1);grants=arena)
        Interior.device_cellular_contract(state,env,c.parameters,c.geometry;resource_workspace=work,
            resource_transient_aliases=(:scratch,:demand,:scale,:returns,:occupancy))
    else
        Interior.device_cellular_contract(state,env,c.parameters,c.geometry)
    end
    cells=Interior.DeviceCellular(state,env,w.backend,c.parameters,c.geometry,c.attempted_steps,contract)
    arrays=t(h.arrays);feed_arrays=t(h.feedstock.arrays)
    feed=(arrays=feed_arrays,parameters=h.feedstock.parameters,argument=merge(feed_arrays,(status=status,)))
    arg=Interior.device_habitat_argument(cells,arrays,h.parameters)
    ownership=(cells=cells,arrays=arrays,feedstock=feed,parameters=h.parameters,argument=arg,
        array_shapes=map(size,arrays),feed_shapes=map(size,feed_arrays))
    habitat=Interior.DeviceCatalyticHabitat(cells,arrays,feed,h.parameters,arg,ownership)
    result=Interior.DeviceConstructiveEcology(cells,habitat)
    w.time=cpu.time;w.seed=cpu.seed
    Interior.validate(w,result)
    result
end
"""Reachable backing allocation bytes, parameterized for CPU or native GPU arrays."""
function retained_bytes(x,isstorage;storage_key=identity,storage_bytes=a->length(a)*sizeof(eltype(a)))
    seen=IdDict{Any,Bool}();bytes=Ref(0)
    function visit(a)
        if isstorage(a)
            key=storage_key(a)
            haskey(seen,key)&&return;seen[key]=true;bytes[]+=storage_bytes(a)
        elseif a isa SubArray || a isa Base.ReshapedArray
            visit(parent(a))
        elseif a isa Tuple || a isa NamedTuple || a isa AbstractDict
            foreach(visit,values(a))
        elseif isstructtype(typeof(a)) && !(a isa Number||a isa Symbol||a isa String||a isa Function)
            haskey(seen,a)&&return;seen[a]=true
            for k in 1:fieldcount(typeof(a));visit(getfield(a,k));end
        end
    end
    visit(x);bytes[]
end
end
