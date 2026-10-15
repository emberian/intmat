# CPU-only exact storage specialization. Numerical kernels/schedules remain Interior's.
module DeviceMemorySpecialization
using Interior
const KA=Interior.KernelAbstractions

"""Logical columns with a checked immutable zero tail; only active columns stored."""
struct ActiveColumns{T,A<:AbstractMatrix{T}} <: AbstractMatrix{T}
    data::A
    columns::Int
end
Base.size(a::ActiveColumns)=(size(a.data,1),a.columns)
Base.IndexStyle(::Type{<:ActiveColumns})=IndexCartesian()
@inline Base.getindex(a::ActiveColumns{T},i::Int,k::Int) where T=k<=size(a.data,2) ? a.data[i,k] : zero(T)
@inline function Base.setindex!(a::ActiveColumns{T},value,i::Int,k::Int) where T
    if k<=size(a.data,2);a.data[i,k]=value
    else;isequal(value,zero(T))||throw(ArgumentError("nonzero write to specialized inactive column"));end
    value
end

struct UniformColumn{T}
    value::T
end
"""Each immutable mask column is either uniform or retained dense, exactly."""
struct MaskColumns{T,C} <: AbstractMatrix{T}
    rows::Int
    columns::C
end
Base.size(a::MaskColumns)=(a.rows,length(a.columns))
Base.IndexStyle(::Type{<:MaskColumns})=IndexCartesian()
@inline maskvalue(a::UniformColumn,i)=a.value
@inline maskvalue(a::AbstractVector,i)=a[i]
@inline Base.getindex(a::MaskColumns,i::Int,k::Int)=maskvalue(a.columns[k],i)
function mask_columns(a)
    cols=ntuple(size(a,2)) do k
        column=view(a,:,k);value=first(column)
        all(x->isequal(x,value),column) ? UniformColumn(value) : copy(column)
    end
    MaskColumns{eltype(a),typeof(cols)}(size(a,1),cols)
end
struct UniformVector{T} <: AbstractVector{T}
    rows::Int
    value::T
end
Base.size(a::UniformVector)=(a.rows,)
Base.getindex(a::UniformVector,i::Int)=a.value

"""Exact CPU proof-of-representation: implicit zero pages, no deleted histories.
GPU execution requires a bounded atomic page map/pool, not this Dict backend.
"""
struct PagedMatrix{T} <: AbstractMatrix{T}
    rows::Int
    columns::Int
    page_rows::Int
    pages::Dict{Int,Matrix{T}}
end
Base.size(a::PagedMatrix)=(a.rows,a.columns)
Base.IndexStyle(::Type{<:PagedMatrix})=IndexCartesian()
@inline function Base.getindex(a::PagedMatrix{T},i::Int,k::Int) where T
    page=get(a.pages,cld(i,a.page_rows),nothing)
    page===nothing ? zero(T) : page[mod1(i,a.page_rows),k]
end
function Base.setindex!(a::PagedMatrix{T},value,i::Int,k::Int) where T
    key=cld(i,a.page_rows);page=get(a.pages,key,nothing)
    if page===nothing
        isequal(value,zero(T)) && return value
        page=zeros(T,a.page_rows,a.columns);a.pages[key]=page
    end
    page[mod1(i,a.page_rows),k]=value
end
function paged(a::AbstractMatrix;page_rows=16)
    out=PagedMatrix{eltype(a)}(size(a,1),size(a,2),page_rows,Dict{Int,Matrix{eltype(a)}}())
    for k in axes(a,2),i in axes(a,1);out[i,k]=a[i,k];end
    out
end
# Only these small wrappers are CPU-backed in this prototype.
KA.get_backend(::Union{ActiveColumns,MaskColumns,UniformVector,PagedMatrix})=CPU()

function zero_tail(a,active)
    all(x->isequal(x,zero(eltype(a))),view(a,:,active+1:size(a,2))) || throw(ArgumentError("inactive species/ledger column contains data"))
    ActiveColumns(copy(view(a,:,1:active)),size(a,2))
end
function specialize(w,p;pages=true,page_rows=16)
    pages && Threads.nthreads()!=1 && throw(ArgumentError("Dict-page storage requires CPU1; no thread-safe GPU allocator claimed"))
    w.backend isa CPU || throw(ArgumentError("storage prototype currently requires CPU1; no GPU adaptation claimed"))
    c=p.cells;h=p.habitat;e=c.environment;V=prod(Int.(c.geometry.dims))
    # Complete diffusion epochs precede clear_grid; scale overwrites demand;
    # clear_occupancy occurs after waste commit. One 3-column scratch arena.
    scratch=copy(e.scratch)
    env=merge(e,(scratch=scratch,demand=view(scratch,:,1:2),scale=view(scratch,:,1:2),
        returns=view(scratch,:,3),occupancy=view(scratch,:,3),
        source_masks=mask_columns(e.source_masks),sink_masks=mask_columns(e.sink_masks)))
    cells=Interior.DeviceCellular(c.state,env,c.backend,c.parameters,c.geometry,c.attempted_steps,
        Interior.device_cellular_contract(c.state,env,c.parameters,c.geometry))
    f=h.feedstock.arrays
    feed=merge(f,(fields=zero_tail(f.fields,2),scratch=ActiveColumns(view(scratch,:,1:2),3),
        source_masks=mask_columns(f.source_masks),sink_masks=mask_columns(f.sink_masks),
        permeability=all(x->isequal(x,first(f.permeability)),f.permeability) ? UniformVector(V,first(f.permeability)) : f.permeability))
    a=h.arrays
    arrays=merge(a,(grid_ledger=zero_tail(a.grid_ledger,11),grid_ledger_low=zero_tail(a.grid_ledger_low,11)))
    if pages
        # Persistent material/work accounting is retained bit-for-bit, per site.
        arrays=merge(arrays,(grid_ledger=ActiveColumns(paged(view(a.grid_ledger,:,1:11);page_rows),12),
            grid_ledger_low=ActiveColumns(paged(view(a.grid_ledger_low,:,1:11);page_rows),12),
            grid_inert=vec(paged(reshape(a.grid_inert,:,1);page_rows)),
            grid_inert_low=vec(paged(reshape(a.grid_inert_low,:,1);page_rows))))
    end
    feedstock=(arrays=feed,parameters=h.feedstock.parameters,argument=merge(feed,(status=env.status,)))
    arg=Interior.device_habitat_argument(cells,arrays,h.parameters)
    own=(cells=cells,arrays=arrays,feedstock=feedstock,parameters=h.parameters,argument=arg,
        array_shapes=map(size,arrays),feed_shapes=map(size,feed))
    habitat=Interior.DeviceCatalyticHabitat(cells,arrays,feedstock,h.parameters,arg,own)
    result=Interior.DeviceConstructiveEcology(cells,habitat)
    # Source object must no longer be used for stepping: persistent state is borrowed.
    Interior.validate(w,result);result
end

"""Count actual retained Array backing bytes once, following views and aliases."""
function retained_array_bytes(x)
    seen=IdDict{Any,Bool}();bytes=Ref(0)
    function visit(a)
        if a isa Array
            haskey(seen,a)&&return;seen[a]=true;bytes[]+=sizeof(a)
        elseif a isa SubArray || a isa Base.ReshapedArray
            visit(parent(a))
        elseif a isa NamedTuple || a isa Tuple
            foreach(visit,values(a))
        elseif a isa AbstractDict
            foreach(visit,values(a))
        elseif isstructtype(typeof(a)) && !(a isa Number || a isa Symbol || a isa String || a isa Function)
            haskey(seen,a)&&return;seen[a]=true
            for k in 1:fieldcount(typeof(a));visit(getfield(a,k));end
        end
    end
    visit(x);bytes[]
end
const TRANSIENT_ENV=(:scratch,:demand,:scale,:returns,:occupancy)
function persistent_snapshot(p)
    x=Interior.device_habitat_snapshot(p)
    names=Tuple(k for k in keys(x.cells.environment) if !(k in TRANSIENT_ENV))
    cells=merge(x.cells,(environment=NamedTuple{names}(x.cells.environment),))
    feed_names=Tuple(k for k in keys(x.feedstock) if k!==:scratch)
    merge(x,(cells=cells,feedstock=NamedTuple{feed_names}(x.feedstock)))
end
end
