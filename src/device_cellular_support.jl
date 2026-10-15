# Borrowed support descriptors are immutable NamedTuples so KA adapts only arrays.
# Iterators are constructed inside kernels; no host callback or temporary edge list.
@inline device_support_descriptor(s,g) = hasproperty(s,:sites) ?
    (sites=s.sites,weights=s.support,count=s.support_count) :
    (indices=s.axis_indices,distance2=s.axis_distance2,normalization=s.support_normalization,
     fallback=s.support_fallback,count=s.support_count,candidates=g.candidates,dims=g.dims,
     radius2=g.support_radius2)
device_support_descriptor(p::DeviceCellular)=device_support_descriptor(p.state,p.geometry)
struct DeviceSupportEdges{D}
    descriptor::D
    cell::Int32
end
@inline device_support_edges(descriptor,i)=DeviceSupportEdges(descriptor,Int32(i))
Base.IteratorSize(::Type{<:DeviceSupportEdges})=Base.SizeUnknown()
Base.eltype(::Type{<:DeviceSupportEdges})=Tuple{Int32,Float32}
# Dense edges advance linearly. Axis edges carry indices rather than dividing a
# flattened cursor for every candidate; the x-fastest candidate order is identical.
@inline function Base.iterate(it::DeviceSupportEdges)
    hasproperty(it.descriptor,:sites) ? Base.iterate(it,Int32(1)) :
        Base.iterate(it,(Int32(1),Int32(1),Int32(1)))
end
@inline function Base.iterate(it::DeviceSupportEdges,cursor::Int32)
    d=it.descriptor;i=it.cell
    cursor>d.count[i] && return nothing
    ((d.sites[i,cursor],d.weights[i,cursor]),cursor+Int32(1))
end
@inline function Base.iterate(it::DeviceSupportEdges,state::NTuple{3,Int32})
    d=it.descriptor;i=it.cell
    d.count[i]==0 && return nothing
    fallback=d.fallback[i]
    if fallback>0
        return state==(Int32(1),Int32(1),Int32(1)) ?
            ((fallback,1f0),(Int32(1),Int32(1),Int32(2))) : nothing
    end
    nx=d.candidates[1];ny=d.candidates[2];nz=d.candidates[3]
    x,y,z=state
    while z<=nz
        r2=(d.distance2[i,x,1]+d.distance2[i,y,2])+d.distance2[i,z,3]
        inside=r2<d.radius2;site=Int32(0);weight=0f0
        if inside
            ix=d.indices[i,x,1];iy=d.indices[i,y,2];iz=d.indices[i,z,3]
            site=ix+d.dims[1]*((iy-Int32(1))+d.dims[2]*(iz-Int32(1)))
            weight=(1f0-r2/d.radius2)^2/d.normalization[i]
        end
        x+=Int32(1)
        if x>nx
            x=Int32(1);y+=Int32(1)
            if y>ny
                y=Int32(1);z+=Int32(1)
            end
        end
        inside && return ((site,weight),(x,y,z))
    end
    nothing
end
