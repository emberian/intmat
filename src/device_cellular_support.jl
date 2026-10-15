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
@inline function Base.iterate(it::DeviceSupportEdges,cursor::Int32=Int32(1))
    d=it.descriptor;i=it.cell
    if hasproperty(d,:sites)
        cursor>d.count[i] && return nothing
        return ((d.sites[i,cursor],d.weights[i,cursor]),cursor+Int32(1))
    end
    d.count[i]==0 && return nothing
    fallback=d.fallback[i]
    if fallback>0
        return cursor==1 ? ((fallback,1f0),Int32(2)) : nothing
    end
    nx=d.candidates[1];ny=d.candidates[2];limit=nx*ny*d.candidates[3]
    while cursor<=limit
        offset=cursor-Int32(1);ix=rem(offset,nx)+Int32(1)
        iy=rem(div(offset,nx),ny)+Int32(1);iz=div(offset,nx*ny)+Int32(1)
        r2=(d.distance2[i,ix,1]+d.distance2[i,iy,2])+d.distance2[i,iz,3]
        cursor+=Int32(1)
        if r2<d.radius2
            x=d.indices[i,ix,1];y=d.indices[i,iy,2];z=d.indices[i,iz,3]
            site=x+d.dims[1]*((y-Int32(1))+d.dims[2]*(z-Int32(1)))
            weight=(1f0-r2/d.radius2)^2/d.normalization[i]
            return ((site,weight),cursor)
        end
    end
    nothing
end
