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

"""Reconstructible scratch for ordered parallel materialized-support construction.

The physical support layout and arithmetic order are unchanged. Serial execution
owns no scratch. Candidate scratch is never scientific checkpoint inventory.
"""
struct DeviceSupportWorkspace{A,B,C}
    arrays::A
    backend::B
    contract::C
end
function DeviceSupportWorkspace(backend::KernelAbstractions.Backend,state,g)
    hasproperty(state,:sites) || throw(ArgumentError("parallel support requires materialized support"))
    dims=(Int(g.capacity),prod(Int.(g.candidates)))
    size(state.sites)==dims || throw(ArgumentError("materialized support capacity differs"))
    z(T,d...)=KernelAbstractions.zeros(backend,T,d...)
    arrays=(sites=z(Int32,dims...),weights=z(Float32,dims...),totals=z(Float32,Int(g.capacity)))
    DeviceSupportWorkspace(arrays,backend,(arrays=arrays,shapes=map(size,arrays),dims=dims))
end
function device_support_workspace_validate(work::Nothing,backend,g)
    true
end
function device_support_workspace_validate(work::DeviceSupportWorkspace,backend,g)
    a=work.arrays;c=work.contract
    work.backend==backend && a===c.arrays && map(size,a)==c.shapes &&
        c.dims==(Int(g.capacity),prod(Int.(g.candidates))) && g.support_strategy==0 &&
        size(a.sites)==size(a.weights)==c.dims && size(a.totals)==(Int(g.capacity),) &&
        eltype(a.sites)===Int32 && eltype(a.weights)===Float32 && eltype(a.totals)===Float32 &&
        all(x->KernelAbstractions.get_backend(x)==backend,values(a)) ||
        throw(ArgumentError("parallel support workspace role/backend/shape changed"))
    true
end
function device_support_workspace_nonalias(work::DeviceSupportWorkspace,state,env,resource,lifecycle)
    a=values(work.arrays)
    borrowed=(values(state)...,values(env)...,resource.grants,
        (x for level in resource.levels for x in (level.high,level.low))...,lifecycle.highwater)
    for i in eachindex(a)
        any(j->Base.mightalias(a[i],a[j]),1:i-1) && throw(ArgumentError("support scratch roles alias"))
        any(x->Base.mightalias(a[i],x),borrowed) && throw(ArgumentError("support scratch aliases another process role"))
    end
    true
end
device_support_workspace_nonalias(::Nothing,args...)=true
device_support_workspace(p)=p.ownership.stages.support
device_support_execution(p)=device_support_workspace(p)===nothing ? :serial : :parallel
function device_support_execution_validate(execution,strategy)
    execution in (:serial,:parallel) || throw(ArgumentError("support_execution must be :serial or :parallel"))
    execution===:parallel && strategy!==:materialized && throw(ArgumentError("parallel support requires materialized strategy"))
    execution
end

@kernel function dc_support_candidates!(s,e,q,g,sites,weights)
    i,c=@index(Global,NTuple)
    if e.status[1]==0 && s.active[i]==1
        t=c-1;dx=rem(t,g.candidates[1]);t=div(t,g.candidates[1]);dy=rem(t,g.candidates[2]);dz=div(t,g.candidates[2])
        lo=ntuple(d->d>g.nd ? Int32(1) : Int32(floor((s.position[i,d]-q.exchange_radius)/g.spacing[d]+.5f0)),3)
        pos=(dc_period(lo[1]+dx,g.dims[1]),dc_period(lo[2]+dy,g.dims[2]),dc_period(lo[3]+dz,g.dims[3]))
        r2=0f0
        for d in 1:Int(g.nd)
            delta=abs((pos[d]-.5f0)*g.spacing[d]-s.position[i,d]);delta=min(delta,g.extent[d]-delta);r2+=delta^2
        end
        inside=r2<q.exchange_radius^2
        sites[i,c]=inside ? dc_index(pos...,g.dims) : Int32(0)
        weights[i,c]=inside ? (1f0-r2/q.exchange_radius^2)^2 : 0f0
    end
end
@kernel function dc_support_compress!(s,e,g,sites,weights,totals)
    i=@index(Global,Linear)
    if e.status[1]==0
        s.support_count[i]=0
        if s.active[i]==1
            count=Int32(0);total=0f0
            for c in 1:size(sites,2)
                site=sites[i,c]
                if site!=0
                    count+=Int32(1);weight=weights[i,c]
                    s.sites[i,count]=site;s.support[i,count]=weight;total+=weight
                end
            end
            if count==0
                pos=ntuple(d->d>g.nd ? Int32(1) : min(g.dims[d],Int32(floor(s.position[i,d]/g.spacing[d]))+1),3)
                count=Int32(1);s.sites[i,1]=dc_index(pos...,g.dims);s.support[i,1]=1f0;total=1f0
            end
            totals[i]=total;s.support_count[i]=count
        end
    end
end
@kernel function dc_support_normalize!(s,e,totals)
    i,c=@index(Global,NTuple)
    if e.status[1]==0 && s.active[i]==1 && c<=s.support_count[i]
        s.support[i,c]/=totals[i]
    end
end

function device_support_launch!(p,::Nothing,s,e,q,g)
    dc_support!(p.backend,128)(s,e,q,g;ndrange=Int(g.capacity))
end
function device_support_launch!(p,work::DeviceSupportWorkspace,s,e,q,g)
    a=work.arrays
    dc_support_candidates!(p.backend,128)(s,e,q,g,a.sites,a.weights;ndrange=size(a.sites))
    dc_support_compress!(p.backend,128)(s,e,g,a.sites,a.weights,a.totals;ndrange=Int(g.capacity))
    dc_support_normalize!(p.backend,128)(s,e,a.totals;ndrange=size(a.sites))
end
