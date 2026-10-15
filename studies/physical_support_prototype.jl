module PhysicalSupportPrototype
using Interior
export AxisSupport, foreach_support_edge, sample_support!, support_storage_bytes, SupportAllocationScratch, allocate_support!

"""Bounded-memory geometry for the existing normalized radial exchange law.

Caches one-dimensional periodic coordinates and squared distances, plus each
cell's normalization. Three-dimensional support edges/weights are regenerated.
Storage is O(N*sum(axis candidates)); exchange work still scales with physical
support volume. Float64 mode follows CellStencil precision; Float32 mode follows
the device's geometry arithmetic. No support radius is changed with resolution.
This is an isolated prototype, not a replacement process or checkpoint format.
"""
struct AxisSupport{T,G}
    grid::G
    radius::T
    volume::T
    counts::NTuple{3,Int}
    indices::Array{Int32,3}  # candidate, axis, cell; CPU-local cell slabs
    distance2::Array{T,3}
    normalization::Vector{T}
    edge_count::Vector{Int32}
    fallback::Vector{Int32}
end
function AxisSupport(grid,positions,radius;precision::Type{T}=Float64) where T
    T in (Float32,Float64) || throw(ArgumentError("Float32/Float64 precision required"))
    size(positions,1)==3 && all(isfinite,positions) || throw(ArgumentError("finite physical positions required"))
    isfinite(radius)&&radius>0 || throw(ArgumentError("positive finite radius required"))
    dims=grid.size;nd=dims[3]==1 ? 2 : 3;n=size(positions,2)
    extent=T.(grid.extent);h=ntuple(d->extent[d]/T(dims[d]),3);r=T(radius)
    isfinite(r)&&r>0 && all(x->isfinite(T(x)),positions) || throw(ArgumentError("support inputs must remain finite and radius positive in selected precision"))
    counts=ntuple(d->d>nd ? 1 : min(dims[d],2ceil(Int,r/h[d])+2),3)
    a=maximum(counts);indices=zeros(Int32,a,3,n);distance2=zeros(T,a,3,n)
    normalization=zeros(T,n);edge_count=zeros(Int32,n);fallback=zeros(Int32,n)
    st=AxisSupport(grid,r,prod(h[d] for d in 1:nd),counts,indices,distance2,normalization,edge_count,fallback)
    r2=r*r
    @inbounds for cell in 1:n
        xyz=ntuple(d->mod(T(positions[d,cell]),extent[d]),3)
        lo=ntuple(d->d>nd ? 1 : floor(Int,(xyz[d]-r)/h[d]+T(.5)),3)
        for d in 1:3,j in 1:counts[d]
            pos=mod1(lo[d]+j-1,dims[d]);indices[j,d,cell]=pos
            delta=d>nd ? zero(T) : abs((T(pos)-T(.5))*h[d]-xyz[d])
            delta=min(delta,extent[d]-delta)
            distance2[j,d,cell]=delta*delta
        end
        total=zero(T);count=Int32(0)
        for iz in 1:counts[3],iy in 1:counts[2],ix in 1:counts[1]
            squared=(distance2[ix,1,cell]+distance2[iy,2,cell])+distance2[iz,3,cell]
            if squared<r2
                weight=(one(T)-squared/r2)^2
                total+=weight;count+=1
            end
        end
        if count==0
            nearest=ntuple(d->d>nd ? 1 : min(dims[d],floor(Int,xyz[d]/h[d])+1),3)
            fallback[cell]=nearest[1]+dims[1]*((nearest[2]-1)+dims[2]*(nearest[3]-1))
            total=one(T);count=1
        end
        normalization[cell]=total;edge_count[cell]=count
    end
    st
end

"""Call `f(site, normalized_weight, i,j,k)` in the support's spatial order."""
@inline function foreach_support_edge(f,st::AxisSupport{T},cell) where T
    dims=st.grid.size
    @inbounds if st.fallback[cell]!=0
        site=Int(st.fallback[cell]);i=mod(site-1,dims[1])+1
        j=mod((site-1)÷dims[1],dims[2])+1;k=(site-1)÷(dims[1]*dims[2])+1
        f(site,one(T),i,j,k);return nothing
    end
    r2=st.radius*st.radius
    @inbounds for iz in 1:st.counts[3],iy in 1:st.counts[2],ix in 1:st.counts[1]
        squared=(st.distance2[ix,1,cell]+st.distance2[iy,2,cell])+st.distance2[iz,3,cell]
        if squared<r2
            weight=(one(T)-squared/r2)^2/st.normalization[cell]
            i=Int(st.indices[ix,1,cell]);j=Int(st.indices[iy,2,cell]);k=Int(st.indices[iz,3,cell])
            site=i+dims[1]*((j-1)+dims[2]*(k-1))
            f(site,weight,i,j,k)
        end
    end
    nothing
end

"""Fuse concentration and central-difference gradient sampling over identical supports.
Outputs are species×cell and species×3×cell, and may be reused across calls.
"""
function sample_support!(concentrations,gradients,fields,st::AxisSupport)
    n=length(st.normalization);species=length(fields);dims=st.grid.size
    size(concentrations)==(species,n) && size(gradients)==(species,3,n) || throw(DimensionMismatch("sample output shape"))
    all(f->size(f)==dims,fields) || throw(DimensionMismatch("field/support grid"))
    fill!(concentrations,0);fill!(gradients,0)
    h=ntuple(d->Float64(st.grid.extent[d])/dims[d],3);volume=Float64(st.volume)
    @inbounds for cell in 1:n
        foreach_support_edge(st,cell) do site,weight,i,j,k
            xyz=(i,j,k)
            for species in eachindex(fields)
                field=fields[species]
                concentrations[species,cell]+=Float64(weight)*Float64(field[site])
            end
            for d in 1:3
                dims[d]==1 && continue
                stride=d==1 ? 1 : d==2 ? dims[1] : dims[1]*dims[2]
                plus=xyz[d]==dims[d] ? site-(dims[d]-1)*stride : site+stride
                minus=xyz[d]==1 ? site+(dims[d]-1)*stride : site-stride
                scale=Float64(weight)/(2h[d]*volume)
                for species in eachindex(fields)
                    field=fields[species]
                    gradients[species,d,cell]+=scale*(Float64(field[plus])-Float64(field[minus]))
                end
            end
        end
        for s in eachindex(fields);concentrations[s,cell]/=volume;end
    end
    concentrations,gradients
end
support_storage_bytes(st::AxisSupport)=sizeof(st.indices)+sizeof(st.distance2)+sizeof(st.normalization)+sizeof(st.edge_count)+sizeof(st.fallback)
"""Grid-sized reusable scratch; no per-cell×support proposals or site dictionary."""
struct SupportAllocationScratch
    largest::Vector{Float64}
    scaled::Vector{Float64}
    actual::Vector{Float64}
    active_sites::Vector{Int32}
end
SupportAllocationScratch(st::AxisSupport)=SupportAllocationScratch(zeros(prod(st.grid.size)),zeros(prod(st.grid.size)),zeros(prod(st.grid.size)),Int32[])

"""Same synchronous proportional capped allocation as CellStencil.
Requests and grants are species×cell masses. No unmet request is rerouted.
Float32 withdrawals round toward the request cap; actual deltas fund grants.
The three-pass largest/scaled sum retains overflow resistance for finite requests.
"""
function allocate_support!(pool,st::AxisSupport,requests,workspace::SupportAllocationScratch)
    n=length(st.normalization);species=length(pool.fields)
    size(requests)==(species,n) || throw(DimensionMismatch("support requests"))
    all(x->isfinite(x)&&x>=0&&isfinite(Float64(x)),requests) || throw(ArgumentError("finite nonnegative Float64 requests required"))
    all(f->size(f)==st.grid.size,pool.fields) || throw(DimensionMismatch("resource/support grid"))
    length(workspace.largest)==prod(st.grid.size) || throw(DimensionMismatch("support scratch"))
    grants=zeros(Float64,species,n)
    for s in 1:species
        empty!(workspace.active_sites)
        @inbounds for cell in 1:n
            requested=Float64(requests[s,cell]);requested==0 && continue
            foreach_support_edge(st,cell) do site,weight,i,j,k
                proposal=requested*Float64(weight)
                if proposal>0
                    if workspace.largest[site]==0
                        push!(workspace.active_sites,Int32(site));workspace.scaled[site]=0
                    end
                    workspace.largest[site]=max(workspace.largest[site],proposal)
                end
            end
        end
        @inbounds for cell in 1:n
            requested=Float64(requests[s,cell]);requested==0 && continue
            foreach_support_edge(st,cell) do site,weight,i,j,k
                proposal=requested*Float64(weight)
                proposal>0 && (workspace.scaled[site]+=proposal/workspace.largest[site])
            end
        end
        field=pool.fields[s]
        @inbounds for site in workspace.active_sites
            before=Float64(field[site]);take=min(before,workspace.largest[site]*workspace.scaled[site])
            after=Float32(max(0.,before-take));before-Float64(after)>take && (after=nextfloat(after))
            workspace.actual[site]=before-Float64(after);field[site]=after
            pool.withdrawn[s]+=workspace.actual[site]
        end
        @inbounds for cell in 1:n
            requested=Float64(requests[s,cell]);requested==0 && continue
            foreach_support_edge(st,cell) do site,weight,i,j,k
                proposal=requested*Float64(weight)
                if proposal>0
                    grants[s,cell]+=workspace.actual[site]*(proposal/workspace.largest[site])/workspace.scaled[site]
                end
            end
        end
        @inbounds for site in workspace.active_sites;workspace.largest[site]=0.;end
    end
    grants
end

end
