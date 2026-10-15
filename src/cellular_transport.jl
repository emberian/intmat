"""
    CellStencil(grid, positions, radius)

Sparse, periodic grid–cell exchange with a fixed physical support radius.
`positions` is 3×N in physical coordinates (wrapped into `[0, extent)`).
Grid sites are cell centers; nz=1 uses only x and y. The radial kernel
`max(0, 1 - (distance/radius)^2)^2` is normalized separately for each cell.
If no center lies inside a sub-grid support, the nearest center receives weight
one: unresolved supports therefore retain grid-scale spatial error.

The arrays and sparse site lookup are reused by `rebuild_cell_stencil!`.
Storage and exchange work are O(number of support edges), with no cell×grid
matrix. Exchange uses mutable scratch, so do not share a stencil concurrently.
"""
struct CellStencil
    grid::Grid
    radius::Float64
    cell_volume::Float64
    offsets::Vector{Int}
    sites::Vector{Int}
    weights::Vector{Float64}
    owners::Vector{Int}
    next_edge::Vector{Int}
    site_order::Vector{Int}
    site_heads::Dict{Int,Int}
    proposals::Vector{Float64}
end

function CellStencil(grid::Grid, positions::AbstractMatrix, radius::Real)
    r=Float64(radius)
    isfinite(r) && r>0 || throw(ArgumentError("exchange radius must be positive and finite"))
    nd=grid.size[3]==1 ? 2 : 3
    volume=prod(Float64(grid.extent[d])/grid.size[d] for d in 1:nd)
    stencil=CellStencil(grid,r,volume,Int[],Int[],Float64[],Int[],Int[],
        Int[],Dict{Int,Int}(),Float64[])
    rebuild_cell_stencil!(stencil,positions)
end

# At most n candidates, even when the periodic support spans the whole domain.
function _cell_stencil_axis(x,r,L,n)
    2r>=L && return 1:n
    h=L/n
    lo=ceil(Int,(x-r)/h+0.5); hi=floor(Int,(x+r)/h+0.5)
    hi-lo+1>=n && return 1:n
    (mod1(i,n) for i in lo:hi)
end

"""Rebuild physical supports after cell movement, birth, or removal, reusing storage."""
function rebuild_cell_stencil!(st::CellStencil, positions::AbstractMatrix)
    size(positions,1)==3 || throw(DimensionMismatch("cell positions must be 3×N"))
    all(isfinite,positions) || throw(ArgumentError("cell positions must be finite"))
    for a in (st.offsets,st.sites,st.weights,st.owners,st.next_edge,st.site_order)
        empty!(a)
    end
    empty!(st.site_heads)
    dims=st.grid.size; extent=Float64.(st.grid.extent)
    nd=dims[3]==1 ? 2 : 3
    linear=LinearIndices(dims)
    function edge!(site,weight,cell)
        push!(st.sites,site); push!(st.weights,weight); push!(st.owners,cell)
        previous=get(st.site_heads,site,0)
        previous==0 && push!(st.site_order,site)
        push!(st.next_edge,previous)
        st.site_heads[site]=length(st.sites)
    end
    push!(st.offsets,1)
    for cell in axes(positions,2)
        xyz=ntuple(d->mod(Float64(positions[d,cell]),extent[d]),3)
        candidates=ntuple(d->d>nd ? (1:1) :
            _cell_stencil_axis(xyz[d],st.radius,extent[d],dims[d]),3)
        start=length(st.sites)+1
        for k in candidates[3], j in candidates[2], i in candidates[1]
            indices=(i,j,k)
            distance=sqrt(sum(1:nd) do d
                delta=abs((indices[d]-0.5)*extent[d]/dims[d]-xyz[d])
                min(delta,extent[d]-delta)^2
            end)
            if distance<st.radius
                weight=(1-(distance/st.radius)^2)^2
                weight>0 && edge!(linear[i,j,k],weight,cell)
            end
        end
        if length(st.sites)<start
            nearest=ntuple(d->d>nd ? 1 : min(dims[d],floor(Int,xyz[d]*dims[d]/extent[d])+1),3)
            edge!(linear[nearest...],1.0,cell)
        end
        total=sum(@view st.weights[start:end])
        for e in start:length(st.sites); st.weights[e]/=total; end
        push!(st.offsets,length(st.sites)+1)
    end
    resize!(st.proposals,length(st.sites))
    # Allocation overwrites this scratch before reading it. Initialize it here
    # so a checkpoint never serializes allocator-dependent newly grown slots.
    fill!(st.proposals,0.)
    st
end

function _check_cell_stencil(pool::ResourcePool,st::CellStencil)
    all(a->size(a)==st.grid.size,pool.fields) || throw(DimensionMismatch("resource/stencil grid"))
    nothing
end

"""Sample concentration (grid-site mass / active-dimensional site volume)."""
function sample_concentrations(pool::ResourcePool,st::CellStencil)
    _check_cell_stencil(pool,st)
    result=zeros(Float64,length(pool.fields),length(st.offsets)-1)
    for cell in axes(result,2), s in axes(result,1)
        for e in st.offsets[cell]:st.offsets[cell+1]-1
            result[s,cell]+=st.weights[e]*Float64(pool.fields[s][st.sites[e]])
        end
        result[s,cell]/=st.cell_volume
    end
    result
end

"""Sample species×3×cell physical concentration gradients.
Periodic central differences of the field concentration are averaged with the
same support weights. This is a local finite-volume approximation, not the
exact derivative of the radial interpolant. Singleton axes have zero gradient.
"""
function sample_concentration_gradients(pool::ResourcePool,st::CellStencil)
    _check_cell_stencil(pool,st)
    result=zeros(Float64,length(pool.fields),3,length(st.offsets)-1)
    dims=st.grid.size; cartesian=CartesianIndices(dims)
    for cell in axes(result,3), e in st.offsets[cell]:st.offsets[cell+1]-1
        index=Tuple(cartesian[st.sites[e]])
        for d in 1:3
            dims[d]==1 && continue
            plus=CartesianIndex(ntuple(k->k==d ? mod1(index[k]+1,dims[k]) : index[k],3))
            minus=CartesianIndex(ntuple(k->k==d ? mod1(index[k]-1,dims[k]) : index[k],3))
            scale=st.weights[e]/(2*(Float64(st.grid.extent[d])/dims[d])*st.cell_volume)
            for s in axes(result,1)
                field=pool.fields[s]
                result[s,d,cell]+=scale*(Float64(field[plus])-Float64(field[minus]))
            end
        end
    end
    result
end

"""
    allocate_cell_resources!(pool, stencil, requests)

Allocate species×cell requested MASS synchronously. Each request is distributed
over its normalized kernel weights; all demands on the same species/site are
scaled proportionally if they exceed supply. Unmet demand is not rerouted to
other sites. This preserves fair competition independently of cell order, up to
floating point roundoff. Float64 grants and `pool.withdrawn` record the actual
Float32 field deltas; withdrawals round toward the request cap.
"""
function allocate_cell_resources!(pool::ResourcePool,st::CellStencil,requests::AbstractMatrix)
    _check_cell_stencil(pool,st)
    size(requests)==(length(pool.fields),length(st.offsets)-1) || throw(DimensionMismatch("cell resource requests"))
    all(x->isfinite(x)&&x>=0&&isfinite(Float64(x)),requests) || throw(ArgumentError("requests must be finite nonnegative Float64 masses"))
    grants=zeros(Float64,size(requests))
    for s in eachindex(pool.fields)
        for e in eachindex(st.sites)
            st.proposals[e]=Float64(requests[s,st.owners[e]])*st.weights[e]
        end
        field=pool.fields[s]
        for site in st.site_order
            # Scaling by the largest proposal avoids overflow for finite requests.
            largest=0.0; e=st.site_heads[site]
            while e!=0
                largest=max(largest,st.proposals[e]); e=st.next_edge[e]
            end
            largest==0 && continue
            scaled=0.0; e=st.site_heads[site]
            while e!=0
                scaled+=st.proposals[e]/largest; e=st.next_edge[e]
            end
            before=Float64(field[site]); take=min(before,largest*scaled)
            after=Float32(max(0.0,before-take))
            before-Float64(after)>take && (after=nextfloat(after))
            actual=before-Float64(after)
            field[site]=after; pool.withdrawn[s]+=actual
            e=st.site_heads[site]
            while e!=0
                grants[s,st.owners[e]]+=actual*(st.proposals[e]/largest)/scaled
                e=st.next_edge[e]
            end
        end
    end
    grants
end

"""Deposit mass over one cell's normalized support; return its actual field delta.
Unrepresentable fractions remain with the caller. `pool.returned` is updated.
"""
function deposit_cell_resource!(pool::ResourcePool,species,st::CellStencil,cell::Integer,amount::Real)
    _check_cell_stencil(pool,st)
    1<=cell<length(st.offsets) || throw(BoundsError(st.offsets,cell))
    mass=Float64(amount)
    isfinite(mass)&&mass>=0 || throw(ArgumentError("deposit must be finite nonnegative"))
    s=resource_index(pool,species); field=pool.fields[s]
    edges=st.offsets[cell]:st.offsets[cell+1]-1
    # Check overflow before changing any site.
    for e in edges
        isfinite(Float32(Float64(field[st.sites[e]])+mass*st.weights[e])) || throw(ArgumentError("resource overflow"))
    end
    actual=0.0
    for e in edges
        remaining=mass-actual
        actual+remaining>mass && (remaining=prevfloat(remaining))
        wanted=min(remaining,mass*st.weights[e])
        actual+=resource_add!(pool,s,st.sites[e],wanted)
    end
    actual
end
