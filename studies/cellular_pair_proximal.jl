# Opt-in operator-splitting study, not part of Interior's default include tree.
mutable struct _CellularRadialCache
    pairs::Vector{Tuple{Int,Int}}
    origin::Matrix{Float64}
    radii::Vector{Float64}
    limit::Float64
    builds::Int
    extent::NTuple{3,Float64}
    dimensions::Int
end
_CellularRadialCache()=_CellularRadialCache(Tuple{Int,Int}[],zeros(3,0),Float64[],0.,0,(0.,0.,0.),0)
function _study_potential_pairs(cache::_CellularRadialCache,grid,x,radii)
    extent=Tuple(Float64.(grid.extent));dimensions=_cellular_ndim(grid)
    if cache.extent!=extent||cache.dimensions!=dimensions||size(cache.origin)!=size(x)||!isequal(cache.radii,radii)||
        !_cellular_inside_pair_skin(grid,x,cache.origin,cache.limit)
        cache.limit=isempty(radii) ? 0. : .2minimum(radii)
        cache.pairs=cellular_pairs(grid,x,radii;range=4cache.limit)
        cache.origin=copy(x);cache.radii=copy(radii);cache.builds+=1
        cache.extent=extent;cache.dimensions=dimensions
    end
    cache.pairs
end
function (cache::_CellularRadialCache)(args...;kwargs...)
    cellular_pair_relax(args...;kwargs...,cache)
end
struct _CellularRadialFlow
    i::Int
    j::Int
    rest::Float64
    stiffness::Float64
    tensile::Bool
end
function _cellular_radial_flow!(grid,ids,x,edge,h)
    i=edge.i;j=edge.j;delta=_cellular_delta(grid,x,i,j);distance=_cellular_length(delta)
    gap=distance-edge.rest
    (edge.tensile ? gap>0 : gap<0) || return
    normal=distance>0 ? delta./distance : (ids[i]<ids[j] ? (1.,0.,0.) : (-1.,0.,0.))
    # Unit mobility for both cells: relative separation relaxes at rate 2k.
    change=-expm1(-2edge.stiffness*h)*(-gap)/2
    nd=_cellular_ndim(grid)
    for d in 1:nd
        L=Float64(grid.extent[d]);move=change*normal[d]
        x[d,i]=mod(x[d,i]-move,L);x[d,j]=mod(x[d,j]+move,L)
    end
end

"""Study-only symmetric exact radial-flow splitting of passive soft contacts.

Each edge receives two half-time exact two-body overdamped flows. The forward
and reverse composition retains unit per-body mobility; it is not a projection
or degree-scaled surrogate. Contact discovery occurs before each sweep, with
fresh discovery for the final true-potential gate. Cross-edge coupling can make
a large split step increase total energy: that step is rejected unchanged by
the caller. This is not a backward-Euler nonlinear solve, so `residual` is NaN.
"""
function cellular_pair_relax(grid,ids,y,radii,bonds::CellBonds,h;
    repulsion=10.,stiffness=2.,max_evaluations=10000,cache=nothing,kwargs...)
    _cellular_validate_geometry(grid,y,radii);_cellular_validate_bonds(bonds)
    isfinite(h)&&h>0&&max_evaluations>=2 || throw(ArgumentError("positive h and two objective evaluations required"))
    length(ids)==length(radii)&&length(unique(ids))==length(ids) || throw(ArgumentError("unique cell IDs required"))
    all(k->isfinite(k)&&k>=0,(repulsion,stiffness)) || throw(ArgumentError("invalid stiffness"))
    index=Dict(id=>i for (i,id) in enumerate(ids));x=Float64.(y)
    nd=_cellular_ndim(grid)
    for i in axes(x,2),d in 1:3;x[d,i]=d>nd ? 0. : mod(x[d,i],Float64(grid.extent[d]));end
    pairs=cache===nothing ? cellular_pairs(grid,x,radii) : _study_potential_pairs(cache,grid,x,radii)
    before=_cellular_energy_forces(grid,ids,x,radii,bonds,index,repulsion,stiffness,Val(false),pairs).energy
    # The complete initial neighbor graph is retained through both half sweeps;
    # an edge may activate as earlier neighboring flows move its endpoints.
    order=issorted(ids) ? eachindex(pairs) : sortperm(pairs;by=ij->minmax(ids[ij[1]],ids[ij[2]]))
    edges=[_CellularRadialFlow(pairs[e][1],pairs[e][2],Float64(radii[pairs[e][1]]+radii[pairs[e][2]]),repulsion,false) for e in order]
    for e in sortperm(bonds.endpoints)
        a,b=bonds.endpoints[e];i=get(index,a,0);j=get(index,b,0)
        (i==0||j==0)&&continue
        push!(edges,_CellularRadialFlow(i,j,bonds.rest_lengths[e],stiffness,true))
    end
    for edge in edges;_cellular_radial_flow!(grid,ids,x,edge,h/2);end
    for edge in Iterators.reverse(edges);_cellular_radial_flow!(grid,ids,x,edge,h/2);end
    final_pairs=cache===nothing ? nothing : _study_potential_pairs(cache,grid,x,radii)
    after=_cellular_energy_forces(grid,ids,x,radii,bonds,index,repulsion,stiffness,Val(false),final_pairs).energy
    (;positions=x,converged=isfinite(after)&&after<=before,residual=NaN,
      iterations=1,cg_iterations=0,evaluations=2,potential_before=before,
      potential_after=after,relaxation_energy_loss=before-after,
      proximal_objective=NaN,displacement=x-y)
end
