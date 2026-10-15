# Study-only coordinate proximal splitting; load after both other study files.
struct _CellularLocalEdge
    neighbor::Int
    rest::Float64
    stiffness::Float64
    tensile::Bool
end
mutable struct _CellularCoordinateCache
    graph::_CellularRadialCache
    midpoint::Bool
end
_CellularCoordinateCache(;midpoint=false)=_CellularCoordinateCache(_CellularRadialCache(),midpoint)
_study_potential_pairs(cache::_CellularCoordinateCache,args...)=_study_potential_pairs(cache.graph,args...)
(cache::_CellularCoordinateCache)(args...;kwargs...)=cellular_coordinate_relax(args...;kwargs...,cache)
@inline function _coordinate_delta(grid,z,x,j)
    nd=_cellular_ndim(grid)
    ntuple(d->d>nd ? 0. : _cellular_mod_period(x[d,j]-z[d]+Float64(grid.extent[d])/2,Float64(grid.extent[d]))-Float64(grid.extent[d])/2,3)
end
function _coordinate_state(grid,ids,i,z,x,anchor,edges,tau)
    gradient=(z-anchor)/tau;hessian=Matrix{Float64}(I,3,3)/tau;energy=0.
    for edge in edges
        j=edge.neighbor;delta=_coordinate_delta(grid,z,x,j);distance=_cellular_length(delta);gap=distance-edge.rest
        (edge.tensile ? gap>0 : gap<0)||continue
        k=edge.stiffness;energy+=.5k*gap^2
        normal=distance>0 ? delta./distance : (ids[i]<ids[j] ? (1.,0.,0.) : (-1.,0.,0.))
        tangent=distance>0 ? k*gap/distance : 0.
        for d in 1:3
            gradient[d]-=k*gap*normal[d]
            for e in 1:3
                hessian[d,e]+=(k-tangent)*normal[d]*normal[e]+(d==e ? tangent : 0.)
            end
        end
    end
    objective=energy+sum(abs2,z-anchor)/(2tau)
    (;energy,gradient,hessian,objective)
end
function _coordinate_energy(grid,z,x,edges)
    energy=0.
    for edge in edges
        distance=_cellular_length(_coordinate_delta(grid,z,x,edge.neighbor));gap=distance-edge.rest
        (edge.tensile ? gap>0 : gap<0)||continue
        energy+=.5edge.stiffness*gap^2
    end
    energy
end
function _coordinate_inside(grid,z,origin,i,limit)
    distance2=0.
    for d in 1:_cellular_ndim(grid)
        L=Float64(grid.extent[d]);delta=_cellular_mod_period(z[d]-origin[d,i]+L/2,L)-L/2
        distance2+=delta^2
    end
    distance2<=limit^2
end
function _coordinate_body(grid,ids,i,x,edges,tau,graph,tolerance,max_iterations)
    anchor=x[:,i];z=copy(anchor);state=_coordinate_state(grid,ids,i,z,x,anchor,edges,tau)
    initial=state.energy;scale=1+maximum(abs,state.gradient);evaluations=1;iterations=0;reason=:iteration_limit
    for _ in 1:max_iterations
        residual=maximum(abs,state.gradient)
        residual<=tolerance*scale&&return (;converged=true,z,residual,evaluations,iterations,initial_energy=initial,reason=:accepted,threshold=tolerance*scale,chart_margin=graph.limit-sqrt(sum(abs2,_coordinate_delta(grid,z,graph.origin,i))))
        matrix=copy(state.hessian);direction=-tau*state.gradient
        if all(isfinite,matrix)
            factor=cholesky(Symmetric(matrix);check=false)
            if !isposdef(factor)
                lower=minimum(matrix[d,d]-sum(abs(matrix[d,e]) for e in 1:3 if e!=d) for d in 1:3)
                shift=max(0.,-lower+inv(tau))
                for d in 1:3;matrix[d,d]+=shift;end
                factor=cholesky(Symmetric(matrix);check=false)
            end
            isposdef(factor)&&(direction=-(factor\state.gradient))
        end
        slope=_implicit_dot(state.gradient,direction)
        if !(slope<0&&isfinite(slope));reason=:direction;break;end
        alpha=1.;accepted=false
        # A local energy sums many edges; its forward rounding error scales
        # with that count. Residual must still improve, and the final global
        # potential gate permits no increase.
        roundoff=8max(1,length(edges))*eps(max(1.,abs(state.objective),abs(initial)))
        for _ in 1:30
            trial_z=z+alpha*direction
            if _coordinate_inside(grid,trial_z,graph.origin,i,graph.limit)
                trial=_coordinate_state(grid,ids,i,trial_z,x,anchor,edges,tau);evaluations+=1
                armijo=trial.objective<=state.objective+1e-4alpha*slope
                refined=trial.objective<=state.objective+roundoff&&
                    _implicit_dot(trial.gradient,trial.gradient)<.99*_implicit_dot(state.gradient,state.gradient)
                if isfinite(trial.objective)&&trial.objective<=initial+roundoff&&(armijo||refined)
                    z=trial_z;state=trial;accepted=true;break
                end
            end
            alpha/=2
        end
        if !accepted;reason=:line_search;break;end
        iterations+=1
    end
    (;converged=false,z,residual=maximum(abs,state.gradient),evaluations,iterations,initial_energy=initial,reason,threshold=tolerance*scale,chart_margin=graph.limit-sqrt(sum(abs2,_coordinate_delta(grid,z,graph.origin,i))))
end

"""Coordinate splitting with local 3D backward-Euler proximal substeps.

A forward/reverse body sweep uses half-time unit-mobility proximal updates of
all incident pair potentials. Every local objective is checked, so other-edge
coupling cannot spend unreported energy. Zero net-force bodies are fixed points.
The default is first-order (local flows use backward Euler). With
`cache.midpoint=true`, midpoint extrapolation and symmetric body order give
second order on a smooth fixed graph; endpoint energy is checked separately.
Roundoff allowances in local sums scale with incident-edge count, while the
final global potential comparison permits no increase. A conservative
cached graph bounds each local chart; failure requests macrostep reduction.
Residual reports the largest accepted local proximal residual, not a global
backward-Euler residual. Graph energy is checked again after the whole sweep.
"""
function cellular_coordinate_relax(grid,ids,y,radii,bonds::CellBonds,h;
    repulsion=10.,stiffness=2.,tolerance=1e-7,max_iterations=80,
    max_evaluations=10000,cache=_CellularCoordinateCache(),kwargs...)
    _cellular_validate_geometry(grid,y,radii);_cellular_validate_bonds(bonds)
    isfinite(h)&&h>0&&max_evaluations>=2||throw(ArgumentError("positive h and two objective evaluations required"))
    n=length(ids);length(radii)==n&&length(unique(ids))==n||throw(ArgumentError("unique cell IDs required"))
    all(k->isfinite(k)&&k>=0,(repulsion,stiffness))||throw(ArgumentError("invalid stiffness"))
    x=Float64.(y);nd=_cellular_ndim(grid)
    for i in 1:n,d in 1:3;x[d,i]=d>nd ? 0. : mod(x[d,i],Float64(grid.extent[d]));end
    graph=cache.graph;pairs=_study_potential_pairs(graph,grid,x,radii)
    # Reserve half the skin for this sweep's local solves.
    if !_cellular_inside_pair_skin(grid,x,graph.origin,graph.limit/2)
        graph.origin=zeros(3,0);pairs=_study_potential_pairs(graph,grid,x,radii)
    end
    index=Dict(id=>i for (i,id) in enumerate(ids))
    before=_cellular_energy_forces(grid,ids,x,radii,bonds,index,repulsion,stiffness,Val(false),pairs).energy
    adjacent=[_CellularLocalEdge[] for _ in 1:n]
    for (i,j) in pairs
        rest=Float64(radii[i]+radii[j])
        push!(adjacent[i],_CellularLocalEdge(j,rest,repulsion,false));push!(adjacent[j],_CellularLocalEdge(i,rest,repulsion,false))
    end
    for e in sortperm(bonds.endpoints)
        a,b=bonds.endpoints[e];i=get(index,a,0);j=get(index,b,0)
        (i==0||j==0)&&continue
        push!(adjacent[i],_CellularLocalEdge(j,bonds.rest_lengths[e],stiffness,true))
        push!(adjacent[j],_CellularLocalEdge(i,bonds.rest_lengths[e],stiffness,true))
    end
    order=sortperm(ids);iterations=0;local_evaluations=0;residual=0.;converged=true;failed_body=0;failed_details=nothing
    for sequence in (order,Iterators.reverse(order)),i in sequence
        tau=h/2
        anchor=x[:,i]
        result=_coordinate_body(grid,ids,i,x,adjacent[i],cache.midpoint ? tau/2 : tau,graph,tolerance,max_iterations)
        iterations+=result.iterations;local_evaluations+=result.evaluations
        residual=max(residual,result.residual)
        if !result.converged;converged=false;failed_body=ids[i];failed_details=result;break;end
        endpoint=cache.midpoint ? 2result.z-anchor : result.z
        if cache.midpoint
            # Midpoint solves the same local BE proximal equation at half tau.
            # Extrapolation needs its own physical endpoint energy/chart gate.
            start_local=result.initial_energy
            end_local=_coordinate_energy(grid,endpoint,x,adjacent[i])
            if !_coordinate_inside(grid,endpoint,graph.origin,i,graph.limit)||end_local>start_local
                converged=false;failed_body=ids[i];failed_details=(;reason=:midpoint_endpoint);break
            end
        end
        for d in 1:nd;x[d,i]=mod(endpoint[d],Float64(grid.extent[d]));end
    end
    final_pairs=_study_potential_pairs(graph,grid,x,radii)
    after=_cellular_energy_forces(grid,ids,x,radii,bonds,index,repulsion,stiffness,Val(false),final_pairs).energy
    reason=!converged ? :local_solve : !(isfinite(after)&&after<=before) ? :global_energy : :accepted
    (;positions=x,converged=converged&&isfinite(after)&&after<=before,residual,reason,failed_body,failed_details,
      iterations,cg_iterations=0,evaluations=2,local_evaluations,
      potential_before=before,potential_after=after,relaxation_energy_loss=before-after,
      proximal_objective=NaN,displacement=x-y)
end
