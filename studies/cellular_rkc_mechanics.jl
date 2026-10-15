# Opt-in stabilized explicit passive mechanics. Load after implicit/pair studies.
# RKC2 recurrence: Sommeijer/Shampine/Verwer MAS-R9715, section2.
# https://ir.cwi.nl/pub/4764/04764D.pdf
function cellular_rkc_coefficients(stages::Int;damping=2/13)
 stages>=2&&isfinite(damping)&&damping>0||throw(ArgumentError("at least two stages and positive damping required"))
 w0=1+damping/stages^2
 t=zeros(stages+1);d=similar(t);dd=similar(t);fill!(d,0);fill!(dd,0)
 t[1]=1;t[2]=w0;d[2]=1
 for j in 2:stages
  t[j+1]=2w0*t[j]-t[j-1]
  d[j+1]=2t[j]+2w0*d[j]-d[j-1]
  dd[j+1]=4d[j]+2w0*dd[j]-dd[j-1]
 end
 b=zeros(stages+1)
 for j in 2:stages;b[j+1]=dd[j+1]/d[j+1]^2;end
 b[1]=b[2]=b[3];w1=d[end]/dd[end]
 (;stages,w0,w1,t,b,stability_interval=(w0+1)/w1)
end
function cellular_rkc_step(force,x,h,coeff)
 stages=coeff.stages;w0=coeff.w0;w1=coeff.w1;t=coeff.t;b=coeff.b
 f0=force(x);previous_previous=zeros(size(x));previous=(b[2]*w1*h).*f0
 for j in 2:stages
  current_force=force(x+previous)
  mu=2b[j+1]*w0/b[j];nu=-b[j+1]/b[j-1]
  mu_tilde=2b[j+1]*w1/b[j];gamma_tilde=-(1-b[j]*t[j])*mu_tilde
  current=mu.*previous.+nu.*previous_previous.+(mu_tilde*h).*current_force.+(gamma_tilde*h).*f0
  previous_previous,previous=previous,current
 end
 x+previous
end
mutable struct _CellularRKCCache
 graph::_CellularRadialCache
 force_evaluations::Int
 energy_evaluations::Int
 accepted_stages::Vector{Int}
 stage_limit::Int
end
_CellularRKCCache(;stage_limit=256)=_CellularRKCCache(_CellularRadialCache(),0,0,Int[],stage_limit)
_study_potential_pairs(cache::_CellularRKCCache,args...)=_study_potential_pairs(cache.graph,args...)
(cache::_CellularRKCCache)(args...;kwargs...)=cellular_rkc_relax(args...;cache,kwargs...)

function _rkc_stage_distance2(z)
 largest=0.
 @inbounds for i in axes(z,2)
  distance=0.
  for d in axes(z,1);distance+=z[d,i]^2;end
  largest=max(largest,distance)
 end
 largest
end
_rkc_failed(x,initial,evaluations,stages,reason)=(positions=x,converged=false,residual=NaN,
 potential_before=initial,potential_after=Inf,evaluations=evaluations,iterations=0,
 cg_iterations=0,reason=reason,stages=stages)

"""Second-order stabilized explicit passive flow, with new numerical semantics.
The negative-real-axis stability interval is not a nonlinear contact theorem.
Every stage rediscovers contacts through the conservative skin; excessive stage
travel or final potential increase rejects the whole trial. Bonds refresh only
at the enclosing study wrapper's macro endpoints. No default method is replaced.
"""
function cellular_rkc_relax(grid,ids,y,radii,bonds::CellBonds,h;
 repulsion=10.,stiffness=2.,max_evaluations=10000,cache=_CellularRKCCache(),bond_columns=nothing,kwargs...)
 _cellular_validate_geometry(grid,y,radii);_cellular_validate_bonds(bonds)
 isfinite(h)&&h>0&&max_evaluations>=2||throw(ArgumentError("positive finite interval and two evaluations required"))
 length(ids)==length(radii)&&length(unique(ids))==length(ids)||throw(ArgumentError("unique cell IDs required"))
 all(k->isfinite(k)&&k>=0,(repulsion,stiffness))||throw(ArgumentError("invalid stiffness"))
 x=Float64.(y);n=length(ids);index=Dict(id=>i for(i,id)in enumerate(ids))
 f=similar(x);rows=zeros(n);pairs=_study_potential_pairs(cache,grid,x,radii)
 initial=_cellular_energy_forces!(grid,ids,x,radii,bonds,index,repulsion,stiffness,Val(true),pairs,f,rows,bond_columns)
 cache.force_evaluations+=1;evaluations=1
 bound=isempty(rows) ? 0. : 2maximum(rows)
 stages=max(2,ceil(Int,sqrt((1+1.1h*bound)/.65)))
 stages<=cache.stage_limit&&max_evaluations>=stages+1||return _rkc_failed(x,initial.energy,evaluations,stages,:work_budget)
 coeff=cellular_rkc_coefficients(stages)
 h*bound<=coeff.stability_interval||return _rkc_failed(x,initial.energy,evaluations,stages,:linear_stability_estimate)
 f0=copy(f);previous_previous=zeros(size(x));previous=(coeff.b[2]*coeff.w1*h).*f0
 limit=isempty(radii) ? Inf : .2minimum(radii)
 for j in 2:stages
  _rkc_stage_distance2(previous)<=limit^2||return _rkc_failed(x,initial.energy,evaluations,stages,:stage_travel)
  z=x+previous;pairs=_study_potential_pairs(cache,grid,z,radii)
  _cellular_energy_forces!(grid,ids,z,radii,bonds,index,repulsion,stiffness,Val(true),pairs,f,rows,bond_columns)
  cache.force_evaluations+=1;evaluations+=1
  mu=2coeff.b[j+1]*coeff.w0/coeff.b[j];nu=-coeff.b[j+1]/coeff.b[j-1]
  mu_tilde=2coeff.b[j+1]*coeff.w1/coeff.b[j];gamma_tilde=-(1-coeff.b[j]*coeff.t[j])*mu_tilde
  current=mu.*previous.+nu.*previous_previous.+(mu_tilde*h).*f.+(gamma_tilde*h).*f0
  previous_previous,previous=previous,current
 end
 _rkc_stage_distance2(previous)<=limit^2||return _rkc_failed(x,initial.energy,evaluations,stages,:stage_travel)
 z=x+previous;pairs=_study_potential_pairs(cache,grid,z,radii)
 after=_cellular_energy_forces(grid,ids,z,radii,bonds,index,repulsion,stiffness,Val(false),pairs,bond_columns).energy
 cache.energy_evaluations+=1;evaluations+=1
 accepted=isfinite(after)&&after<=initial.energy
 accepted&&push!(cache.accepted_stages,stages)
 (;positions=z,converged=accepted,residual=NaN,potential_before=initial.energy,potential_after=after,
  evaluations,iterations=0,cg_iterations=0,reason=accepted ? :accepted : :potential_increase,stages)
end
