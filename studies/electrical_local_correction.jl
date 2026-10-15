# Isolated temporal-reuse experiment. No default electrical methods overwritten.
include(joinpath(@__DIR__, "branch_solver_prototype.jl"))

"""Red-black smoothing of a bounded patch; external values remain the warm field.
The global residual and ordinary MG fallback, not patch boundaries, certify it.
"""
function electrical_patch!(u,mask,I;radius=8,sweeps=16)
 n=size(u,1);lo=ntuple(d->max(2,I[d]-radius),3);hi=ntuple(d->min(n-1,I[d]+radius),3)
 updates=0
 for sweep in 1:sweeps,color in 0:1,k in lo[3]:hi[3],j in lo[2]:hi[2]
  first=lo[1]+mod(color-lo[1]-j-k,2)
  @inbounds for i in first:2:hi[1]
   if !mask[i,j,k]
    u[i,j,k]=(u[i-1,j,k]+u[i+1,j,k]+u[i,j-1,k]+u[i,j+1,k]+u[i,j,k-1]+u[i,j,k+1])*(1f0/6f0)
    updates+=1
   end
  end
 end
 updates
end
"""Check exactly the stencils touched by the marked site and patch writes."""
function electrical_patch_residual(u,mask,I,radius)
 n=size(u,1);lo=ntuple(d->max(2,I[d]-radius-1),3);hi=ntuple(d->min(n-1,I[d]+radius+1),3)
 r=0f0
 @inbounds for k in lo[3]:hi[3],j in lo[2]:hi[2],i in lo[1]:hi[1]
  mask[i,j,k] && continue
  value=u[i-1,j,k]+u[i+1,j,k]+u[i,j-1,k]+u[i,j+1,k]+u[i,j,k-1]+u[i,j,k+1]-6u[i,j,k]
  r=max(r,abs(value)/6)
 end
 r
end
function electrical_patch_solve!(levels,I;radius=8,sweeps=16,tol=2e-5,maxcycles=200,
 previous_bound=nothing)
 f=first(levels);updates=electrical_patch!(f.u,f.occupied,I;radius,sweeps)
 # A bound is reusable ONLY when this site and patch were the only writes to
 # fine u/mask since the previous certified constant-coefficient harmonic solve.
 # Every untouched stencil is bitwise unchanged. Unknown edits require a full scan.
 r=if previous_bound!==nothing && isfinite(previous_bound) && previous_bound<=tol
  max(previous_bound,Float64(electrical_patch_residual(f.u,f.occupied,I,radius)))
 else
  Float64(residual!(f.residual,f.u,f.rhs,f.occupied))
 end
 r<=tol && return (cycles=0,residual=r,patch_updates=updates)
 s=mgsolve!(levels;tol,maxcycles)
 (cycles=s.cycles,residual=s.residual,patch_updates=updates)
end

function local_correction_study(;n=33,events=32,tol=2e-5,radius=8,sweeps=16)
 solve_tol=.8tol
 mask=specimen(n);base=process(n,mask);mgsolve!(hierarchy(base.potential,base.occupied);tol=solve_tol)
 control=deepcopy(base);candidate=deepcopy(base)
 a=hierarchy(control.potential,control.occupied);b=hierarchy(candidate.potential,candidate.occupied)
 rows=[];matched=true;bound=Interior.potential_residual(candidate)
 for event in 1:events
  I=select_site(control,(1.,1.,1.));J=select_site(candidate,(1.,1.,1.));matched &= I==J
  # Apply the same event in both arms; report independent-choice agreement.
  Interior.electrical_mark!(control,I);Interior.electrical_mark!(candidate,I)
  control.events+=1;candidate.events+=1;mark_hierarchy!(a,I);mark_hierarchy!(b,I)
  reference=mgsolve!(a;tol=solve_tol);trial=electrical_patch_solve!(b,I;radius,sweeps,tol=solve_tol,previous_bound=bound);bound=trial.residual
  independent=Interior.potential_residual(candidate)
  independent<=tol || error("global residual gate failed")
  push!(rows,(event=event,site=Tuple(I),choice_match=I==J,base_cycles=reference.cycles,
   patch_cycles=trial.cycles,patch_updates=trial.patch_updates,residual=independent,
   max_difference=maximum(abs.(control.potential.-candidate.potential))))
  control.reached_ground && break
 end
 (n=n,tolerance=tol,radius=radius,sweeps=sweeps,events=length(rows),all_choices_match=matched,
  base_cycles=sum(r.base_cycles for r in rows),patch_cycles=sum(r.patch_cycles for r in rows),
  patch_updates=sum(r.patch_updates for r in rows),rows=rows)
end
if abspath(PROGRAM_FILE)==@__FILE__
 n=length(ARGS)>0 ? parse(Int,ARGS[1]) : 33
 for (radius,sweeps) in ((4,8),(8,16),(12,24))
  row=local_correction_study(;n,events=32,radius,sweeps)
  println(JSON3.write(row));flush(stdout)
 end
end
