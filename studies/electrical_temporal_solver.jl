# Optional solve coordinator for fresh constant-coefficient harmonic producers.
# Not loaded by Interior or any running archive.
include(joinpath(@__DIR__,"electrical_patch_kernel.jl"))
include(joinpath(@__DIR__,"branch_metal_multigrid.jl"))

"""Owns a residual certificate for explicitly tracked fine-field edits.
All external writes to any level, RHS, mask, boundary or operator data MUST call
`electrical_temporal_invalidate!`. Direct array writes cannot be detected without
an expensive scan. Only one known Dirichlet addition may use a patch; other edits
use ordinary MG. The underlying solver is homogeneous constant-coefficient
Laplace, not a variable-conductivity solver.
"""
mutable struct ElectricalTemporalSolver{L,N,P}
 levels::L
 norm::N
 patch_norm::P
 tolerance::Float64
 workgroup::NTuple{3,Int}
 maxcycles::Int
 prior_bound::Float64
 certificate::Bool
 patch_enabled::Bool
 pending_site::Union{Nothing,CartesianIndex{3}}
 pending_count::Int
 patch_pending::Bool
 generation::UInt64
 certified_generation::UInt64
 invalidation_reason::String
end
function ElectricalTemporalSolver(levels;tol=1.6e-5,workgroup=(32,4,2),maxcycles=200)
 isfinite(tol)&&tol>0&&maxcycles>0 || throw(ArgumentError("positive solver controls required"))
 ElectricalTemporalSolver(levels,DeviceMGNorm(first(levels)),similar(first(levels).u,Float32,1),
  Float64(tol),Tuple(Int.(workgroup)),Int(maxcycles),Inf,false,true,nothing,0,false,
  UInt64(0),UInt64(0),"initial state requires global solve")
end
function electrical_temporal_invalidate!(solver::ElectricalTemporalSolver,reason)
 solver.certificate=false;solver.patch_enabled=false;solver.patch_pending=false
 solver.pending_site=nothing;solver.pending_count=0;solver.generation+=1
 solver.invalidation_reason=String(reason)
 solver
end
function electrical_temporal_mark!(solver::ElectricalTemporalSolver,I::CartesianIndex{3})
 n=size(first(solver.levels).u,1)
 all(d->1<I[d]<n,1:3)||throw(ArgumentError("new conductor site must be interior"))
 eligible=solver.certificate&&solver.patch_enabled&&solver.pending_count==0
 device_mark!(solver.levels,I)
 solver.pending_count+=1;solver.pending_site=I
 solver.patch_pending=eligible&&solver.pending_count==1
 solver.certificate=false;solver.generation+=1
 solver
end
function _electrical_temporal_publish!(solver,stats)
 solver.prior_bound=Float64(stats.residual);solver.certificate=true
 solver.certified_generation=solver.generation
 solver.pending_site=nothing;solver.pending_count=0;solver.patch_pending=false
 solver.invalidation_reason=""
 stats
end
function electrical_temporal_solve!(solver::ElectricalTemporalSolver)
 if solver.certificate&&solver.certified_generation==solver.generation
  return (cycles=0,residual=solver.prior_bound,patch_used=false,cached=true)
 end
 try
  stats=if solver.patch_pending&&solver.patch_enabled&&solver.pending_count==1
   f=first(solver.levels)
   patch=electrical_device_patch8!(f.u,f.occupied,solver.pending_site;norm=solver.patch_norm)
   bound=max(solver.prior_bound,patch.residual)
   if bound<=solver.tolerance
    (cycles=0,residual=bound,patch_used=true,cached=false)
   else
    ordinary=device_solve!(solver.levels;tol=solver.tolerance,maxcycles=solver.maxcycles,
     workgroup=solver.workgroup,norm_workspace=solver.norm)
    merge(ordinary,(patch_used=true,cached=false))
   end
  else
   ordinary=device_solve!(solver.levels;tol=solver.tolerance,maxcycles=solver.maxcycles,
    workgroup=solver.workgroup,norm_workspace=solver.norm)
   merge(ordinary,(patch_used=false,cached=false))
  end
  _electrical_temporal_publish!(solver,stats)
 catch
  electrical_temporal_invalidate!(solver,"failed solve requires global retry")
  rethrow()
 end
end
"""Explicit re-enable after unknown edits. Caller must restore the supported
constant-coefficient grounded harmonic operator and synchronize every coarse mask.
The fine RHS check rejects a known incompatible Poisson problem. Certification
always performs an ordinary global solve before local reuse can resume.
"""
function electrical_temporal_recertify!(solver::ElectricalTemporalSolver;homogeneous_operator=false)
 electrical_temporal_invalidate!(solver,"explicit recertification")
 homogeneous_operator===true || throw(ArgumentError("explicit supported-operator assertion required"))
 all(iszero,first(solver.levels).rhs)||throw(ArgumentError("patch reuse requires homogeneous fine RHS"))
 stats=electrical_temporal_solve!(solver)
 solver.patch_enabled=true
 stats
end
