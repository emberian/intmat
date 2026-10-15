#!/usr/bin/env julia
# Replay exact frozen numerical sources; instrument only the pre-event85 capture.
using Interior, SHA, Serialization
using Interior: JSON3
include(joinpath(@__DIR__,"..","studies","electrical_cdf_diagnostic.jl"))
const ELECTRICAL_DIVERGENCE_CONTEXT=Ref{Any}()
function capture_electrical_divergence!(a,b,c,output,rows)
 context=ELECTRICAL_DIVERGENCE_CONTEXT[]
 expected=JSON3.read(read(context.original_results,String)).rows
 length(rows)==84||error("wrong common prefix length")
 all(i->rows[i].baseline_site==Tuple(expected[i].baseline_site)&&rows[i].candidate_site==Tuple(expected[i].candidate_site),1:84)||error("replay differs from frozen original prefix")
 a.p.events==b.p.events==84||error("wrong event counters")
 a.p.occupied==b.p.occupied&&a.p.frontier==b.p.frontier||error("common geometry differs")
 sa,da=ArchiveLaw.frontier_data(a.p,first(a.levels),a.temp,a.dose,a.food,AMDGPU.ROCArray)
 sb,db=ArchiveLaw.frontier_data(b.p,first(b.levels),b.temp,b.dose,b.food,AMDGPU.ROCArray)
 sa==sb||error("frontier ordering differs")
 common=Dict{Symbol,Array{Float32,3}}()
 for key in (:temp,:dose,:food,:solid)
  left=Array(getproperty(a,key));right=Array(getproperty(b,key))
  reinterpret(UInt32,vec(left))==reinterpret(UInt32,vec(right))||error("common material field differs: $key")
  common[key]=left
 end
 ua=Array(first(a.levels).u);ub=Array(first(b.levels).u)
 law=electrical_weight_law(joinpath(context.source,"studies","electrothermal_archive.jl"))
 write(joinpath(output,"generated-propensity.jl"),law.generated_source)
 audit_a=electrical_cdf_audit(law,a.p,sa,da[1],da[3],da[4],c)
 audit_b=electrical_cdf_audit(law,b.p,sb,db[1],db[3],db[4],c)
 audit_a.site==Tuple(expected[85].baseline_site)&&audit_b.site==Tuple(expected[85].candidate_site)||error("actual event85 choices not reproduced")
 input=(step=84,next_event=85,seed=a.p.seed,eta=a.p.eta,deposit_per_cell=a.p.deposit_per_cell,
  occupied=copy(a.p.occupied),frontier=sa,baseline_potential=ua,candidate_potential=ub,
  common_fields=common,baseline_frontier_data=da,candidate_frontier_data=db,recipe=c)
 serialize(joinpath(output,"divergence-input.jls"),input)
 # Publish capture before refinement; interruption cannot lose the actual inputs.
 capture=(step=84,matched_prefix_choices=84,baseline=audit_a,candidate=audit_b,
  comparison=electrical_cdf_compare(audit_a,audit_b),source_archive_sha256=law.source_sha256,
  generated_propensity_sha256=law.generated_sha256)
 open(io->JSON3.pretty(io,capture),joinpath(output,"capture.json"),"w")
 references=Any[]
 for (label,initial) in (("baseline",ua),("candidate",ub))
  levels=TemporalPatch.device_levels(TemporalPatch.hierarchy(copy(initial),copy(a.p.occupied)),AMDGPU.ROCArray)
  workspace=TemporalPatch.DeviceMGNorm(first(levels))
  for tolerance in (1e-6,3e-7)
   stats=TemporalPatch.device_solve!(levels;tol=tolerance,maxcycles=500,workgroup=(32,4,2),norm_workspace=workspace)
   u=Array(first(levels).u);a.p.potential.=u
   independent=Interior.potential_residual(a.p)
   audit=electrical_cdf_audit(law,a.p,sa,Float32[u[I] for I in sa],da[3],da[4],c)
   push!(references,(initial_iterate=label,tolerance=tolerance,cycles=stats.cycles,
    residual=stats.residual,independent_residual=independent,audit=audit,
    baseline_comparison=electrical_cdf_compare(audit_a,audit),candidate_comparison=electrical_cdf_compare(audit_b,audit),
    baseline_potential_difference=maximum(abs.(ua.-u)),candidate_potential_difference=maximum(abs.(ub.-u))))
   serialize(joinpath(output,"reference-$(label)-$(tolerance).jls"),u)
  end
 end
 report=(scope="common geometry/material state84; actual event85 CDF margins and tighter Float32 harmonic references, not exact continuum truth",
  source=context.source,original_results=context.original_results,replay_choices=84,
  actual=capture,references=references)
 open(io->JSON3.pretty(io,report),joinpath(output,"results.json"),"w")
 println(JSON3.write((captured_step=84,original_baseline_choice=audit_a.site,original_candidate_choice=audit_b.site,
  tight_choices=[r.audit.site for r in references],original_cdf_distance=capture.comparison.cdf_maximum_difference)))
 report
end
function replay_electrical_divergence(source,original_results,output)
 source=abspath(source)
 original=read(joinpath(source,"scripts","study_electrical_temporal_coupling.jl"),String)
 needle=" for event in 1:events\n"
 length(findall(needle,original))==1||error("frozen driver capture marker changed")
 injected=replace(original,needle=>needle*"  event==85&&return capture_electrical_divergence!(a,b,c,output,rows)\n";count=1)
 ispath(output)&&error("refusing existing output");mkpath(output)
 write(joinpath(output,"instrumented-coupling.jl"),injected)
 write(joinpath(output,"driver-provenance.json"),JSON3.write((original_driver_sha256=bytes2hex(sha256(original)),
  instrumented_driver_sha256=bytes2hex(sha256(injected)),injection="return capture before event85; no preceding numerical edits")))
 # The original driver owns fresh-output creation; retain instrumentation adjacent.
 actual_output=joinpath(output,"capture")
 ELECTRICAL_DIVERGENCE_CONTEXT[]=(source=source,original_results=original_results)
 Base.include_string(Main,injected,joinpath(source,"scripts","study_electrical_temporal_coupling.jl"))
 entry=Base.invokelatest(getproperty,Main,:electrical_coupled_study)
 Base.invokelatest(entry,"amd",joinpath(source,"experiments","art","amber-confluence.toml"),actual_output;n=257,events=85)
end
if abspath(PROGRAM_FILE)==@__FILE__
 length(ARGS)==3||error("usage: FROZEN_SOURCE ORIGINAL_RESULTS FRESH_OUTPUT")
 replay_electrical_divergence(ARGS...)
end
