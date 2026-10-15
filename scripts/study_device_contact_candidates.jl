#!/usr/bin/env julia
using Interior,KernelAbstractions,Random,SHA,Serialization,Statistics
using Interior: JSON3
Base.include(Interior,joinpath(@__DIR__,"..","studies","device_contact_candidates.jl"))
include("device_habitat_pilot_fixture.jl")
function contact_candidate_fixture(dense)
 w,p=device_habitat_pilot_fixture(CPU();size=12,founders=dense ? 128 : 64,capacity=512,
  support_strategy=:materialized,energy_type=Float64)
 if dense
  for i in eachindex(p.cells.state.active),d in 1:3
   p.cells.state.active[i]==1&&(p.cells.state.position[i,d]=2f0+.25f0*(p.cells.state.position[i,d]-2f0))
  end
 end
 w,p
end
function candidate_pair_error(reference,candidate)
 a=reference.cells.state;b=candidate.cells.state;ids=all(a.ids.==b.ids)&&all(a.active.==b.active)
 exact=isequal(reference,candidate)
 aa=Interior.device_habitat_audit(reference);bb=Interior.device_habitat_audit(candidate)
 float_errors=Dict{String,Float64}()
 for (label,x,y) in (("cells",reference.cells.state,candidate.cells.state),("environment",reference.cells.environment,candidate.cells.environment),
  ("habitat",reference.arrays,candidate.arrays),("feedstock",reference.feedstock,candidate.feedstock))
  for k in keys(x)
   if eltype(x[k])<:AbstractFloat
    float_errors[label*"/"*string(k)]=maximum(abs,Float64.(x[k]).-Float64.(y[k]);init=0.)
   end
  end
 end
 (exact=exact,same_ids=ids,same_bonds=isequal(a.bonds,b.bonds)&&isequal(a.bond_ids,b.bond_ids),
  work_difference=bb.work-aa.work,material_difference=bb.material-aa.material,potential_difference=bb.potential-aa.potential,
  max_float_errors=float_errors)
end
function contact_candidate_study(out;steps=6,samples=5)
 ispath(out)&&error("fresh output required");mkpath(out)
 study=read(joinpath(@__DIR__,"..","studies","device_contact_candidates.jl"),String)
 serialize(joinpath(out,"provenance.jls"),(core=Interior.arena_sources(),study=study,driver=read(@__FILE__,String),fixture=read(joinpath(@__DIR__,"device_habitat_pilot_fixture.jl"),String)))
 records=Any[]
 for dense in (false,true)
  base=contact_candidate_fixture(dense)
  first=Interior.device_habitat_audit(Interior.device_habitat_snapshot(base[2]))
  # Warm both complete coupled schedules on disposable states.
  for mode in (:ordinary,:coverage)
   w,p=deepcopy(base);cache=Interior.device_contact_candidate_cache(p.cells;strict_order=false)
   carrier=Interior.DeviceContactCandidateCarrier(p.cells,cache)
   mode===:ordinary ? Interior.device_cellular_step!(w,p.cells,.001;anabolism=p.habitat) : Interior.device_candidate_step!(w,carrier,.001;anabolism=p.habitat)
  end
  pairs=Any[]
  for sample in 1:samples
   rows=Dict{String,Any}();snapshots=Dict{Symbol,Any}()
   for mode in (isodd(sample) ? (:ordinary,:coverage) : (:coverage,:ordinary))
    w,p=deepcopy(base);cache=Interior.device_contact_candidate_cache(p.cells;strict_order=false)
    carrier=Interior.DeviceContactCandidateCarrier(p.cells,cache)
    GC.gc()
    timing=@timed for _ in 1:steps
     mode===:ordinary ? Interior.device_cellular_step!(w,p.cells,.001;anabolism=p.habitat) : Interior.device_candidate_step!(w,carrier,.001;anabolism=p.habitat)
    end
    snapshot=Interior.device_habitat_snapshot(p);snapshots[mode]=snapshot
    audit=Interior.device_habitat_audit(snapshot)
    rows[string(mode)]=(seconds=timing.time,bytes=timing.bytes,gc_seconds=timing.gctime,audit=audit,
     mass_residual=audit.material-first.material-(audit.supplied-first.supplied)+(audit.removed-first.removed),
     energy_residual=audit.potential-first.potential+audit.work-first.work+audit.heat-first.heat-
      (audit.supplied_potential-first.supplied_potential)+(audit.removed_potential-first.removed_potential),
     final_macro_cache_hits=sum(cache.hits),final_macro_fallbacks=sum(cache.fallbacks),final_macro_invalid=only(cache.invalid),
     final_macro_overflow_rows=count(>(0),cache.overflow),maximum_cached_candidates=maximum(cache.count),skin=cache.skin)
   end
   errors=candidate_pair_error(snapshots[:ordinary],snapshots[:coverage])
   rows["errors"]=errors;rows["sample"]=sample;rows["paired_speedup"]=rows["ordinary"].seconds/rows["coverage"].seconds
   push!(pairs,rows);println(JSON3.write((dense=dense,pair=rows)));flush(stdout)
   sample==1&&serialize(joinpath(out,dense ? "dense-states.jls" : "sparse-states.jls"),snapshots)
  end
  push!(records,(dense=dense,founders=dense ? 128 : 64,capacity=512,grid=12,steps=steps,dt=.001,
   paired_median_speedup=median(p["paired_speedup"] for p in pairs),samples=pairs))
 end
 result=(scope="CPU KernelAbstractions complete coupled cellular/catalytic/sorptive schedule; no GPU throughput claim",julia=string(VERSION),machine=Sys.MACHINE,
  core=Interior.arena_sources(),study_sha256=bytes2hex(sha256(study)),records=records)
 open(io->JSON3.pretty(io,result),joinpath(out,"results.json"),"w");result
end
if abspath(PROGRAM_FILE)==@__FILE__;contact_candidate_study(only(ARGS));end
