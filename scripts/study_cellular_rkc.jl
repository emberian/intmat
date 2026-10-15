#!/usr/bin/env julia
# Mechanics-only numerical experiment, never edits an arena/checkpoint.
using Interior, Serialization, SHA, Statistics
using Interior: JSON3
const RKC_STUDY_FILES=("studies/cellular_implicit_mechanics.jl","studies/cellular_pair_proximal.jl","studies/cellular_rkc_mechanics.jl","studies/cellular_rkc_funding.jl")
const RKC_STUDY_BYTES=Dict(f=>read(joinpath(@__DIR__,"..",f)) for f in RKC_STUDY_FILES)
for f in RKC_STUDY_FILES;Base.include_string(Interior,String(copy(RKC_STUDY_BYTES[f])),joinpath(@__DIR__,"..",f));end
function rkc_fixture_study(fixture_path,reference_path,old_adaptive_path,out)
 ispath(out)&&error("refusing existing output");mkpath(out)
 for(f,bytes)in RKC_STUDY_BYTES;write(joinpath(out,basename(f)),bytes);end
 write(joinpath(out,"driver.jl"),read(@__FILE__))
 fixture=deserialize(fixture_path);reference=deserialize(reference_path);old=deserialize(old_adaptive_path)
 funding_mode=get(ENV,"RKC_FUNDING","original");funding_mode in ("original","cached")||error("unknown funding mode")
 funded=funding_mode=="cached" ? Interior.cellular_rkc_funded_mechanics! : Interior.cellular_implicit_mechanics!
 grid=fixture.grid;ids=fixture.ids;radii=fixture.radii;q=fixture.parameters;duration=parse(Float64,string(fixture.dt))
 kwargs=(motility=fixture.motility,adhesion=fixture.adhesion,repulsion=q.repulsion,stiffness=q.stiffness,
  break_strain=q.break_strain,max_degree=q.max_degree,motor_cost=q.motor_cost)
 function difference(x,y)
  nd=grid.size[3]==1 ? 2 : 3;sum2=0.;largest=0.
  for i in axes(x,2)
   d2=0.
   for d in 1:nd
    L=Float64(grid.extent[d]);delta=mod(x[d,i]-y[d,i]+L/2,L)-L/2
    d2+=delta^2
   end
   sum2+=d2;largest=max(largest,sqrt(d2))
  end
  (rms=sqrt(sum2/(nd*length(ids))),maximum=largest)
 end
 x=copy(fixture.positions);b=deepcopy(fixture.bonds);energy=copy(fixture.energy)
 baseline=Interior.cellular_mechanics!(grid,ids,x,radii,b,duration;energy,kwargs...)
 exact_old_default=isequal(x,old.x)&&isequal(energy,old.energy)&&isequal(b.endpoints,old.b.endpoints)&&isequal(b.rest_lengths,old.b.rest_lengths)
 baseline_comparison=(exact=exact_old_default,position_difference=difference(x,old.x),energy_maximum_difference=maximum(abs.(energy-old.energy)),
  same_bonds=b.endpoints==old.b.endpoints&&b.rest_lengths==old.b.rest_lengths,substeps=baseline.substeps,old_substeps=old.substeps,
  duration=duration,fixture_dt_type=string(typeof(fixture.dt)),machine=Sys.MACHINE)
 println(JSON3.write(baseline_comparison));flush(stdout)
 open(io->JSON3.pretty(io,baseline_comparison),joinpath(out,"baseline-comparison.json"),"w")
 (baseline_comparison.same_bonds&&baseline_comparison.position_difference.maximum<1e-10&&baseline_comparison.energy_maximum_difference<1e-10)||error("retained default differs materially; investigate before reference use")
 rows=Any[(name="adaptive",splits=1,substeps=baseline.substeps,position_difference=difference(x,reference.x),
  work=sum(baseline.spent_motor_work),energy_l1_difference=sum(abs.(energy-reference.energy)),ledger=baseline.energy_residual)]
 for splits in parse.(Int,split(get(ENV,"RKC_CASE_SPLITS","1,2,4,8,16"),","))
  x=copy(fixture.positions);b=deepcopy(fixture.bonds);energy=copy(fixture.energy);cache=Interior._CellularRKCCache()
  result=funded(grid,ids,x,radii,b,duration;energy,kwargs...,
   relaxer=cache,symmetric_motor=true,max_step=duration/splits,max_evaluations=1024)
  a=Set(b.endpoints);z=Set(reference.b.endpoints)
  row=(name="rkc2",splits=splits,substeps=result.substeps,retries=result.retries,force_evaluations=cache.force_evaluations,
   endpoint_energy_evaluations=cache.energy_evaluations,accepted_stages=copy(cache.accepted_stages),broadphase_builds=cache.graph.builds,
   position_difference=difference(x,reference.x),work=sum(result.spent_motor_work),energy_l1_difference=sum(abs.(energy-reference.energy)),
   bond_jaccard=isempty(union(a,z)) ? 1. : length(intersect(a,z))/length(union(a,z)),ledger=result.energy_residual)
  push!(rows,row);serialize(joinpath(out,"rkc-$(splits).jls"),(x=x,b=b,energy=energy,result=result))
  println(JSON3.write(row));flush(stdout)
 end
 samples=parse(Int,get(ENV,"RKC_BENCH_SAMPLES","0"));timings=Any[]
 function measured(method)
  x=copy(fixture.positions);b=deepcopy(fixture.bonds);energy=copy(fixture.energy);cache=Interior._CellularRKCCache()
  GC.gc()
  timing=@timed if method===:adaptive
   Interior.cellular_mechanics!(grid,ids,x,radii,b,duration;energy,kwargs...)
  else
   funded(grid,ids,x,radii,b,duration;energy,kwargs...,
    relaxer=cache,symmetric_motor=true,max_step=duration/16,max_evaluations=1024)
  end
  (seconds=timing.time,bytes=timing.bytes,gc_seconds=timing.gctime,substeps=timing.value.substeps)
 end
 if samples>0
  measured(:adaptive);measured(:rkc2)
  for sample in 1:samples
   a,b=if isodd(sample);(measured(:adaptive),measured(:rkc2))
    else;right=measured(:rkc2);left=measured(:adaptive);(left,right)
   end
   push!(timings,(sample=sample,adaptive=a,rkc16=b,paired_speedup=a.seconds/b.seconds))
  end
 end
 report=(scope="mechanics-only equal outer interval; local times not performance claims, controlled host timing explicitly separate",population=length(ids),duration=duration,funding_mode=funding_mode,reference_work=reference.work,timings=timings,
  median_paired_speedup=isempty(timings) ? nothing : median(r.paired_speedup for r in timings),
  baseline_comparison=baseline_comparison,exact_old_default=exact_old_default,fixture_sha256=bytes2hex(sha256(read(fixture_path))),
  reference_sha256=bytes2hex(sha256(read(reference_path))),old_adaptive_sha256=bytes2hex(sha256(read(old_adaptive_path))),
  loaded_default_sources=Interior.arena_sources(),study_sources=Dict(f=>bytes2hex(sha256(bytes)) for (f,bytes) in RKC_STUDY_BYTES),rows=rows)
 open(io->JSON3.pretty(io,report),joinpath(out,"results.json"),"w")
 report
end
if abspath(PROGRAM_FILE)==@__FILE__
 length(ARGS)==4||error("usage: FIXTURE REFERENCE_1024 OLD_ADAPTIVE FRESH_OUTPUT")
 rkc_fixture_study(ARGS...)
end
