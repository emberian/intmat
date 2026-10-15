#!/usr/bin/env julia
# New closed, dense initial condition, not a migration of a retained trajectory.
using Interior, TOML, SHA, Serialization, Statistics
using Interior: JSON3
const RKC_ECOLOGY_FILES=("cellular_implicit_mechanics.jl","cellular_pair_proximal.jl","cellular_rkc_mechanics.jl","cellular_rkc_funding.jl","cellular_rkc_ecology.jl")
for file in RKC_ECOLOGY_FILES;Base.include(Interior,joinpath(@__DIR__,"..","studies",file));end
# Refine only mechanics: physiology and its remaining-outer-interval funding
# horizon still run once per common outer step.
function install_fine_mechanics()
 source=read(joinpath(@__DIR__,"..","src","cellular_contacts.jl"),String)
 bytes2hex(sha256(source))==Interior.arena_sources()["cellular_contacts.jl"]||error("loaded source mismatch")
 lo=findfirst("function cellular_mechanics!(",source)
 hi=findnext("\nend\n",source,last(lo))
 body=source[first(lo):last(hi)]
 length(findall("dt::Real;",body))==1&&length(findall("h=min(remaining,speed>0",body))==1||error("reference source markers changed")
 body=replace(body,"function cellular_mechanics!("=>"function _rkc_fine_reference!(",
   "dt::Real;"=>"dt::Real,refinement::Integer;",
   "h=min(remaining,speed>0"=>"h=min(Float64(dt)/refinement,remaining,speed>0")
 occursin("h=min(Float64(dt)/refinement",body)||error("reference cap marker changed")
 Base.include_string(Interior,body,"rkc-fine-reference-study")
 body
end
const RKC_REFERENCE_SOURCE=install_fine_mechanics()
struct RKCReferencePolicy{N} end
(p::RKCReferencePolicy{N})(grid,ids,positions,radii,bonds,dt;kwargs...) where N =
 Interior._rkc_fine_reference!(grid,ids,positions,radii,bonds,dt,N;kwargs...)
function rkc_dense_initial(population)
 cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","arena","constructive-volume.toml"))
 cfg["grid"]["shape"]=[10,10,10];cfg["grid"]["extent"]=[5.,5.,5.]
 cfg["cellular"]["population"]=population;cfg["cellular"]["capacity"]=4population
 for name in ("resources","feedstock")
  cfg[name]["source_rates"].*=0
  haskey(cfg[name],"sink_rates")&&(cfg[name]["sink_rates"].*=0)
 end
 build_arena(cfg)
end
function rkc_full_study(out;population=192,steps=10)
 ispath(out)&&error("fresh output required");mkpath(out)
 source=Dict(file=>read(joinpath(@__DIR__,"..","studies",file),String) for file in RKC_ECOLOGY_FILES)
 serialize(joinpath(out,"provenance.jls"),(core=Interior.arena_sources(),study=source,reference=RKC_REFERENCE_SOURCE,driver=read(@__FILE__,String)))
 initial=rkc_dense_initial(population)
 rows=Any[];states=Dict{String,Any}()
 # All arms share the same seeded biological state and outer clock. End-to-end
 # timing includes every physiology/transport stage; deepcopy, diagnostics and
 # checkpoint observation are outside the timer.
 arms=(("adaptive",Interior.cellular_mechanics!), ("rkc16",Interior.RKCMechanicalPolicy()),
  ("fine256",RKCReferencePolicy{256}()),("fine512",RKCReferencePolicy{512}()),("fine1024",RKCReferencePolicy{1024}()),
  ("rkc8",Interior.RKCMechanicalPolicy(splits=8)),("rkc32",Interior.RKCMechanicalPolicy(splits=32)))
 for (name,policy) in arms
  a=deepcopy(initial);physical=a.processes
  a.processes=map(p->p isa ConstructiveEcology ? Interior.RKCExperimentalEcology(p;policy) : p,physical)
  model=only(p for p in physical if p isa ConstructiveEcology)
  mass=Interior.arena_inventory(a.world,physical);energy=Interior.arena_energy(physical)
  # Warm code on a disposable state, not the measured biological history.
  warm=deepcopy(a);arena_step!(warm);GC.gc()
  timing=@timed for _ in 1:steps;arena_step!(a);end
  state=(ids=copy(model.cells.ids),positions=copy(model.cells.positions),energy=copy(model.cells.energy),
   hidden=copy(model.cells.hidden),metabolites=copy(model.cells.metabolites),biomass=copy(model.cells.biomass),
   resource_fields=deepcopy(model.cells.resources.fields),enzymes=deepcopy(model.habitat.enzymes),
   bonds=deepcopy(model.cells.bonds),work=model.cells.work_spent,births=model.cells.births,deaths=model.cells.deaths)
  serialize(joinpath(out,name*".jls"),state);states[name]=state
  row=(name=name,seconds=timing.time,bytes=timing.bytes,gc_seconds=timing.gctime,steps=steps,time=a.world.time,population=length(state.ids),
   births=state.births,deaths=state.deaths,work=state.work,
   mass_residual=Interior.arena_inventory(a.world,physical)-mass,energy_residual=Interior.arena_energy(physical)-energy,
   finite=all(isfinite,state.positions)&&all(isfinite,state.energy),minimum_funds=minimum(state.energy),
   policy_calls=policy isa Interior.RKCMechanicalPolicy ? policy.calls : 0)
  push!(rows,row);println(JSON3.write(row));flush(stdout)
 end
 reference=states["fine1024"];errors=Dict{String,Any}()
 for name in ("adaptive","rkc8","rkc16","rkc32","fine256","fine512")
  state=states[name];common=intersect(state.ids,reference.ids);index=Dict(id=>i for(i,id)in enumerate(state.ids));refindex=Dict(id=>i for(i,id)in enumerate(reference.ids))
  sum2=0.;largest=0.
  for id in common
   d2=sum((mod(state.positions[d,index[id]]-reference.positions[d,refindex[id]]+2.5,5.)-2.5)^2 for d in 1:3)
   sum2+=d2;largest=max(largest,sqrt(d2))
  end
  e=Set(state.bonds.endpoints);r=Set(reference.bonds.endpoints)
  errors[name]=(common_ids=length(common),same_ids=state.ids==reference.ids,rms=sqrt(sum2/(3length(common))),maximum=largest,
   bond_jaccard=isempty(union(e,r)) ? 1. : length(intersect(e,r))/length(union(e,r)),work_error=state.work-reference.work,
   resource_relative_l2=sqrt(sum(sum(abs2,Float64.(state.resource_fields[k]).-reference.resource_fields[k]) for k in eachindex(state.resource_fields))/sum(sum(abs2,a) for a in reference.resource_fields)))
 end
 timings=Any[]
 samples=parse(Int,get(ENV,"RKC_FULL_SAMPLES","0"))
 for sample in 1:samples
  pair=Dict{String,Any}()
  for arm in (isodd(sample) ? (1,2) : (2,1))
   name,original_policy=arms[arm];policy=deepcopy(original_policy)
   a=deepcopy(initial);a.processes=map(p->p isa ConstructiveEcology ? Interior.RKCExperimentalEcology(p;policy) : p,a.processes)
   GC.gc();timing=@timed for _ in 1:steps;arena_step!(a);end
   pair[name]=(seconds=timing.time,bytes=timing.bytes,gc_seconds=timing.gctime)
  end
  pair["sample"]=sample;pair["paired_speedup"]=pair["adaptive"].seconds/pair["rkc16"].seconds
  push!(timings,pair);println(JSON3.write(pair));flush(stdout)
 end
 result=(fixture="fresh seeded $(population)-founder dense closed constructive ecology; no old-state migration",config=initial.config,
  rows=rows,errors=errors,timings=timings,median_paired_speedup=isempty(timings) ? nothing : median(t["paired_speedup"] for t in timings),source=Interior.arena_sources(),study_sha256=Dict(f=>bytes2hex(sha256(s)) for(f,s)in source))
 open(io->JSON3.pretty(io,result),joinpath(out,"results.json"),"w");println(JSON3.write(errors));result
end
if abspath(PROGRAM_FILE)==@__FILE__;rkc_full_study(only(ARGS));end
