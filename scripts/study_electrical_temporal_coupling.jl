#!/usr/bin/env julia
# Independent fresh coupled histories; no production state is opened for writing.
using Interior, SHA, Statistics, Serialization, TOML
using Interior: JSON3, KernelAbstractions
const ELECTRICAL_COUPLED_PATHS=("studies/electrical_temporal_solver.jl","studies/electrical_patch_kernel.jl",
 "studies/electrical_local_correction.jl","studies/branch_solver_prototype.jl",
 "studies/branch_metal_multigrid.jl","studies/electrothermal_archive.jl",
 "scripts/study_electrical_temporal_coupling.jl")
const ELECTRICAL_COUPLED_BYTES=Dict(f=>read(joinpath(@__DIR__,"..",f)) for f in ELECTRICAL_COUPLED_PATHS)
module TemporalPatch
 include(joinpath(@__DIR__,"..","studies","electrical_temporal_solver.jl"))
end
module ArchiveLaw
 include(joinpath(@__DIR__,"..","studies","electrothermal_archive.jl"))
end
function electrical_coupled_study(backend,recipe,output;kwargs...)
 backend=="amd"&&(@eval using AMDGPU)
 Base.invokelatest(_electrical_coupled_study,backend,recipe,output;kwargs...)
end
function _electrical_coupled_study(backend_name,recipe,output;n=33,events=256)
 ispath(output)&&error("refusing existing output");mkpath(output)
 c=TOML.parsefile(recipe);tol=Float64(c["tolerance"]);dt=Float32(c["dt"])
 write(joinpath(output,"recipe.toml"),read(recipe))
 hardware=Dict{String,Any}()
 if backend_name=="amd"
  device=AMDGPU.device();props=AMDGPU.HIP.properties(device)
  hardware["device"]=AMDGPU.HIP.name(device);hardware["shared_memory_per_block"]=Int(props.sharedMemPerBlock)
  hardware["max_threads_per_block"]=Int(props.maxThreadsPerBlock)
  hardware["shared_memory_per_block"]>=36*1024&&hardware["max_threads_per_block"]>=256||error("patch hardware limit")
 end
 upload=backend_name=="cpu" ? identity : backend_name=="amd" ? AMDGPU.ROCArray : error("cpu|amd required")
 n=ArchiveLaw.branch_grid_size(n);group=Tuple(Int.(get(c,"solver_workgroup",[32,4,2])))
 function arm()
  mask=ArchiveLaw.root_mask(n;root_scale=get(c,"root_scale",1.),stub_scale=get(c,"stub_scale",1.))
  p=ArchiveLaw.process(n,mask);p.seed=UInt64(c["seed"]);p.eta=c["eta"];p.deposit_per_cell=.25f0
  levels=TemporalPatch.device_levels(TemporalPatch.hierarchy(p.potential,p.occupied),upload)
  solver=TemporalPatch.ElectricalTemporalSolver(levels;tol=.8tol,workgroup=group)
  TemporalPatch.electrical_temporal_solve!(solver)
  temp=upload(zeros(Float32,n,n,n));scratch=similar(temp);dose=upload(zeros(Float32,n,n,n))
  food=upload(ones(Float32,n,n,n));solid=upload(zeros(Float32,n,n,n));b=KernelAbstractions.get_backend(temp)
  for I in findall(mask);ArchiveLaw.archive_deposit!(b,1)(food,solid,LinearIndices(mask)[I],.25f0;ndrange=1);end
  (p=p,levels=levels,solver=solver,temp=temp,scratch=scratch,dose=dose,food=food,solid=solid)
 end
 a=arm();b=arm();h=ntuple(_->2f0/n,3)
 dt*(Float32(c["cooling"])+6Float32(c["diffusion"])/h[1]^2)<.9 || error("thermal positivity guard")
 function step(state,patched)
  sites,data=ArchiveLaw.frontier_data(state.p,first(state.levels),state.temp,state.dose,state.food,upload)
  choice=ArchiveLaw.choose_archive(state.p,sites,data,c;allow_arrest=true)
  choice===nothing&&return (state=state,site=nothing,cycles=0,residual=state.solver.prior_bound,patch_used=false)
  I,_=choice;ArchiveLaw.Interior.electrical_mark!(state.p,I);state.p.events+=1
  if patched;TemporalPatch.electrical_temporal_mark!(state.solver,I)
  else;TemporalPatch.device_mark!(state.levels,I)
  end
  backend=KernelAbstractions.get_backend(state.temp)
  ArchiveLaw.archive_deposit!(backend,1)(state.food,state.solid,LinearIndices(state.p.occupied)[I],.25f0;ndrange=1)
  position=ntuple(d->Float32((I[d]-.5)*2/n-1),3)
  ArchiveLaw.archive_heat!(backend,256)(state.scratch,state.temp,state.dose,state.food,state.solid,h,dt,position,
   Float32(c["heat"]),Float32(c["width"]),Float32(c["diffusion"]),Float32(c["cooling"]),
   Float32(c["threshold"]),Float32(c["conversion"]);ndrange=(n,n,n))
  state=merge(state,(temp=state.scratch,scratch=state.temp))
  stats=patched ? TemporalPatch.electrical_temporal_solve!(state.solver) :
   TemporalPatch.device_solve!(state.levels;tol=.8tol,workgroup=group,norm_workspace=state.solver.norm)
  (state=state,site=Tuple(I),cycles=stats.cycles,residual=stats.residual,patch_used=patched&&stats.patch_used)
 end
 function timed_step(state,patched)
  backend=KernelAbstractions.get_backend(state.temp);KernelAbstractions.synchronize(backend)
  start=time_ns();r=step(state,patched);KernelAbstractions.synchronize(backend)
  merge(r,(seconds=(time_ns()-start)*1e-9,))
 end
 rows=Any[];first_difference=nothing;max_residual=0.
 for event in 1:events
  ra,rb=if isodd(event);(timed_step(a,false),timed_step(b,true))
   else;right=timed_step(b,true);left=timed_step(a,false);(left,right)
  end
  a=ra.state;b=rb.state
  first_difference===nothing&&ra.site!=rb.site&&(first_difference=event)
  if event%16==0||event==events||ra.site===nothing||rb.site===nothing
   for state in (a,b)
    state.p.potential.=Array(first(state.levels).u)
    residual=Interior.potential_residual(state.p);max_residual=max(max_residual,residual)
    residual<=tol||error("independent global residual gate failed")
   end
  end
  push!(rows,(event=event,baseline_site=ra.site,candidate_site=rb.site,
   baseline_cycles=ra.cycles,candidate_cycles=rb.cycles,baseline_seconds=ra.seconds,candidate_seconds=rb.seconds,
   baseline_residual=ra.residual,candidate_residual=rb.residual,patch_used=rb.patch_used))
  event%32==0&&(println(JSON3.write(last(rows)));flush(stdout))
  (ra.site===nothing||rb.site===nothing)&&break
 end
 function state_record(state)
  fields=(potential=Array(first(state.levels).u),temperature=Array(state.temp),dose=Array(state.dose),
   food=Array(state.food),solid=Array(state.solid),occupied=copy(state.p.occupied))
  mass=sum(Float64,fields.food)+sum(Float64,fields.solid)
  minimum(fields.food)>=0&&abs(mass-n^3)/n^3<1e-6||error("material ledger failed")
  (fields=fields,events=state.p.events,mass_relative_error=(mass-n^3)/n^3,
   minimum_food=minimum(fields.food),solid_total=sum(Float64,fields.solid),dose_maximum=maximum(fields.dose))
 end
 fa=state_record(a);fb=state_record(b)
 serialize(joinpath(output,"final-states.jls"),(baseline=fa,candidate=fb))
 sources=Dict{String,String}()
 for (f,bytes) in ELECTRICAL_COUPLED_BYTES
  sources[f]=bytes2hex(sha256(bytes));write(joinpath(output,basename(f)),bytes)
 end
 field_differences=Dict(string(key)=>maximum(abs(Float64(x)-Float64(y)) for (x,y) in zip(getfield(fa.fields,key),getfield(fb.fields,key)))
  for key in (:potential,:temperature,:dose,:food,:solid))
 steady=length(rows)>16 ? rows[17:end] : rows
 report=(scope="fresh independent coupled histories; local timing under contention is not a GPU speed claim",
  backend=backend_name,hardware=hardware,julia=string(VERSION),threads=Threads.nthreads(),n=n,requested_events=events,completed_events=length(rows),recipe=c,
  coarse_fixture="actual recipe coefficients and root envelope at overridden resolution; no physical rescaling",
  warmup_events_excluded_from_steady_timing=min(16,length(rows)),
  steady_baseline_seconds=sum(r.baseline_seconds for r in steady),steady_candidate_seconds=sum(r.candidate_seconds for r in steady),
  final_field_differences=field_differences,
  first_choice_difference=first_difference,matched_choices=count(r->r.baseline_site==r.candidate_site,rows),
  occupied_jaccard=count(fa.fields.occupied .& fb.fields.occupied)/count(fa.fields.occupied .| fb.fields.occupied),
  maximum_independent_residual=max_residual,baseline_mass_relative_error=fa.mass_relative_error,
  candidate_mass_relative_error=fb.mass_relative_error,baseline_solid_total=fa.solid_total,candidate_solid_total=fb.solid_total,
  baseline_cycles=sum(r.baseline_cycles for r in rows),candidate_cycles=sum(r.candidate_cycles for r in rows),
  source_sha256=sources,rows=rows)
 open(io->JSON3.pretty(io,report),joinpath(output,"results.json"),"w")
 report
end
if abspath(PROGRAM_FILE)==@__FILE__
 length(ARGS)>=3||error("usage: cpu|amd RECIPE FRESH_OUTPUT [N=33] [EVENTS=256]")
 electrical_coupled_study(ARGS[1],ARGS[2],ARGS[3];n=length(ARGS)>3 ? parse(Int,ARGS[4]) : 33,
  events=length(ARGS)>4 ? parse(Int,ARGS[5]) : 256)
end
