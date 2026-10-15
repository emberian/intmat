# Separate device prototype/driver. Never loaded by the existing CPU run.
include(joinpath(@__DIR__,"branch_solver_prototype.jl"))
include(joinpath(@__DIR__,"branch_metal_multigrid.jl"))
using Serialization, SHA, TOML
function device_frontier!(p,level,upload)
 sites=findall(p.frontier);indices=Int32[LinearIndices(p.potential)[I] for I in sites]
 gpu_indices=upload(indices);gpu_values=similar(level.u,Float32,length(sites))
 backend=KernelAbstractions.get_backend(level.u)
 branch_gather_kernel!(backend,256)(gpu_values,level.u,gpu_indices;ndrange=length(sites))
 values=Array(gpu_values)
 for (I,value) in zip(sites,values);p.potential[I]=value;end
 sites
end
function metal_benchmark(ns;upload)
 results=[]
 for n in ns
  mask=specimen(n);p=process(n,mask);cpu=hierarchy(p.potential,p.occupied);device=device_levels(cpu,upload)
  # Warm compilation separately from timing.
  device_solve!(device;tol=1e-6);copyto!(first(device).u,process(n,mask).potential)
  time=@elapsed stat=device_solve!(device;tol=1e-6)
  u=Array(first(device).u);p.potential.=u;core=Interior.potential_residual(p)
  comparison=nothing
  if n<=65
   p.potential.=0;p.potential[p.occupied].=1;mgsolve!(cpu;tol=1e-6)
   comparison=maximum(abs.(p.potential.-u))
  end
  J=CartesianIndex((n+1)÷2+1,(n+1)÷2,(n+1)÷2);device_mark!(device,J)
  warm=@elapsed ws=device_solve!(device;tol=1e-6)
  row=(n=n,cold_seconds=time,cycles=stat.cycles,core_residual=core,max_cpu_difference=comparison,warm_seconds=warm,warm_cycles=ws.cycles,warm_residual=ws.residual)
  push!(results,row);println(JSON3.write(row));flush(stdout)
 end
 results
end
function metal_growth_parity(;n=33,events=100,upload)
 mask=falses(n,n,n);mask[ntuple(_->(n+1)÷2,3)...]=true
 p=process(n,mask);q=process(n,mask);p.tolerance=q.tolerance=3e-7;p.eta=q.eta=4.
 cpu=hierarchy(p.potential,p.occupied);device=device_levels(hierarchy(q.potential,q.occupied),upload)
 mgsolve!(cpu;tol=3e-7,maxcycles=500);device_solve!(device;tol=3e-7,maxcycles=500)
 matches=0;maxdiff=0.;elapsed=@elapsed for event in 1:events
  device_frontier!(q,first(device),upload)
  I=select_site(p,(1.,1.,1.));J=select_site(q,(1.,1.,1.))
  I==J || error("device growth diverged at$event: $I != $J")
  matches+=1;Interior.electrical_mark!(p,I);p.events+=1;mark_hierarchy!(cpu,I);mgsolve!(cpu;tol=3e-7,maxcycles=500)
  Interior.electrical_mark!(q,J);q.events+=1;device_mark!(device,J);device_solve!(device;tol=3e-7,maxcycles=500)
  maxdiff=max(maxdiff,maximum(abs.(p.potential.-Array(first(device).u))))
  p.reached_ground && break
 end
 row=(n=n,matched_events=matches,max_potential_difference=maxdiff,seconds=elapsed,grounded=p.reached_ground)
 println(JSON3.write(row));row
end
function branch_geometry_snapshot(path,p)
 sites=findall(p.occupied);ids=Dict(I=>i for (i,I) in enumerate(sites));edges=Tuple{Int,Int}[]
 for (i,I) in enumerate(sites),axis in 1:3
  J=I+CartesianIndex(ntuple(d->d==axis ? 1 : 0,3))
  haskey(ids,J) && push!(edges,(i,ids[J]))
 end
 degree=zeros(Int,length(sites))
 for (a,b) in edges;degree[a]+=1;degree[b]+=1;end
 report=(event=p.events,shape=size(p.potential),positions=Tuple.(sites),edges=edges,tips=count(==(1),degree),branch_vertices=count(>=(3),degree),bounds=ntuple(d->(minimum(I->I[d],sites),maximum(I->I[d],sites)),3))
 open(io->JSON3.pretty(io,report),path,"w")
end
function save_branch_checkpoint(path,p,level,history,cfg)
 p.potential.=Array(level.u)
 # Uniform untouched reservoir is reconstructed exactly from occupied mask and .8 transfer.
 temporary=path*".new"
 open(temporary,"w") do io;serialize(io,(process=p,history=history,recipe=cfg));end
 mv(temporary,path;force=true)
end
include(joinpath(@__DIR__,"branch_export_support.jl"))
function metal_native_branch(recipe;upload,resume=nothing,batch_events=0)
 cfg=TOML.parsefile(recipe);cfg["slug"]*= "-metal";cfg["backend"]="Metal"
 n=branch_grid_size(cfg["n"]);batch_events>=0 || error("batch event count must be nonnegative")
 out=joinpath(dirname(@__DIR__),"outputs","refinement",cfg["slug"])
 isfile(joinpath(out,"run.json")) && error("completed specimen already exists")
 out_exists=ispath(out)
 if out_exists
  resume===nothing && error("partial output already exists; supply its checkpoint explicitly or choose a fresh slug")
  own_checkpoint=joinpath(out,"checkpoint.jls")
  isfile(own_checkpoint) && isfile(resume) && realpath(resume)==realpath(own_checkpoint) || error("existing partial output requires its own checkpoint; use a fresh slug for a different checkpoint")
 end
 history=Any[]
 if resume!==nothing
  saved=open(deserialize,resume)
  for key in ("n","seed","eta","tolerance","extent")
   saved.recipe[key]==cfg[key] || error("checkpoint recipe mismatch: $key")
  end
  p=saved.process;history=saved.history
  size(p.potential)==(n,n,n) || error("checkpoint shape differs from recipe")
  p.events<=cfg["events"] || error("requested horizon precedes checkpoint event")
  p.potential[p.occupied].=1
 else
  mask=falses(n,n,n);mask[ntuple(_->(n+1)÷2,3)...]=true
  p=process(n,mask);p.seed=UInt64(cfg["seed"]);p.eta=cfg["eta"];p.tolerance=cfg["tolerance"]
 end
 mkpath(dirname(out));out_exists || mkdir(out)
 open(io->TOML.print(io,cfg),joinpath(out,"recipe.toml"),"w")
 levels=device_levels(hierarchy(p.potential,p.occupied),upload)
 device_solve!(levels;tol=.9p.tolerance)
 stop_at=batch_events>0 ? min(cfg["events"],p.events+batch_events) : cfg["events"]
 elapsed=@elapsed while p.events<stop_at && !p.reached_ground
  device_frontier!(p,first(levels),upload)
  I=select_site(p,Tuple(Float64.(cfg["extent"])./n))
  Interior.electrical_mark!(p,I);p.events+=1;device_mark!(levels,I)
  t=@elapsed solved=device_solve!(levels;tol=.9p.tolerance)
  p.last_residual=solved.residual;p.last_iterations=solved.cycles;p.elapsed=p.events*p.event_interval
  if p.events%25==0 || p.reached_ground || p.events==stop_at
   p.potential.=Array(first(levels).u);core=Interior.potential_residual(p)
   core<=p.tolerance || error("independent core residual $core exceeds tolerance")
   p.last_residual=core
   sites=findall(p.occupied); span=ntuple(d->maximum(I->I[d],sites)-minimum(I->I[d],sites)+1,3)
   row=(event=p.events,time=p.elapsed,occupied=count(p.occupied),span_cells=span,residual=core,vcycles=solved.cycles,solve_seconds=t,reached_ground=p.reached_ground)
   push!(history,row);println(JSON3.write(row));flush(stdout)
   if p.events%200==0 || p.events==stop_at || p.reached_ground
    branch_geometry_snapshot(joinpath(out,"geometry-"*lpad(string(p.events),4,'0')*".json"),p)
   end
   save_branch_checkpoint(joinpath(out,"checkpoint.jls"),p,first(levels),history,cfg)
   open(io->JSON3.pretty(io,Dict("shape"=>[n,n,n],"history"=>history,"recipe"=>cfg)),joinpath(out,"progress.json"),"w")
  end
 end
 if p.events<cfg["events"] && !p.reached_ground
  save_branch_checkpoint(joinpath(out,"checkpoint.jls"),p,first(levels),history,cfg)
  println("PAUSED ",out," event=",p.events);return
 end
 p.potential.=Array(first(levels).u);p.last_residual=Interior.potential_residual(p)
 w=World(Grid(n,n,n;extent=Tuple(cfg["extent"]));seed=cfg["seed"])
 w.time=p.elapsed;w.fields.deposit[p.occupied].=p.deposit_per_cell;w.fields.nutrient[p.occupied].-=p.deposit_per_cell
 export_native_branch(w,p,cfg,out,history,elapsed)
end
function branch_metal_cli(args=ARGS)
 usage="usage: study refinement_branches_metal grow RECIPE.toml [CHECKPOINT.jls] | benchmark [N ...]"
 isempty(args) && error(usage)
 mode=first(args)
 if mode=="grow"
  length(args) in (2,3) || error(usage)
  branch_grid_size(TOML.parsefile(args[2])["n"])
 elseif mode=="benchmark"
  ns=length(args)>1 ? parse.(Int,args[2:end]) : [33,65,257]
  foreach(branch_grid_size,ns)
 else
  error(usage)
 end
 @eval using Metal
 Metal.functional() || error("Metal unavailable; run the CPU study or use a working Metal environment")
 upload=Metal.MtlArray
 if mode=="grow"
  metal_native_branch(args[2];upload,resume=length(args)==3 ? args[3] : nothing,batch_events=parse(Int,get(ENV,"BRANCH_BATCH_EVENTS","0")))
 else
  report=Dict("benchmarks"=>metal_benchmark(ns;upload),"growth_parity"=>metal_growth_parity(;upload))
  mkpath(joinpath(@__DIR__,"..","outputs","refinement"))
  open(io->JSON3.pretty(io,report),joinpath(@__DIR__,"..","outputs","refinement","branch-metal-benchmark.json"),"w")
 end
end
if abspath(PROGRAM_FILE)==(@__FILE__) || (isdefined(Main,:INTERIOR_RUNNING_STUDY) && Main.INTERIOR_RUNNING_STUDY=="refinement_branches_metal")
 branch_metal_cli()
end
