include(joinpath(@__DIR__,"branch_solver_prototype.jl"))
using TOML, SHA, Printf, Serialization
function fine_branch(recipe;pilot_events=0)
 cfg=TOML.parsefile(recipe);n=branch_grid_size(cfg["n"]);pilot_events>=0 || error("pilot event count must be nonnegative")
 out=joinpath(dirname(@__DIR__),"outputs","refinement",cfg["slug"])
 ispath(out) && error("output already exists; choose a fresh recipe slug (CPU study does not resume partial outputs)")
 mkpath(dirname(out));mkdir(out);cp(recipe,joinpath(out,"recipe.toml"))
 mask=falses(n,n,n);root=CartesianIndex(ntuple(_->(n+1)÷2,3));mask[root]=true
 p=process(n,mask);p.seed=UInt64(cfg["seed"]);p.eta=cfg["eta"];p.tolerance=cfg["tolerance"]
 w=World(Grid(n,n,n;extent=Tuple(cfg["extent"]));seed=cfg["seed"])
 w.fields.nutrient[root]-=p.deposit_per_cell;w.fields.deposit[root]+=p.deposit_per_cell
 levels=hierarchy(p.potential,p.occupied);mgsolve!(levels;tol=p.tolerance)
 history=Any[];total=pilot_events>0 ? pilot_events : cfg["events"]
 elapsed=@elapsed for event in 1:total
  site=select_site(p,w.grid.spacing)
  w.fields.nutrient[site]>=p.deposit_per_cell || error("local nutrient budget")
  Interior.electrical_occupy!(w,p,site);p.events+=1;mark_hierarchy!(levels,site)
  t=@elapsed solved=mgsolve!(levels;tol=p.tolerance)
  # Use the independent production residual formula for the stopping contract.
  residual=Interior.potential_residual(p,w.grid.spacing)
  while residual>p.tolerance
   vcycle!(levels);residual=Interior.potential_residual(p,w.grid.spacing)
  end
  p.last_residual=residual;p.last_iterations=solved.cycles;p.elapsed=p.events*p.event_interval;w.time=p.elapsed
  if event%25==0 || p.reached_ground || event==total
   row=(event=event,time=w.time,occupied=count(p.occupied),residual=residual,vcycles=solved.cycles,solve_seconds=t,reached_ground=p.reached_ground)
   push!(history,row);println(JSON3.write(row));flush(stdout)
   open(joinpath(out,"progress.json"),"w") do io;JSON3.pretty(io,Dict("shape"=>[n,n,n],"history"=>history,"recipe"=>cfg));end
  end
  p.reached_ground && break
 end
 open(joinpath(out,"checkpoint.jls"),"w") do io;serialize(io,(process=p,deposit=w.fields.deposit,nutrient=w.fields.nutrient,time=w.time));end
 initial_mass=Float64(n)^3;final_mass=sum(Float64,w.fields.nutrient)+sum(Float64,w.fields.deposit)
 final_mass==initial_mass || error("mass conservation")
 record_state(w,(p,),joinpath(out,"state");step=p.events)
 sites=findall(p.occupied);ids=Dict(I=>i for (i,I) in enumerate(sites));edges=0
 open(joinpath(out,"channels.obj"),"w") do io
  println(io,"# Fresh native fine-grid conductor, exact occupied face adjacency")
  for I in sites;println(io,"v ",join(ntuple(d->2*(I[d]-.5)/n-1,3)," "));end
  for (i,I) in enumerate(sites),axis in 1:3
   J=I+CartesianIndex(ntuple(d->d==axis ? 1 : 0,3))
   if haskey(ids,J);println(io,"l ",i," ",ids[J]);edges+=1;end
  end
 end
 # Bounded optical fields; no optical dilation or interpolation is used for the graph.
 fields=(density=fill(.02f0,n,n,n),activity=.02f0.*w.fields.deposit,structure=w.fields.deposit,mask=ones(Float32,n,n,n),potential=p.potential)
 export_volume(fields,joinpath(out,"volume"))
 degrees=[count(J->haskey(ids,J),Interior.electrical_neighbors(I)) for I in sites]
 report=Dict("title"=>cfg["title"],"description"=>cfg["description"],"shape"=>[n,n,n],"extent"=>cfg["extent"],"spacing"=>collect(w.grid.spacing),"time"=>w.time,"params"=>cfg,"metrics"=>Dict("events"=>p.events,"occupied"=>length(sites),"edges"=>edges,"branch_vertices"=>count(>=(3),degrees),"tips"=>count(==(1),degrees),"reached_ground"=>p.reached_ground,"potential_residual"=>p.last_residual,"mass_drift"=>final_mass-initial_mass,"elapsed_seconds"=>elapsed),"solver"=>"geometric correction multigrid; red-black smoothing; fine-grid Dirichlet constraints; independently checked core harmonic residual","growth"=>"same per-face field^eta probability and stateless seed/event hash as DielectricGrowth; each occupied cell consumes .8 nutrient","history"=>history,"source_hash"=>bytes2hex(sha256(read(joinpath(dirname(@__DIR__),"src","electrical.jl")))),"prototype_hash"=>bytes2hex(sha256(read(joinpath(@__DIR__,"branch_solver_prototype.jl")))),"study_hash"=>bytes2hex(sha256(read(@__FILE__))),"optics"=>Dict("channel_radius"=>cfg["optical_radius"],"diameter"=>2cfg["optical_radius"],"geometry"=>"exact face-adjacency native graph; optical tube thickness"))
 open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w")
 args=["--palette","glass-ice","--conductive-graph","--graph-radius",string(cfg["optical_radius"]),"--absorption","0.1","--scattering","0.02","--emission","1.5","--orbit","0.6","--elevation","0.35","--studio-background","slate"]
 open(io->JSON3.pretty(io,Dict("id"=>cfg["slug"],"family"=>"branches","title"=>cfg["title"],"description"=>cfg["description"],"volume"=>"volume","renderer"=>"volume","args"=>args,"size"=>2200,"samples"=>128)),joinpath(out,"render-job.json"),"w")
 println("READY ",out);report
end
if abspath(PROGRAM_FILE)==(@__FILE__) || (isdefined(Main,:INTERIOR_RUNNING_STUDY) && Main.INTERIOR_RUNNING_STUDY=="refinement_branches")
 length(ARGS)==1 || error("usage: study refinement_branches RECIPE.toml")
 fine_branch(ARGS[1];pilot_events=parse(Int,get(ENV,"BRANCH_PILOT_EVENTS","0")))
end
