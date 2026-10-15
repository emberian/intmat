using Interior, JSON3, TOML, SHA, LinearAlgebra, Printf
BLAS.set_num_threads(1)
const ROOT=dirname(@__DIR__)
function bake(recipe)
 cfg=TOML.parsefile(recipe);slug=cfg["slug"];out=joinpath(ROOT,"outputs","collection-20261015",slug)
 isfile(joinpath(out,"run.json")) && return
 mkpath(out);cp(recipe,joinpath(out,"recipe.toml");force=true)
 n=cfg["n"]; w=World(Grid(n,n,n;extent=Tuple(cfg["extent"]));seed=cfg["seed"])
 # Local nutrient availability gates the same field-driven frontier rule.
 for I in CartesianIndices(w.fields.nutrient)
  x,y,z=(Tuple(I).-((n+1)/2,))./(n/2)
  allowed=cfg["material"]=="sheet" ? abs(z-.14sin(5x))<.11 : cfg["material"]=="fork" ? (abs(y)<.12 || abs(y-.55x)<.10 || abs(y+.55x)<.10) : true
  w.fields.nutrient[I]=allowed ? 1f0 : 0f0
 end
 m0=sum(Float64,w.fields.nutrient)+sum(Float64,w.fields.deposit)
 p=Interior.DielectricGrowth(w;seed=cfg["seed"],root=Tuple(cfg["root"]),eta=cfg["eta"],tolerance=2e-5,max_iterations=1600,relaxation=1.86)
 history=Any[];elapsed=@elapsed for i in 1:cfg["events"]
  previous=p.events;step!(w,p,.05f0)
  if p.events%25==0 || p.reached_ground || i==cfg["events"] || previous==p.events
   push!(history,(events=p.events,time=w.time,occupied=count(p.occupied),residual=p.last_residual,iterations=p.last_iterations,reached_ground=p.reached_ground));@printf("%s event=%d residual=%.3g\n",slug,p.events,p.last_residual);flush(stdout)
  end
  (p.reached_ground || previous==p.events) && break
 end
 record_state(w,(p,),joinpath(out,"state");step=p.events)
 sites=findall(p.occupied);ids=Dict(I=>i for (i,I) in enumerate(sites));edges=0
 open(joinpath(out,"channels.obj"),"w") do io
  println(io,"# Exact conductive face adjacency; optical radius is styling")
  for I in sites;println(io,"v ",join(ntuple(d->2*(I[d]-.5)/n-1,3)," "));end
  for (i,I) in enumerate(sites),a in 1:3
   J=I+CartesianIndex(ntuple(d->d==a ? 1 : 0,3))
   if haskey(ids,J);println(io,"l ",i," ",ids[J]);edges+=1;end
  end
 end
 # Preserve all conductors to the grounded box, instead of clipping tips to a sphere.
 f=Interior.electrical_optical_fields(w,p;halo=.03f0,carrier=.028f0)
 fields=merge(f,(mask=ones(Float32,n,n,n),activity=.025f0.*f.activity,potential=copy(p.potential),nutrient=copy(w.fields.nutrient),deposit=copy(w.fields.deposit),frontier=Float32.(p.frontier)))
 export_volume(fields,joinpath(out,"volume"))
 mf=sum(Float64,w.fields.nutrient)+sum(Float64,w.fields.deposit)
 degree=[count(J->haskey(ids,J),Interior.electrical_neighbors(I)) for I in sites]
 meta=Dict("title"=>cfg["title"],"description"=>cfg["description"],"shape"=>[n,n,n],"time"=>w.time,"params"=>cfg,"metrics"=>Dict("events"=>p.events,"occupied"=>length(sites),"edges"=>edges,"branch_vertices"=>count(>=(3),degree),"tips"=>count(==(1),degree),"potential_residual"=>p.last_residual,"initial_mass"=>m0,"final_mass"=>mf,"mass_drift"=>mf-m0,"reached_ground"=>p.reached_ground,"elapsed_seconds"=>elapsed),"history"=>history,"random_state"=>Dict("algorithm"=>"SplitMix64 stateless seed + event counter","seed"=>string(p.seed),"event_counter"=>p.events),"source_hash"=>bytes2hex(sha256(read(joinpath(ROOT,"src","electrical.jl")))),"study_hash"=>bytes2hex(sha256(read(@__FILE__))),"model"=>"Lattice dielectric breakdown; field-weighted face-connected growth; nutrient-gated frontier; conductor one / box zero Dirichlet boundaries","optics"=>Dict("channel_radius"=>cfg["radius"],"channel_diameter"=>2cfg["radius"],"coordinates"=>"normalized cube [-1,1]; anisotropic simulation extent retained in params, shown as normalized lattice","carrier"=>"polished box; dilute amber; exact graph filaments; optical halo only"))
 open(joinpath(out,"run.json"),"w") do io;JSON3.pretty(io,meta);end
 args=["--palette",cfg["palette"],"--conductive-graph","--graph-radius",string(cfg["radius"]),"--polished-carrier","--carrier-shape","box","--absorption","0.4","--scattering","0.07","--emission","1.5","--orbit",string(cfg["orbit"]),"--elevation","0.35","--surface-roughness","0.045"]
 open(joinpath(out,"render-job.json"),"w") do io;JSON3.pretty(io,Dict("id"=>slug,"family"=>"branches","title"=>cfg["title"],"description"=>cfg["description"],"volume"=>"volume","renderer"=>"volume","args"=>args,"size"=>1600,"samples"=>96));end
 println("READY ",out);flush(stdout)
end
for recipe in (isempty(ARGS) ? sort(readdir(joinpath(ROOT,"experiments","collection","branches");join=true)) : ARGS);bake(recipe);end
