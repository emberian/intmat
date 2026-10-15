# Artifact exporter shared by device prototype; generated from the validated CPU export path.
function export_native_branch(w,p,cfg,out,history,elapsed)
 n=cfg["n"];initial_mass=Float64(n)^3;final_mass=sum(Float64,w.fields.nutrient)+sum(Float64,w.fields.deposit)
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
 report=Dict("title"=>cfg["title"],"description"=>cfg["description"],"shape"=>[n,n,n],"extent"=>cfg["extent"],"spacing"=>collect(w.grid.spacing),"time"=>w.time,"params"=>cfg,"metrics"=>Dict("events"=>p.events,"occupied"=>length(sites),"edges"=>edges,"branch_vertices"=>count(>=(3),degrees),"tips"=>count(==(1),degrees),"reached_ground"=>p.reached_ground,"potential_residual"=>p.last_residual,"mass_drift"=>final_mass-initial_mass,"potential_min"=>minimum(p.potential),"potential_max"=>maximum(p.potential),"elapsed_seconds"=>elapsed),"solver"=>"KA/Metal geometric correction multigrid; red-black smoothing; exact fine-grid Dirichlet constraints; independent core residual checked at checkpoints","growth"=>"same per-face field^eta probability and stateless seed/event hash as DielectricGrowth; each occupied cell consumes .8 nutrient","history"=>history,"source_hash"=>bytes2hex(sha256(read(joinpath(dirname(@__DIR__),"src","electrical.jl")))),"prototype_hash"=>bytes2hex(sha256(read(joinpath(@__DIR__,"branch_solver_prototype.jl")))),"device_solver_hash"=>bytes2hex(sha256(read(joinpath(@__DIR__,"branch_metal_multigrid.jl")))),"study_hash"=>bytes2hex(sha256(read(joinpath(@__DIR__,"refinement_branches_metal.jl")))),"optics"=>Dict("channel_radius"=>cfg["optical_radius"],"diameter"=>2cfg["optical_radius"],"geometry"=>"exact face-adjacency native graph; optical tube thickness"))
 open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w")
 args=["--palette","glass-ice","--conductive-graph","--graph-radius",string(cfg["optical_radius"]),"--absorption","0.1","--scattering","0.02","--emission","1.5","--orbit","0.6","--elevation","0.35","--studio-background","slate"]
 open(io->JSON3.pretty(io,Dict("id"=>cfg["slug"],"family"=>"branches","title"=>cfg["title"],"description"=>cfg["description"],"volume"=>"volume","renderer"=>"volume","args"=>args,"size"=>2200,"samples"=>128)),joinpath(out,"render-job.json"),"w")
 println("READY ",out);report
end
