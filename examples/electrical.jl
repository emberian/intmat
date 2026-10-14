using Interior, JSON3, Printf
# Standalone inclusion permits use before the shared module registration lands.
isdefined(Interior,:DielectricGrowth) || Base.include(Interior,joinpath(@__DIR__,"..","src","electrical.jl"))
out=length(ARGS)>0 ? abspath(ARGS[1]) : joinpath(@__DIR__,"..","outputs","electrical")
n=length(ARGS)>1 ? parse(Int,ARGS[2]) : 40
requested_events=length(ARGS)>2 ? parse(Int,ARGS[3]) : 280
mkpath(out);w=World(Grid(n,n,n);seed=29)
p=Interior.DielectricGrowth(w;seed=29,eta=3.0,event_interval=.05,tolerance=2e-5)
initial_mass=Float64(sum(Float64,w.fields.nutrient)+sum(Float64,w.fields.deposit))
frames=Any[]
function capture()
 name=@sprintf("frame-%04d",p.events)
 for field in (:deposit,:activity,:nutrient)
  open(joinpath(out,"$name-$field.f32"),"w") do io;write(io,Array(getproperty(w.fields,field)));end
 end
 open(joinpath(out,"$name-potential.f32"),"w") do io;write(io,p.potential);end
 push!(frames,Dict("name"=>name,"events"=>p.events,"time"=>w.time,"occupied"=>count(p.occupied),"potential_residual"=>p.last_residual,"reached_ground"=>p.reached_ground))
end
capture()
elapsed=@elapsed for i in 1:requested_events
 step!(w,p,.05f0)
 if p.events%25==0 || i==requested_events || p.reached_ground
  capture();@printf("event %d / %d occupied=%d residual=%.3g reached_ground=%s\n",p.events,requested_events,count(p.occupied),p.last_residual,p.reached_ground)
 end
 p.reached_ground && break
end
fields=Interior.electrical_optical_fields(w,p)
export_volume(fields,joinpath(out,"volume"))
meta=Dict("title"=>"A field writes a tree / dielectric breakdown surrogate","shape"=>collect(w.grid.size),"extent"=>collect(w.grid.extent),"frames"=>frames,"events"=>p.events,"requested_events"=>requested_events,"seed"=>29,"eta"=>p.eta,"event_interval"=>p.event_interval,"deposit_per_cell"=>p.deposit_per_cell,"solver_tolerance"=>p.tolerance,"solver_max_iterations"=>p.max_iterations,"potential_residual"=>p.last_residual,"reached_ground"=>p.reached_ground,"initial_mass"=>initial_mass,"final_mass"=>sum(Float64,w.fields.nutrient)+sum(Float64,w.fields.deposit),"elapsed_seconds"=>elapsed,"random_state"=>Dict("algorithm"=>"SplitMix64 hash of seed + event * golden ratio constant","seed_u64"=>string(p.seed),"event_counter"=>p.events),"model"=>"Quasi-static lattice dielectric-breakdown; driven conductor at1, grounded box at0; fixed eta; no learned rule or calibrated voltage","optics"=>"Spherical viewing crop, dilute amber carrier, direct deposited-channel emission with one-cell six-neighbor halo")
open(joinpath(out,"run.json"),"w") do io;JSON3.pretty(io,meta);end
render_volume(fields,joinpath(out,"specimen.ppm");settings=RenderSettings(width=360,height=360,steps=100,emission=6.,scattering=.12,absorption=.65,exposure=1.5))
println("Saved electrical tree: ",out)
