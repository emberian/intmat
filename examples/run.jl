using Interior, TOML, Printf
using Interior: JSON3

configpath=length(ARGS)>0 ? abspath(ARGS[1]) : joinpath(@__DIR__,"choir.toml")
out=length(ARGS)>1 ? abspath(ARGS[2]) : joinpath(@__DIR__,"..","outputs","choir")
config=TOML.parsefile(configpath);mkpath(out)
backendname=length(ARGS)>2 ? ARGS[3] : "cpu"
backend=if backendname=="metal"
    @eval using Metal
    Metal.MetalBackend()
elseif backendname=="cpu"
    CPU()
else
    error("backend must be cpu or metal")
end
scene=build_scene(config;backend);w=scene.world;sim=config["simulation"]
dt=Float32(sim["dt"]);steps=Int(sim["steps"]);every=Int(get(sim,"capture_every",50))
every>0 && steps>=0 || error("positive capture_every and nonnegative steps required")
processes=scene.processes
mass()=Float64(sum(Float64,Array(w.fields.nutrient))+sum(Float64,Array(w.fields.deposit)))+(scene.particles===nothing ? 0.0 : sum(Float64,scene.particles.storage))
initial_mass=mass();frames=Any[]
state_dir=mktempdir(out;prefix="states-",cleanup=false)
function capture(index)
    name=@sprintf("frame-%04d",index)
    for field in (:deposit,:activity,:nutrient)
        open(joinpath(out,"$name-$field.f32"),"w") do io
            write(io,Array(getproperty(w.fields,field)))
        end
    end
    state_path=record_state(w,processes,joinpath(state_dir,name);step=index)
    push!(frames,Dict("state"=>relpath(joinpath(state_dir,name,"state.json"),out),"name"=>name,"time"=>w.time,"mass"=>mass(),"solid"=>sum(Float64,Array(w.fields.deposit)),"activity_min"=>minimum(w.fields.activity),"activity_max"=>maximum(w.fields.activity)))
end
capture(0)
elapsed=@elapsed for i in 1:steps
    step!(w,processes,dt)
    if i%every==0 || i==steps
        all(a->all(isfinite,a),values(w.fields)) || error("Nonfinite simulation at step $i")
        capture(i)
        @printf("step %d/%d  t=%.3f  solid=%.3f  mass residual=%.3g\n",i,steps,w.time,sum(w.fields.deposit),mass()-initial_mass)
    end
end
metadata=Dict("title"=>get(config,"title","Interior"),"shape"=>collect(w.grid.size),"extent"=>collect(w.grid.extent),"axis_order"=>"x fastest; cell-centered; periodic domain","frames"=>frames,"config"=>config,"state_directory"=>relpath(state_dir,out),"elapsed_seconds"=>elapsed,"initial_mass"=>initial_mass,"final_mass"=>mass(),"relative_mass_error"=>(mass()-initial_mass)/initial_mass,"julia_version"=>string(VERSION),"backend"=>backendname,"model"=>"Excitable grid plus fixed recurrent particle policies; no training or evolution")
open(joinpath(out,"run.json"),"w") do io;JSON3.pretty(io,metadata);end
cp(configpath,joinpath(out,"config.toml");force=true)
fields=optical_fields(w)
export_volume(fields,joinpath(out,"volume"))
renderopts=(;(Symbol(k)=>(k in ("width","height","steps") ? Int(v) : Float64(v)) for (k,v) in get(config,"render",Dict()))...)
render_volume(fields,joinpath(out,"specimen.ppm");settings=RenderSettings(;renderopts...))
println("Saved simulation, volume and render to ",out)
