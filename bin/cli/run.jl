using Interior, TOML, Printf
using Interior: JSON3

function parse_run_args(args)
    positional=String[]; step_cap=nothing; render_enabled=true
    i=1
    while i<=length(args)
        arg=args[i]
        if arg=="--steps"
            i+=1; i<=length(args) || error("missing --steps value")
            step_cap=parse(Int,args[i]); step_cap>=0 || error("steps must be nonnegative")
        elseif arg=="--no-render"; render_enabled=false
        elseif startswith(arg,"--"); error("unknown option $arg")
        else; push!(positional,arg)
        end
        i+=1
    end
    positional,step_cap,render_enabled
end
positional,step_cap,render_enabled=parse_run_args(ARGS)
length(positional)<=3 || error("run accepts PROFILE OUTPUT BACKEND")
configpath=length(positional)>0 ? abspath(positional[1]) : joinpath(@__DIR__,"..","..","experiments","materials","choir.toml")
out=length(positional)>1 ? abspath(positional[2]) : joinpath(@__DIR__,"..","..","outputs","choir")
# Refuse old recordings before parsing a recipe, loading a backend or allocating.
# Empty precreated directories are accepted; existing files and symlinks are not.
if ispath(out) || islink(out)
    isdir(out) && !islink(out) && isempty(readdir(out)) ||
        error("output must be a fresh or empty directory: $out")
end
config=TOML.parsefile(configpath)
mkpath(out)
backendname=length(positional)>2 ? positional[3] : "cpu"
backend=if backendname=="metal"
    @eval using Metal
    Metal.MetalBackend()
elseif backendname=="cpu"
    CPU()
else
    error("backend must be cpu or metal")
end
scene=build_scene(config;backend);w=scene.world;sim=config["simulation"]
dt=Float32(sim["dt"]);steps=step_cap===nothing ? Int(sim["steps"]) : step_cap;every=Int(get(sim,"capture_every",50))
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
if render_enabled
fields=optical_fields(w)
export_volume(fields,joinpath(out,"volume"))
renderopts=(;(Symbol(k)=>(k in ("width","height","steps") ? Int(v) : Float64(v)) for (k,v) in get(config,"render",Dict()))...)
render_volume(fields,joinpath(out,"specimen.ppm");settings=RenderSettings(;renderopts...))
end
println("Saved simulation to ",out)
