using Interior, TOML, Printf
using Interior: JSON3
# Load an optional backend before entering the run function (Julia world-age boundary).
if !isempty(ARGS)
    let idx=findfirst(==("--backend"),ARGS)
        requested=idx===nothing ? get(TOML.parsefile(ARGS[1]),"backend","cpu") : (idx<length(ARGS) ? ARGS[idx+1] : "")
        requested=="metal" && (@eval using Metal)
    end
end
function main()
isempty(ARGS) && error("usage: bigrun.jl PROFILE OUTPUT [--backend cpu|metal] [--steps N] [--dry-run]")
profile=abspath(ARGS[1]); dry="--dry-run" in ARGS
length(ARGS)>=2 && !startswith(ARGS[2],"--") || error("explicit fresh OUTPUT directory required")
out=abspath(ARGS[2]); config=TOML.parsefile(profile); sim=config["simulation"]
function option(name,default)
    i=findfirst(==(name),ARGS); i===nothing ? default : (i<length(ARGS) ? ARGS[i+1] : error("missing value for $name"))
end
backendname=option("--backend",get(config,"backend","cpu"))
backend=if backendname=="metal"; Metal.MetalBackend()
elseif backendname=="cpu"; CPU()
else; error("backend must be cpu or metal") end
kind=get(config,"kind","dense"); kind in ("dense","phase") || error("unsupported kind $kind")
n=Int(sim["resolution"]); n>1 || error("resolution must exceed one")
extent=Float32.(Tuple(sim["extent"])); fraction=Float64(get(sim,"stability_fraction",.8))
0<fraction<1 || error("stability_fraction must be between zero and one")
get(sim,"dt_policy","")=="minimum_process_stability_limit" || error("unsupported timestep policy")
# Dry-run uses a small grid with identical requested spacing, avoiding huge allocations.
allocated_n=dry ? min(n,8) : n
allocated_extent=dry ? ntuple(d->extent[d]*allocated_n/n,3) : extent
w=World(Grid(allocated_n,allocated_n,allocated_n;extent=allocated_extent);seed=Int(get(sim,"seed",17)),backend)
opts(name)=(;(Symbol(k)=>Float32(v) for (k,v) in get(config,name,Dict()) if k!="enabled")...)
phase=nothing
processes=if kind=="phase"
    get(sim,"initialization","analytic common physical modes")=="analytic common physical modes" || error("unsupported phase initialization")
    get(get(config,"phase",Dict()),"mean",0)==0 && get(get(config,"phase",Dict()),"noise",0)==0 || error("analytic refinement initialization requires phase mean=noise=0")
    phase=PhaseSeparation(w;noise=0,opts("phase")...)
    host=Array(phase.phi)
    for k in 1:allocated_n,j in 1:allocated_n,i in 1:allocated_n
        # Defined in physical normalized coordinates, identical across refinements.
        x=(i-.5f0)/allocated_n; y=(j-.5f0)/allocated_n; z=(k-.5f0)/allocated_n
        host[i,j,k]=.03f0*(sin(2f0*Float32(pi)*2x)+sin(2f0*Float32(pi)*3y)+sin(2f0*Float32(pi)*2z))
    end
    copyto!(phase.phi,host)
    deposition=get(config,"deposition",Dict())
    get(deposition,"enabled",false) ? (phase,Deposition(;opts("deposition")...)) : (phase,)
else
    seed_activity!(w)
    (WaveForcing(;opts("wave")...),Excitable(;opts("excitable")...),
     Diffusion(:nutrient,Float32(get(get(config,"transport",Dict()),"rate",.01));material_feedback=Float32(get(get(config,"transport",Dict()),"material_feedback",2))),
     Deposition(;opts("deposition")...))
end
limit=minimum(stability_limit(w,p) for p in processes); dt=Float32(fraction*limit)
isfinite(dt)&&dt>0 || error("no positive stable timestep")
t_end=Float64(sim["t_end"]); isfinite(t_end)&&t_end>0 || error("t_end must be positive and finite")
planned_steps=ceil(Int,t_end/Float64(dt)); cap=parse(Int,option("--steps",string(planned_steps)))
cap>0 || error("step cap must be positive"); steps=min(cap,planned_steps)
capture_interval=Float64(sim["capture_interval"]); observation_interval=Float64(get(sim,"observation_interval",Inf))
isfinite(capture_interval) && capture_interval>0 && observation_interval>0 || error("capture interval must be finite and positive; observation interval must be positive")
requested=String.(sim["capture_fields"])
function selected(name)
    symbol=Symbol(name)
    hasproperty(w.fields,symbol) && return getproperty(w.fields,symbol)
    phase!==nothing && hasproperty(phase,symbol) && return getproperty(phase,symbol)
    error("unknown capture field $name")
end
for name in requested; selected(name); end
buffers=kind=="phase" ? 20 : 16
captures=ceil(Int,steps*Float64(dt)/capture_interval)+1
# record_state stores eight World fields; World scratch is omitted.
observation_fields=kind=="phase" ? 12 : 8
observations=isfinite(observation_interval) ? ceil(Int,steps*Float64(dt)/observation_interval)+1 : 2
budget=(kind=kind,resolution=n,backend=backendname,extent=collect(extent),dt=dt,initial_stability_limit=limit,
    planned_steps=planned_steps,actual_step_cap=steps,scheduled_duration=steps*Float64(dt),resident_field_bytes=buffers*n^3*4,
    capture_bytes=length(requested)*n^3*4,estimated_raw_output_bytes=(captures*length(requested)+observations*observation_fields)*n^3*4,
    note="Timestep fixed after construction and checked every step; state growth can reject it. Raw disk budget excludes metadata, optional optics and renderer. Observations are exports, not restart checkpoints.",
    dry_run=dry,dry_run_spacing_exact=true)
println(JSON3.write(budget)); flush(stdout)
dry && return
ispath(out) && error("output already exists; choose a fresh directory: $out")
mkpath(out); cp(profile,joinpath(out,"profile.toml"))
open(joinpath(out,"budget.json"),"w") do io; JSON3.pretty(io,budget); end
function capture(step)
    all(a->all(isfinite,a),values(w.fields)) || error("nonfinite world at step $step")
    phase===nothing || all(isfinite,phase.phi) || error("nonfinite phase at step $step")
    prefix=@sprintf("frame-%08d",step)
    for name in requested
        open(joinpath(out,"$prefix-$name.f32"),"w") do io; write(io,Array(selected(name))); end
    end
    mass=sum(Float64,Array(w.fields.nutrient))+sum(Float64,Array(w.fields.deposit))
    item=(step=step,scheduled_time=step*Float64(dt),world_time=w.time,mass=mass,deposit=sum(Float64,Array(w.fields.deposit)),prefix=prefix)
    open(joinpath(out,"frames.jsonl"),"a") do io; println(io,JSON3.write(item)); end
    println(JSON3.write(item)); flush(stdout)
end
capture(0); record_state(w,processes,joinpath(out,"initial-state");step=0)
next_capture=capture_interval; next_observe=observation_interval
wall_cadence=Float64(get(sim,"progress_seconds",15)); wall_cadence>0 || error("progress_seconds must be positive")
last_progress=time_ns()
elapsed=@elapsed for step in 1:steps
    step!(w,processes,dt)
    time=step*Float64(dt)
    now=time_ns()
    if (now-last_progress)/1e9>=wall_cadence
        println(JSON3.write((progress_step=step,total_steps=steps,scheduled_time=time,target_time=steps*Float64(dt))))
        flush(stdout); last_progress=now
    end
    if time>=next_capture || step==steps
        capture(step); next_capture=(floor(time/capture_interval)+1)*capture_interval
    end
    if time>=next_observe && step<steps
        record_state(w,processes,joinpath(out,@sprintf("state-%08d",step));step)
        next_observe=(floor(time/observation_interval)+1)*observation_interval
    end
end
record_state(w,processes,joinpath(out,"final-state");step=steps)
render=get(config,"render",Dict())
if get(render,"enabled",false)
    fields=phase===nothing ? optical_fields(w) : phase_optical_fields(w,phase)
    export_volume(fields,joinpath(out,"volume"))
    settings=RenderSettings(; (Symbol(k)=>(k in ("width","height","steps") ? Int(v) : Float64(v)) for (k,v) in render if k!="enabled")...)
    render_volume(fields,joinpath(out,"specimen.ppm");settings)
end
open(joinpath(out,"run.json"),"w") do io; JSON3.pretty(io,(budget=budget,elapsed_seconds=elapsed,final_world_time=w.time,final_scheduled_time=steps*Float64(dt),completed_steps=steps)); end
println("Saved ",out)

end
main()
