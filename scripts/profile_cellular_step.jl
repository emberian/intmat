#!/usr/bin/env julia
# Read-only inputs: strict checkpoints from the exact pinned Interior source.
# Run with --project=PINNED_SOURCE. All clones happen outside measured regions.
# Replay and host resource controls: docs/cellular-performance-results.json.
# Use fresh output directories; never point this helper at live soak directories.
using Interior, JSON3, Profile, InteractiveUtils, SHA, Serialization, TOML

function cellular(a)
    only(p for p in a.processes if p isa CellularEcology)
end
function state_summary(a)
    p=cellular(a)
    (;step=a.steps,time=a.world.time,population=length(p.ids),bonds=length(p.bonds.endpoints),
      grid=a.world.grid.size,processes=[string(nameof(typeof(p))) for p in a.processes],
      mechanical_substeps=p.mechanical===nothing ? 0 : p.mechanical.substeps,
      parameter_field_type=string(fieldtype(typeof(p),:parameters)),
      parameter_value_type=string(typeof(p.parameters)))
end
function measure(f)
    t=@timed f()
    (;seconds=t.time,bytes=t.bytes,gc_seconds=t.gctime)
end
function measure_steps!(a,n)
    t=@timed for _ in 1:n;arena_step!(a);end
    (;seconds=t.time,bytes=t.bytes,gc_seconds=t.gctime)
end
function component_times(a)
    result=Any[]
    push!(result,(stage="validate_composition",result=measure(()->Interior.validate_composition(a.processes))))
    for p in a.processes
        push!(result,(stage="check_dt/"*string(nameof(typeof(p))),result=measure(()->Interior.check_dt(a.world,p,a.dt))))
    end
    for p in a.processes
        push!(result,(stage="advance/"*string(nameof(typeof(p))),result=measure(()->Interior.advance!(a.world,p,Float32(a.dt)))))
    end
    result
end
function sustained_run!(a,seconds)
    before=state_summary(a);steps=a.steps;GC.gc()
    # POSIX CLOCKS_PER_SEC is 1e6; this lane is run on Linux hbox.
    cpu_before=ccall(:clock,Clong,())/1e6
    started=time_ns();deadline=started+round(UInt64,seconds*1e9)
    t=@timed while time_ns()<deadline;arena_step!(a);end
    cpu_seconds=ccall(:clock,Clong,())/1e6-cpu_before
    (;before,after=state_summary(a),steps=a.steps-steps,wall_seconds=t.time,cpu_seconds,
      steps_per_second=(a.steps-steps)/t.time,bytes=t.bytes,gc_seconds=t.gctime,
      artificial_sleep=false)
end
function profile_case(name,path,out;advance=0,allow_source_change=false,quick=false,sustained=false)
    base=restore_arena(path;allow_source_change)
    for _ in 1:advance;arena_step!(base);end
    before=state_summary(base)
    println(JSON3.write((;case=name,state=before,checkpoint=path)));flush(stdout)
    warm=deepcopy(base)
    for _ in 1:3;arena_step!(warm);end
    measure_steps!(deepcopy(base),8) # Compile the measurement wrapper outside samples.
    measurements=Any[]
    for _ in 1:5
        a=deepcopy(base);GC.gc()
        push!(measurements,measure_steps!(a,8))
    end
    # Reproduce the full process order, recording each existing call separately.
    component_times(deepcopy(base)) # Warm each instrumentation wrapper separately.
    components=component_times(deepcopy(base))
    if !quick
        a=deepcopy(base);GC.gc();Profile.clear()
        @profile for _ in 1:100;arena_step!(a);end
        open(joinpath(out,name*"-profile-flat.txt"),"w") do io
            Profile.print(io;format=:flat,sortedby=:count,mincount=4,C=false)
        end
        open(joinpath(out,name*"-profile-tree.txt"),"w") do io
            Profile.print(io;format=:tree,mincount=8,C=false)
        end
        open(joinpath(out,name*"-warntype.txt"),"w") do io
            p=cellular(base)
            if isdefined(Interior,:_ecology_step!)
                code_warntype(io,Interior._ecology_step!,Tuple{typeof(base.world),typeof(p),Float32,typeof(p.parameters)};optimize=false)
            else
                code_warntype(io,Interior.ecology_step!,Tuple{typeof(base.world),typeof(p),Float32};optimize=false)
            end
        end
    end
    # Separate, identical same-state continuation for baseline/candidate parity.
    a=deepcopy(base)
    for _ in 1:100;arena_step!(a);end
    statefile=joinpath(out,name*"-after100.jls")
    open(io->serialize(io,a),statefile,"w")
    throughput=sustained ? sustained_run!(deepcopy(base),3.0) : nothing
    result=(;case=name,checkpoint=path,checkpoint_sha256=bytes2hex(sha256(read(path))),
        derived_steps=advance,state=before,steps_per_sample=8,measurements,components,
        profile_steps=quick ? 0 : 100,continuation_steps=100,continuation_end=state_summary(a),
        continuation_sha256=bytes2hex(sha256(read(statefile))),allow_source_change,
        loaded_sources=Interior.arena_sources(),throughput)
    open(io->JSON3.write(io,result),joinpath(out,name*".json"),"w")
    println(JSON3.write(result));flush(stdout)
end

1<=length(ARGS)<=5 || error("usage: profile_cellular_step.jl PROFILE_DIRECTORY [--quick] [--allow-source-change] [--sustain] [--large]")
out=abspath(first(ARGS));mkpath(out)
flags=ARGS[2:end]
all(x->x in ("--quick","--allow-source-change","--sustain","--large"),flags) || error("unknown profiling flag")
Profile.init(n=5_000_000,delay=0.0005)
cases=if "--large" in flags
    path=joinpath(out,"large1024.jls")
    if !isfile(path)
        source=dirname(dirname(pathof(Interior)))
        config=TOML.parsefile(joinpath(source,"experiments/arena/cellular-soak-broad.toml"))
        config["grid"]["shape"]=[256,192,1]
        config["grid"]["extent"]=[128.,96.,1.]
        config["cellular"]["population"]=1024
        # One broad stripe remains one stripe over the larger physical domain.
        config["resources"]["niche_waves"]=[1,0,0]
        checkpoint!(build_arena(config),path)
        open(io->JSON3.write(io,config),joinpath(out,"large1024-config.json"),"w")
    end
    (("large1024","large1024.jls",0),)
else
    (("mature-sparse","mature-sparse.jls",0),
                          ("mature-bonded","mature-bonded.jls",0),
                          ("early-clustered","early-initial.jls",475))
end
for (name,file,advance) in cases
    profile_case(name,joinpath(out,file),out;advance,
        quick="--quick" in flags,allow_source_change="--allow-source-change" in flags,
        sustained="--sustain" in flags)
end
