#!/usr/bin/env julia
# Source-pinned selective-energy gate; no long production launch.
include(joinpath(@__DIR__,"..","studies","gate_device_lifecycle_prefix.jl"))
function energy_precision_residual(final,initial)
    (material=final.material-initial.material-final.supplied+initial.supplied+final.removed-initial.removed,
     potential=final.potential-initial.potential+final.work-initial.work+final.heat-initial.heat-final.supplied_potential+initial.supplied_potential+final.removed_potential-initial.removed_potential)
end
function run_energy_precision_gate(profile,outdir)
    ispath(outdir) && error("fresh output required");mkpath(outdir)
    started=time();cfg=TOML.parsefile(profile);dt=Float32(cfg["native"]["dt"])
    cfg["native"]["energy_type"]=="Float64" || error("explicit Float64 profile required")
    # Full-state equality, restored continuation, failed-latch and slot reuse.
    prefix=run_prefix_gate(profile,joinpath(outdir,"prefix"))
    GC.gc()
    config32=deepcopy(cfg);config32["native"]["energy_type"]="Float32"
    w32,p32,_=native_constructed_garden(backend,config32)
    w64,p64,funding=native_constructed_garden(backend,cfg)
    for (w,p) in ((w32,p32),(w64,p64));prefix_gate_steps!(w,p,20,dt);end
    count=100;samples32=Float64[];samples64=Float64[]
    for sample in 1:3
        for mode in (isodd(sample) ? (32,64) : (64,32))
            w,p,rows=mode==32 ? (w32,p32,samples32) : (w64,p64,samples64)
            t=time_ns();prefix_gate_steps!(w,p,count,dt);push!(rows,(time_ns()-t)/1e9)
        end
    end
    # Continue the actual F64 trajectory with the exact-law prefix optimization.
    p32=nothing;w32=nothing;GC.gc()
    carrier=Interior.DeviceLifecyclePrefixCarrier(p64.cells)
    target=parse(Int,get(ENV,"ENERGY_GATE_TREND_STEPS","10000"))
    1<=target<=12000 || error("bounded trend required")
    history=NamedTuple[]
    function observe!()
        final=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p64))
        final.finite || error("nonfinite native energy gate")
        residual=energy_precision_residual(final,funding.initial)
        push!(history,(step=p64.cells.attempted_steps,time=w64.time,residual,physical=final))
        open(io->JSON3.pretty(io,history),joinpath(outdir,"trend.json"),"w")
        rtol=Float64(cfg["run"]["ledger_relative_tolerance"])
        abs(residual.material)<=max(1e-4,rtol*(funding.initial.material+final.supplied)) || error("material closure failed")
        abs(residual.potential)<=max(1e-4,rtol*(funding.initial.potential+final.supplied_potential)) || error("potential closure failed")
    end
    observe!();trend_start=p64.cells.attempted_steps;trend_seconds=Float64[]
    while p64.cells.attempted_steps-trend_start<target && time()-started<250
        n=min(1000,target-(p64.cells.attempted_steps-trend_start))
        t=time_ns();prefix_gate_steps!(w64,p64,n,dt,carrier);push!(trend_seconds,(time_ns()-t)/1e9)
        observe!()
    end
    cp=Interior.device_habitat_checkpoint(w64,p64)
    serialize(joinpath(outdir,"checkpoint.jls"),cp)
    wr=World(w64.grid,backend,NamedTuple(),NamedTuple(),0.,w64.seed)
    restored=Interior.device_habitat_restore(wr,cp)
    Interior.device_habitat_snapshot(restored)==Interior.device_habitat_snapshot(p64) || error("trend strict restart differs")
    result=(source_commit=get(ENV,"DEVICE_SOURCE_COMMIT","unspecified"),runner_sha256=bytes2hex(Interior.SHA.sha256(read(@__FILE__))),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
        profile_sha256=bytes2hex(Interior.SHA.sha256(read(profile))),prefix,
        timing_comparison="matched initial recipe with different battery precision; no Float32/Float64 trajectory-equivalence claim",timing_warm_steps=20,timing_sample_steps=count,
        float32_seconds=samples32,float64_seconds=samples64,
        float32_ms_per_step=1000median(samples32)/count,float64_ms_per_step=1000median(samples64)/count,
        float64_over_float32=median(samples64)/median(samples32),
        trend_steps=p64.cells.attempted_steps-trend_start,trend_target=target,trend_seconds,
        trend_completed=p64.cells.attempted_steps-trend_start==target,
        history,strict_trend_restore_exact=true,energy_bytes=sizeof(eltype(p64.cells.state.energy))*length(p64.cells.state.energy),
        elapsed_seconds=time()-started)
    open(io->JSON3.pretty(io,result),joinpath(outdir,"gate.json"),"w")
    println("ENERGY_GATE_COMPLETE ratio=",result.float64_over_float32," trend_steps=",result.trend_steps," residual=",last(history).residual);flush(stdout)
    result
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("PROFILE FRESH_OUTPUT required")
    run_energy_precision_gate(ARGS...)
end
