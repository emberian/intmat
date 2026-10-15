#!/usr/bin/env julia
# Intrusive event profile, not a throughput benchmark. Exclusive AMD window only.
include("device_cellular_setup.jl")
include("device_habitat_pilot_fixture.jl")
Base.include(Interior,joinpath(@__DIR__,"..","studies","device_support_fastpaths.jl"))
isdefined(Interior,:DeviceProfileCarrier) || Base.include(Interior,joinpath(@__DIR__,"..","studies","device_launch_profile.jl"))
backend isa AMDGPU.ROCBackend || error("AMD event profiling requires ROCBackend")
length(ARGS)==1 || error("fresh output directory required")
outdir=ARGS[1];ispath(outdir) && error("fresh output directory required");mkpath(outdir)
profile_steps=parse(Int,get(ENV,"DEVICE_PROFILE_STEPS","3"))
support_strategy=:materialized
fastpath=false
settle_steps=parse(Int,get(ENV,"DEVICE_PROFILE_SETTLE","5"))
1<=profile_steps<=10 && 0<=settle_steps<=50 || error("invalid bounded profile step count")
function native_profile_step!(w,p)
    if fastpath
        Interior.device_support_fast_step!(w,p.cells,.001;anabolism=p.habitat,check_status=false)
    else
        Interior.device_cellular_step!(w,p.cells,.001;anabolism=p.habitat,check_status=false)
    end
end
profile_cells(p)=fastpath ? Interior.DeviceSupportFastCarrier(p.cells) : p.cells
function profile_fixture(grid,capacity;founders=512)
    w,p=device_habitat_pilot_fixture(backend;size=grid,founders,capacity,support_strategy)
    for _ in 1:settle_steps
        native_profile_step!(w,p)
    end
    Interior.device_cellular_status!(p.cells)
    w,p
end
# Warm the exact carrier/event observer on a disposable coupled state.
ww,wp=profile_fixture(8,8192;founders=2)
warm_plan=Interior.device_profile_plan(ww,wp.cells,.001;anabolism=wp.habitat)
warm_pool=Interior.device_profile_event_pool(AMDGPU,length(warm_plan))
warm_observer=Interior.DeviceProfileEvents(warm_plan,warm_pool,AMDGPU)
Interior.device_profile_step!(ww,wp.cells,warm_observer,.001;anabolism=wp.habitat)
Interior.device_profile_finish(warm_observer);Interior.device_cellular_status!(wp.cells)
ww=nothing;wp=nothing;warm_pool=nothing;warm_observer=nothing;GC.gc();AMDGPU.synchronize()
function run_profile_phase(grid,capacity,mode)
    w,p=profile_fixture(grid,capacity)
    # All planning, pool allocation and host method warming precede wall timing.
    plan=Interior.device_profile_plan(w,profile_cells(p),.001;anabolism=p.habitat)
    observer=mode===:profile ? Interior.DeviceProfileEvents(plan,
        Interior.device_profile_event_pool(AMDGPU,length(plan)*profile_steps),AMDGPU;steps=profile_steps) : nothing
    begin_event=AMDGPU.HIP.HIPEvent(AMDGPU.stream();do_record=false,timing=true)
    end_event=AMDGPU.HIP.HIPEvent(AMDGPU.stream();do_record=false,timing=true)
    AMDGPU.synchronize()
    wall_start=time_ns();enqueue_start=wall_start
    if mode===:profile
        for _ in 1:profile_steps
            Interior.device_profile_step!(w,profile_cells(p),observer,.001;anabolism=p.habitat)
        end
        enqueue_seconds=(time_ns()-enqueue_start)/1e9
        Interior.device_profile_wait!(observer)
    else
        AMDGPU.HIP.record(begin_event)
        for _ in 1:profile_steps
            native_profile_step!(w,p)
        end
        enqueue_seconds=(time_ns()-enqueue_start)/1e9
        AMDGPU.HIP.record(end_event);AMDGPU.HIP.synchronize(end_event)
    end
    wall_seconds=(time_ns()-wall_start)/1e9
    Interior.device_cellular_status!(p.cells)
    summary=mode===:profile ? Interior.device_profile_finish(observer) : nothing
    stream_seconds=mode===:profile ? summary.stream_seconds : Float64(AMDGPU.HIP.elapsed(begin_event,end_event))
    snapshot=Interior.device_habitat_snapshot(p)
    audit=Interior.device_habitat_audit(snapshot)
    record=(;mode,wall_seconds,enqueue_seconds,stream_seconds,launches_per_step=length(plan),
        active=audit.population,capacity,grid,support_strategy,fastpath,settle_steps,profile_steps,audit,summary)
    record,snapshot
end
function profile_parity(a,b;physical_only=false)
    for key in (:active,:ids,:parents,:lineages,:generations,:bonds,:bond_ids)
        a.cells.state[key]==b.cells.state[key] || error("profile changed identity/contact $key")
    end
    a.cells.environment.metadata==b.cells.environment.metadata || error("profile changed counters")
    errors=Dict{String,Float64}()
    for (label,x,y) in (("cells",a.cells.state,b.cells.state),("environment",a.cells.environment,b.cells.environment),
            ("habitat",a.arrays,b.arrays),("feedstock",a.feedstock,b.feedstock))
        for key in keys(x)
            physical_only && key in (:sites,:support,:support_count,:axis_indices,:axis_distance2,:support_normalization,:support_fallback) && continue
            if eltype(x[key])<:AbstractFloat
                errors[label*"/"*string(key)]=maximum(abs,Float64.(x[key]).-Float64.(y[key]);init=0.)
            elseif key===:enzyme_bits
                errors[label*"/enzyme"]=maximum(abs,Float64.(reinterpret(Float32,x[key])).-Float64.(reinterpret(Float32,y[key]));init=0.)
            end
        end
    end
    maximum(values(errors))==0.0 || error("profile changed numerical trajectory")
    errors
end
# A common physical snapshot across representations excludes support-storage roles.
reference_state=nothing
for (strategy,zero_guard) in ((:materialized,false),(:axis,false),(:materialized,true),(:axis,true))
    global support_strategy=strategy
    global fastpath=zero_guard
    grid=128;capacity=8192
    # Compile both substituted kernels and carrier before timing, on a tiny state.
    tw,tp=profile_fixture(8,8192;founders=2)
    native_profile_step!(tw,tp)
    tplan=Interior.device_profile_plan(tw,profile_cells(tp),.001;anabolism=tp.habitat)
    tobs=Interior.DeviceProfileEvents(tplan,Interior.device_profile_event_pool(AMDGPU,length(tplan)),AMDGPU)
    Interior.device_profile_step!(tw,profile_cells(tp),tobs,.001;anabolism=tp.habitat)
    Interior.device_profile_finish(tobs);Interior.device_cellular_status!(tp.cells)
    tw=nothing;tp=nothing;tobs=nothing;GC.gc();AMDGPU.synchronize()
    native_record,native_state=run_profile_phase(grid,capacity,:native)
    GC.gc();AMDGPU.synchronize()
    profiled_record,profiled_state=run_profile_phase(grid,capacity,:profile)
    errors=profile_parity(native_state,profiled_state)
    cross_errors=reference_state===nothing ? nothing : profile_parity(reference_state,native_state;physical_only=true)
    if reference_state===nothing
        global reference_state=native_state
    end
    result=(runtime=device_cellular_runtime(),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,
        habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
        profiler_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"..","studies","device_launch_profile.jl")))),
        runner_sha256=bytes2hex(Interior.SHA.sha256(read(@__FILE__))),
        fastpath_source=Interior.DEVICE_SUPPORT_FASTPATH_SOURCE,
        native=native_record,profiled=profiled_record,errors,cross_errors,
        event_record_wall_ratio=profiled_record.wall_seconds/native_record.wall_seconds,
        caveats=("Intrusive event recording changes host/stream gaps; compare measured overhead before attributing small stages.",
            "Mixed sense/regulate kernels include physical cues/plasticity as well as support loops.",
            "These short settled512-founder states do not establish dense-contact scaling.",
            "256cubed full detached parity was deliberately omitted under the8GiB host cap."))
    open(io->JSON3.pretty(io,result),joinpath(outdir,"paired-$(strategy)-$(zero_guard).json"),"w")
    println("PAIRED_PROFILE_COMPLETE strategy=",strategy," fastpath=",zero_guard," active=",profiled_record.active)
    native_state=nothing;profiled_state=nothing;GC.gc();AMDGPU.synchronize()
end
