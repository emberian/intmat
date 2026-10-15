#!/usr/bin/env julia
# Paired whole-native-step experiment on the caller's exact garden recipe.
include(joinpath(@__DIR__,"..","scripts","device_cellular_setup.jl"))
include(joinpath(@__DIR__,"..","scripts","native_constructed_garden.jl"))
include(joinpath(@__DIR__,"..","scripts","device_habitat_fixture.jl"))
isdefined(Interior,:DeviceLifecyclePrefixCarrier) || Base.include(Interior,joinpath(@__DIR__,"device_lifecycle_prefix.jl"))
using Serialization,Statistics
hashfile(path)=bytes2hex(Interior.SHA.sha256(read(path)))
function prefix_gate_compare(a,b;tolerance=2e-5)
    # Same GPU atomic scheduling variability allowed by the ordinary native gate.
    for key in (:active,:ids,:parents,:lineages,:generations,:bonds,:bond_ids,:new_parent)
        a.cells.state[key]==b.cells.state[key] || error("identity/lifecycle mismatch: $key")
    end
    a.cells.environment.metadata[1:6]==b.cells.environment.metadata[1:6] || error("counter mismatch")
    errors=Dict{String,Float64}()
    for (label,x,y) in (("cells",a.cells.state,b.cells.state),("environment",a.cells.environment,b.cells.environment),
            ("habitat",a.arrays,b.arrays),("feedstock",a.feedstock,b.feedstock))
        for key in keys(x)
            if eltype(x[key])<:AbstractFloat
                errors[label*"/"*string(key)]=maximum(abs,Float64.(x[key]).-Float64.(y[key]);init=0.)
            elseif key===:enzyme_bits
                errors[label*"/enzyme"]=maximum(abs,Float64.(reinterpret(Float32,x[key])).-Float64.(reinterpret(Float32,y[key]));init=0.)
            end
        end
    end
    maximum(values(errors);init=0.)<=tolerance || error("numerical comparison exceeds $tolerance")
    errors
end
function prefix_gate_steps!(w,p,count,dt,carrier=nothing;lifecycle=:dense)
    for _ in 1:count
        if carrier===nothing
            Interior.device_cellular_step!(w,p.cells,dt;anabolism=p.habitat,check_status=false,lifecycle)
        else
            Interior.device_lifecycle_prefix_step!(w,carrier,dt;anabolism=p.habitat,check_status=false)
        end
    end
    Interior.device_cellular_status!(p.cells)
    w.time=Float64(p.cells.attempted_steps)*Float64(dt)
end
function prefix_gate_slot_reuse(backend;energy_type=Float32)
    w,p=constructive_device_fixture(backend;energy_type)
    mass=Array(p.cells.state.mass);mass[1]=3f0;mass[2]=.025f0;copyto!(p.cells.state.mass,mass)
    weights=Array(p.cells.state.weights);h=Int(p.cells.geometry.hidden)
    weights[1,:,h+10].=0;weights[1,end,h+10]=5f0;copyto!(p.cells.state.weights,weights)
    original_id=Array(p.cells.state.ids)[2]
    cp=Interior.device_habitat_checkpoint(w,p)
    w2=World(w.grid;backend,seed=w.seed);fill!(w2.fields.nutrient,0)
    p2=Interior.device_habitat_restore(w2,cp);carrier=Interior.DeviceLifecyclePrefixCarrier(p2.cells)
    prefix_gate_steps!(w,p,1,.01f0);prefix_gate_steps!(w2,p2,1,.01f0,carrier)
    a=Interior.device_habitat_snapshot(p);b=Interior.device_habitat_snapshot(p2)
    errors=prefix_gate_compare(a,b)
    a.cells.environment.metadata[3]>=1 && a.cells.environment.metadata[4]>=1 || error("probe did not exercise birth+death")
    a.cells.state.ids[2]!=original_id && a.cells.state.new_parent[2]==1 || error("dead slot was not reused by founder1")
    (;births=a.cells.environment.metadata[3],deaths=a.cells.environment.metadata[4],
        old_slot2_id=original_id,new_slot2_id=a.cells.state.ids[2],highwater=only(Array(carrier.highwater)),errors)
end
function run_prefix_gate(profile,outdir)
    ispath(outdir) && error("fresh output required")
    warm=parse(Int,get(ENV,"PREFIX_GATE_WARM_STEPS","10"))
    count=parse(Int,get(ENV,"PREFIX_GATE_SAMPLE_STEPS","20"))
    samples=parse(Int,get(ENV,"PREFIX_GATE_SAMPLES","3"))
    continuation=parse(Int,get(ENV,"PREFIX_GATE_CONTINUATION_STEPS","8"))
    seconds=parse(Float64,get(ENV,"PREFIX_GATE_WALL_SECONDS","270"))
    1<=warm<=100 && 1<=count<=200 && 1<=samples<=5 && 1<=continuation<=20 && 0<seconds<=300 || error("bounded controls required")
    cfg=TOML.parsefile(profile);dt=Float32(cfg["native"]["dt"])
    mkpath(outdir);start=time();deadline=start+seconds
    provenance=(source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
        experiment_source=Interior.DEVICE_LIFECYCLE_PREFIX_SOURCE,runtime=device_cellular_runtime(),
        profile_sha256=hashfile(profile),runner_sha256=hashfile(@__FILE__),
        source_commit=get(ENV,"DEVICE_SOURCE_COMMIT","unspecified"),loaded_project_sources=Interior.arena_sources(),
        setup_sha256=hashfile(joinpath(@__DIR__,"..","scripts","device_cellular_setup.jl")),
        fixture_sha256=hashfile(joinpath(@__DIR__,"..","scripts","native_constructed_garden.jl")),
        slot_reuse_fixture_sha256=hashfile(joinpath(@__DIR__,"..","scripts","device_habitat_fixture.jl")),
        highwater_reconstruction="production owned workspace; constructor founder count; strict restore maximum nonzero persistent ID/active/new_parent slot; shared dense/prefix kernels maintain device monotonic watermark",
        controls=(warm_steps=warm,sample_steps=count,samples,continuation,wall_seconds=seconds),
        comparison="identity/contact/new_parent and counters exact; all floating roles <=2e-5 (ordinary native gate tolerance); unused free and nondeterministic integer scratch not compared")
    open(io->JSON3.pretty(io,provenance),joinpath(outdir,"provenance.json"),"w")
    function checkdeadline()
        time()<deadline || error("bounded gate deadline reached; no throughput result")
    end
    # Compile/exercise actual death followed by same-step dead-slot reuse separately.
    slot_reuse=prefix_gate_slot_reuse(backend;energy_type=get(cfg["native"],"energy_type","Float32")=="Float64" ? Float64 : Float32);checkdeadline()
    w,p,funding=native_constructed_garden(backend,cfg)
    initial=Interior.device_habitat_checkpoint(w,p)
    w2=World(w.grid,backend,NamedTuple(),NamedTuple(),0.,w.seed)
    p2=Interior.device_habitat_restore(w2,initial);carrier=Interior.DeviceLifecyclePrefixCarrier(p2.cells)
    initial=nothing;GC.gc();checkdeadline()
    prefix_gate_steps!(w,p,warm,dt);prefix_gate_steps!(w2,p2,warm,dt,carrier)
    warm_errors=prefix_gate_compare(Interior.device_habitat_snapshot(p),Interior.device_habitat_snapshot(p2));checkdeadline()
    ordinary=Float64[];prefix=Float64[]
    for sample in 1:samples
        # Reverse paired order on alternate samples to expose drift/order effects.
        for mode in (isodd(sample) ? (:ordinary,:prefix) : (:prefix,:ordinary))
            checkdeadline();t=time_ns()
            if mode===:ordinary
                prefix_gate_steps!(w,p,count,dt);push!(ordinary,(time_ns()-t)/1e9)
            else
                prefix_gate_steps!(w2,p2,count,dt,carrier);push!(prefix,(time_ns()-t)/1e9)
            end
        end
    end
    final_errors=prefix_gate_compare(Interior.device_habitat_snapshot(p),Interior.device_habitat_snapshot(p2));checkdeadline()
    checkpoint=Interior.device_habitat_checkpoint(w2,p2)
    serialize(joinpath(outdir,"checkpoint.jls"),checkpoint)
    highwater=only(Array(carrier.highwater))
    # Restore strict source-pinned base state, then reconstruct study-only geometry.
    wr=World(w.grid,backend,NamedTuple(),NamedTuple(),0.,w.seed)
    restored=Interior.device_habitat_restore(wr,checkpoint)
    expected=(cells=checkpoint.cells.snapshot,arrays=checkpoint.arrays,feedstock=checkpoint.feedstock,
        parameters=checkpoint.parameters,feed_parameters=checkpoint.feed_parameters)
    Interior.device_habitat_snapshot(restored)==expected || error("strict checkpoint restoration changed state")
    restored_carrier=Interior.DeviceLifecyclePrefixCarrier(restored.cells)
    only(Array(restored_carrier.highwater))==highwater || error("checkpoint highwater reconstruction differs")
    checkpoint=nothing;expected=nothing;GC.gc();checkdeadline()
    prefix_gate_steps!(w2,p2,continuation,dt,carrier);prefix_gate_steps!(wr,restored,continuation,dt,restored_carrier)
    continuation_errors=prefix_gate_compare(Interior.device_habitat_snapshot(p2),Interior.device_habitat_snapshot(restored))
    Array(carrier.highwater)==Array(restored_carrier.highwater) || error("continued highwater differs")
    final=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p2));final.finite || error("nonfinite final state")
    material_error=final.material-funding.initial.material-final.supplied+funding.initial.supplied+final.removed-funding.initial.removed
    potential_error=final.potential-funding.initial.potential+final.work-funding.initial.work+final.heat-funding.initial.heat-final.supplied_potential+funding.initial.supplied_potential+final.removed_potential-funding.initial.removed_potential
    rtol=cfg["run"]["ledger_relative_tolerance"]
    abs(material_error)<=max(1e-4,rtol*(funding.initial.material+final.supplied)) || error("material closure failed")
    abs(potential_error)<=max(1e-4,rtol*(funding.initial.potential+final.supplied_potential)) || error("potential closure failed")
    fill!(restored.cells.environment.status,Int32(8))
    before=Interior.device_habitat_snapshot(restored;check_status=false);h=copy(restored_carrier.highwater)
    Interior.device_lifecycle_prefix_step!(wr,restored_carrier,dt;anabolism=restored.habitat,check_status=false)
    after=Interior.device_habitat_snapshot(restored;check_status=false)
    before.cells.state==after.cells.state && before.cells.environment==after.cells.environment && before.arrays==after.arrays && before.feedstock==after.feedstock && Array(h)==Array(restored_carrier.highwater) || error("latched launch changed arrays")
    rejected=false
    try Interior.device_habitat_checkpoint(wr,restored) catch err;err isa ArgumentError || rethrow();rejected=true;end
    rejected || error("partial checkpoint accepted")
    result=(;provenance,slot_reuse,funding,final,material_error,potential_error,warm_errors,final_errors,continuation_errors,
        ordinary_seconds=ordinary,prefix_seconds=prefix,ordinary_ms_per_step=1000median(ordinary)/count,
        prefix_ms_per_step=1000median(prefix)/count,speedup=median(ordinary)/median(prefix),highwater,
        capacity=Int(p2.cells.geometry.capacity),strict_restore_exact=true,highwater_reconstruction_exact=true,
        failure_latch_preserved=true,partial_checkpoint_rejected=true,elapsed_seconds=time()-start)
    open(io->JSON3.pretty(io,result),joinpath(outdir,"gate.json"),"w")
    println("PREFIX_GATE_COMPLETE speedup=",result.speedup," highwater=",highwater," capacity=",result.capacity);flush(stdout)
    result
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("PROFILE FRESH_OUTPUT required")
    run_prefix_gate(ARGS...)
end
