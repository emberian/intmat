#!/usr/bin/env julia
# Bounded pilot only: inspect budgets, throughput and headroom before any long run.
include("device_cellular_setup.jl")
include("device_habitat_pilot_fixture.jl")
include("device_habitat_capture.jl")
using Serialization
using Interior.KernelAbstractions: @kernel,@index
steps=parse(Int,get(ENV,"DEVICE_HABITAT_PILOT_STEPS","200"))
nvox=parse(Int,get(ENV,"DEVICE_HABITAT_SIZE","128"))
founders=parse(Int,get(ENV,"DEVICE_HABITAT_FOUNDERS","512"))
capacity=parse(Int,get(ENV,"DEVICE_CAPACITY","8192"))
variant=get(ENV,"DEVICE_HABITAT_VARIANT","closed")
cadence=parse(Int,get(ENV,"DEVICE_HABITAT_CAPTURE","50"))
coarse=parse(Int,get(ENV,"DEVICE_HABITAT_COARSE",string(min(nvox,32))))
budget=parse(Int,get(ENV,"DEVICE_HABITAT_BYTES","2000000000"))
restart_steps=parse(Int,get(ENV,"DEVICE_RESTART_STEPS","20"))
0<steps<=1000 && all(>(0),(nvox,founders,capacity,cadence,coarse,budget,restart_steps)) && nvox%coarse==0 || error("invalid bounded pilot limits")
out=isempty(ARGS) ? "outputs/device-habitat-pilot" : ARGS[1]
ispath(out) && error("fresh output required")
mkpath(dirname(abspath(out)));mkdir(out);mkdir(joinpath(out,"frames"))
function pilot_audit(p)
    x=Interior.device_habitat_snapshot(p);a=Interior.device_habitat_audit(x)
    remaining=sum(Float64.(x.cells.environment.resource_ledger[:,1]).+Float64.(x.cells.environment.resource_ledger_low[:,1]))+
        sum(Float64.(x.feedstock.resource_ledger[:,1]).+Float64.(x.feedstock.resource_ledger_low[:,1]))
    a,remaining
end
function pilot_difference(a,b)
    (identities_equal=a.state.active==b.state.active&&a.state.ids==b.state.ids,
        position=maximum(abs.(a.state.position.-b.state.position)),hidden=maximum(abs.(a.state.hidden.-b.state.hidden)),
        metabolites=maximum(abs.(a.state.metabolites.-b.state.metabolites)),fields=maximum(abs.(a.environment.fields.-b.environment.fields)))
end
function main()
    pilot_started=time_ns()
    construction=@elapsed w,p=device_habitat_pilot_fixture(backend;size=nvox,founders,capacity,variant)
    dt=Float32(.001)
    Float64(dt)<=Interior.stability_limit(w,p) || error("pilot dt exceeds coupled stability limit")
    initial,remaining=pilot_audit(p)
    forecast=floor(Int,(initial.material+remaining)/Float64(p.cells.parameters.minimum_mass))+1
    capacity>=forecast || error("capacity $capacity is below funded-live-mass upper bound $forecast")
    bytes=resident_bytes(p);maxframes=cld(steps,cadence)+2
    forecast_frame_bytes=coarse^3*9*4+capacity*200
    ceil(Int,1.02bytes)+1_000_000+maxframes*forecast_frame_bytes<budget || error("checkpoint+frame forecast exceeds capture budget")
    buffer=Interior.KernelAbstractions.zeros(backend,Float32,coarse^3,9)
    used=pilot_frame!(joinpath(out,"frames","frame-000000.jls"),w,p,buffer,coarse)
    samples=Any[];stepseconds=0.;captureseconds=0.
    allocated=@allocated loopseconds=@elapsed for k in 1:steps
        stepseconds+=@elapsed Interior.device_cellular_step!(w,p.cells,dt;anabolism=p.habitat,check_status=false)
        if k%8==0 || k==steps
            stats=Interior.device_cellular_statistics(p.cells)
            w.time=Float64(dt)*stats.completed_steps
            push!(samples,(step=k,time=w.time,cells=stats))
        end
        if k%cadence==0 || k==steps
            captureseconds+=@elapsed used+=pilot_frame!(joinpath(out,"frames","frame-"*lpad(k,6,'0')*".jls"),w,p,buffer,coarse)
            used+bytes<budget || error("capture budget reached")
        end
    end
    final,_=pilot_audit(p)
    material_error=final.material-initial.material-final.supplied+initial.supplied+final.removed-initial.removed
    potential_error=final.potential-initial.potential-final.supplied_potential+initial.supplied_potential+
        final.removed_potential-initial.removed_potential+final.work-initial.work+final.heat-initial.heat
    final.finite && abs(material_error)<max(1e-5,3e-6initial.material) &&
        abs(potential_error)<max(1e-4,5e-6initial.potential) || error("pilot accounting gate failed")
    checkpoint=Interior.device_habitat_checkpoint(w,p)
    serialize(joinpath(out,"checkpoint.jls"),checkpoint);used+=filesize(joinpath(out,"checkpoint.jls"))
    used+1_000_000<budget || error("serialized capture budget exceeded")
    rw=World(w.grid,backend,NamedTuple(),NamedTuple(),0.,w.seed)
    restored=Interior.device_habitat_restore(rw,checkpoint)
    rs=Interior.device_habitat_snapshot(restored)
    rs.cells==checkpoint.cells.snapshot && rs.arrays==checkpoint.arrays && rs.feedstock==checkpoint.feedstock || error("restore bits differ")
    rs=nothing;checkpoint=nothing;GC.gc()
    for _ in 1:restart_steps
        Interior.device_cellular_step!(w,p.cells,dt;anabolism=p.habitat,check_status=false)
        Interior.device_cellular_step!(rw,restored.cells,dt;anabolism=restored.habitat,check_status=false)
    end
    w.time+=Float64(dt)*restart_steps;rw.time+=Float64(dt)*restart_steps
    a=Interior.device_cellular_snapshot(p.cells);b=Interior.device_cellular_snapshot(restored.cells)
    continuation=pilot_difference(a,b)
    continuation.identities_equal && all(<=(1e-4),(continuation.position,continuation.hidden,continuation.metabolites,continuation.fields)) || error("continuation tolerance exceeded")
    result=(source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
        recipe_source_sha256=(runner=bytes2hex(Interior.SHA.sha256(read(@__FILE__))),
            fixture=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"device_habitat_pilot_fixture.jl")))),
            capture=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"device_habitat_capture.jl"))))),runtime=device_cellular_runtime(),
        recipe=(nvox=nvox,extent=p.cells.geometry.extent,founders=founders,capacity=capacity,variant=variant,
            seed=w.seed,exchange_radius=p.cells.parameters.exchange_radius,density=p.cells.parameters.density,
            minimum_mass=p.cells.parameters.minimum_mass,dt=dt,steps=steps,coarse=coarse,cadence=cadence),initial,final,remaining_source_budget=remaining,
        funded_population_upper_bound=forecast,resident_array_bytes=bytes,stored_bytes=used,capture_budget_bytes=budget,
        source_commit=get(ENV,"DEVICE_SOURCE_COMMIT","unspecified"),restart_steps=restart_steps,construction_seconds=construction,
        loop_seconds=loopseconds,main_full_steps_per_second=steps/loopseconds,whole_pilot_seconds=(time_ns()-pilot_started)/1e9,
        step_enqueue_seconds=stepseconds,capture_seconds=captureseconds,host_allocation_bytes=allocated,
        timing_note="pilot includes first-use compilation; not a warmed throughput benchmark; batch8 observations and capture recorded separately",
        material_error,potential_error,continuation,checkpoint_bits_exact=true,samples)
    open(io->JSON3.pretty(io,result),joinpath(out,"pilot.json"),"w");println(JSON3.write(result))
end
try
    main()
catch err
    open(io->JSON3.pretty(io,(error=sprint(showerror,err),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN)),joinpath(out,"failure.json"),"w")
    rethrow()
end
