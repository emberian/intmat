#!/usr/bin/env julia
include("device_cellular_setup.jl")
include("native_constructed_garden.jl")
Base.include(Interior,joinpath(@__DIR__,"..","studies","device_resource_fused.jl"))
using Statistics,Serialization,KernelAbstractions
function fused_gate_steps!(w,p,carrier,n,dt)
    for _ in 1:n
        if carrier===nothing
            Interior.device_cellular_step!(w,p.cells,dt;anabolism=p.habitat,check_status=false)
        else
            Interior.device_resource_fused_step!(w,carrier,dt;anabolism=p.habitat,check_status=false)
        end
    end
    Interior.device_cellular_status!(p.cells)
    w.time=Float64(p.cells.attempted_steps)*Float64(dt)
end
function fused_gate_compare(a,b)
    for key in (:active,:ids,:parents,:lineages,:generations,:bonds,:bond_ids)
        a.cells.state[key]==b.cells.state[key] || error("identity/contact mismatch: $key")
    end
    a.cells.environment.metadata[1:6]==b.cells.environment.metadata[1:6] || error("counter mismatch")
    errors=Dict{String,Float64}()
    for (name,x,y) in (("cells",a.cells.state,b.cells.state),("environment",a.cells.environment,b.cells.environment),
                      ("habitat",a.arrays,b.arrays),("feedstock",a.feedstock,b.feedstock))
        for key in keys(x)
            if eltype(x[key])<:AbstractFloat
                errors[name*"/"*string(key)]=maximum(abs,Float64.(x[key]).-Float64.(y[key]);init=0.)
            elseif key===:enzyme_bits
                errors[name*"/enzyme"]=maximum(abs,Float64.(reinterpret(Float32,x[key])).-Float64.(reinterpret(Float32,y[key]));init=0.)
            end
        end
    end
    maximum(values(errors))<=2e-5 || error("field/ledger comparison exceeded tolerance")
    errors
end
function fused_gate_endpoint_probe(backend)
    z(T,d...)=KernelAbstractions.zeros(backend,T,d...)
    fields=z(Float32,1,3);fill!(fields,1f0)
    source_masks=z(Float32,1,3);fill!(source_masks,1f0)
    sink_masks=z(Float32,1,3);fill!(sink_masks,1f0)
    quota=z(Float32,3,2);fill!(quota,1f0)
    ledger=z(Float32,3,6);host=zeros(Float32,3,6);host[:,1:2].=2f0^27;copyto!(ledger,host)
    e=(;fields,source_masks,sink_masks,quota,resource_ledger=ledger,
       resource_ledger_low=z(Float32,3,6),status=z(Int32,1))
    q=(source_rates=ntuple(_->2f0^28,3),sink_rates=(0f0,0f0,0f0))
    high=z(Float32,1,6);low=z(Float32,1,6)
    Interior.dc_resource_fused_first!(backend,128)(e,q,.25f0,high,low;ndrange=(1,3))
    Interior.dc_resource_tree_finish!(backend,128)(e,high,low;ndrange=3)
    KernelAbstractions.synchronize(backend)
    h=Array(ledger);l=Array(e.resource_ledger_low)
    only(Array(e.status))==0 || error("endpoint precision probe failed")
    for k in 1:3
        BigFloat(h[k,3])+BigFloat(l[k,3])==BigFloat(2)^26-1 || error("endpoint credit lost unit")
        BigFloat(h[k,1])+BigFloat(l[k,1])==BigFloat(2)^26+1 || error("remaining budget lost unit")
    end
    (before=1.,requested=2.0^26,actual=2.0^26-1,remaining=2.0^26+1,passed=true)
end
function run_resource_fused_gate(profile,outdir)
    ispath(outdir)&&error("fresh output required");mkpath(outdir)
    started=time();cfg=TOML.parsefile(profile);dt=Float32(cfg["native"]["dt"])
    shapes=parse.(Int,split(get(ENV,"FUSED_GATE_SHAPES","64"),','))
    all(n->n in (64,128)||(n==8&&backend isa CPU),shapes) || error("bounded64/128 shapes required (CPU smoke8 allowed)")
    count=parse(Int,get(ENV,"FUSED_GATE_SAMPLE_STEPS","200"))
    1<=count<=200 || error("bounded sample steps required")
    endpoint_probe=fused_gate_endpoint_probe(backend);rows=Any[]
    for n in shapes,fed in (false,true)
        time()-started<220 || break
        c=deepcopy(cfg);c["grid"]["shape"]=[n,n,n]
        #128³ uses a bounded funded fixture, not the long-run recipe's1500+1500 budget.
        n==128 && (c["cellular"]["capacity"]=8192;c["feedstock"]["source_budgets"]=[100.,100.])
        !fed && (c["feedstock"]["source_rates"]=[0.,0.])
        wa,a,funding=native_constructed_garden(backend,c)
        wb,b,_=native_constructed_garden(backend,c)
        carrier=Interior.DeviceResourceFusedCarrier(b.cells)
        GC.gc()
        fused_gate_steps!(wa,a,nothing,40,dt);fused_gate_steps!(wb,b,carrier,40,dt)
        ordinary=Float64[];fused=Float64[]
        # No full snapshots or audits between timing samples.
        for sample in 1:3
            for mode in (isodd(sample) ? (:ordinary,:fused) : (:fused,:ordinary))
                w,p,owner,times=mode===:ordinary ? (wa,a,nothing,ordinary) : (wb,b,carrier,fused)
                t=time_ns();fused_gate_steps!(w,p,owner,count,dt);push!(times,(time_ns()-t)/1e9)
            end
        end
        errors=fused_gate_compare(Interior.device_habitat_snapshot(a),Interior.device_habitat_snapshot(b))
        saved=Interior.device_resource_fused_checkpoint(wb,b,carrier)
        wr=World(wb.grid,backend,NamedTuple(),NamedTuple(),0.,wb.seed)
        restored,rc=Interior.device_resource_fused_restore(wr,saved)
        restart_errors=fused_gate_compare(Interior.device_habitat_snapshot(b),Interior.device_habitat_snapshot(restored))
        fused_gate_steps!(wb,b,carrier,3,dt);fused_gate_steps!(wr,restored,rc,3,dt)
        continuation_errors=fused_gate_compare(Interior.device_habitat_snapshot(b),Interior.device_habitat_snapshot(restored))
        final=Interior.device_habitat_audit(Interior.device_habitat_snapshot(b))
        material_error=final.material-funding.initial.material-final.supplied+funding.initial.supplied+final.removed-funding.initial.removed
        potential_error=final.potential-funding.initial.potential+final.work-funding.initial.work+final.heat-funding.initial.heat-final.supplied_potential+funding.initial.supplied_potential+final.removed_potential-funding.initial.removed_potential
        rtol=c["run"]["ledger_relative_tolerance"]
        final.finite && abs(material_error)<=max(1e-4,rtol*(funding.initial.material+final.supplied)) &&
            abs(potential_error)<=max(1e-4,rtol*(funding.initial.potential+final.supplied_potential)) || error("physical closure failed")
        fill!(restored.cells.environment.status,Int32(8));before=Interior.device_habitat_snapshot(restored;check_status=false)
        Interior.device_resource_fused_step!(wr,rc,dt;anabolism=restored.habitat,check_status=false)
        after=Interior.device_habitat_snapshot(restored;check_status=false)
        before.cells.state==after.cells.state && before.cells.environment==after.cells.environment && before.arrays==after.arrays && before.feedstock==after.feedstock || error("failure latch changed physical state")
        rejected=false
        try Interior.device_resource_fused_checkpoint(wr,restored,rc) catch err;err isa ArgumentError||rethrow();rejected=true;end
        rejected||error("failed checkpoint accepted")
        push!(rows,(shape=n,fed,capacity=Int(b.cells.geometry.capacity),support_strategy="materialized",energy_bits=b.cells.geometry.energy_bits,
            ordinary_seconds=ordinary,fused_seconds=fused,ordinary_ms_per_step=1000median(ordinary)/count,
            fused_ms_per_step=1000median(fused)/count,speedup=median(ordinary)/median(fused),
            errors,restart_errors,continuation_errors,material_error,potential_error,physical=final,
            unused_grants_bytes=sizeof(eltype(carrier.ownership.stages.resource.grants))*length(carrier.ownership.stages.resource.grants),
            retained_memory_saved_bytes=0,strict_source_wrapper=true,failure_latch_preserved=true,partial_checkpoint_rejected=true))
        open(io->JSON3.pretty(io,(source=Interior.DEVICE_RESOURCE_FUSED_SOURCE,endpoint_probe,rows)),joinpath(outdir,"partial.json"),"w")
        a=nothing;b=nothing;carrier=nothing;saved=nothing;restored=nothing;rc=nothing
        p=nothing;owner=nothing;before=nothing;after=nothing;GC.gc()
    end
    result=(source_commit=get(ENV,"DEVICE_SOURCE_COMMIT","unspecified"),source=Interior.DEVICE_RESOURCE_FUSED_SOURCE,
        runner_sha256=bytes2hex(Interior.SHA.sha256(read(@__FILE__))),profile_sha256=bytes2hex(Interior.SHA.sha256(read(profile))),
        controls=(warm_steps=40,sample_steps=count,samples=3,shapes),endpoint_probe,rows,
        completed_cases=length(rows),expected_cases=2*length(shapes),elapsed_seconds=time()-started,
        memory_note="V×6 allocation remains retained; only avoided grant traffic is measured. Compact arena still needsV×3 later.")
    open(io->JSON3.pretty(io,result),joinpath(outdir,"gate.json"),"w")
    println("RESOURCE_FUSED_GATE_COMPLETE cases=",length(rows));flush(stdout)
    result
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("PROFILE FRESH_OUTPUT required")
    run_resource_fused_gate(ARGS...)
end
