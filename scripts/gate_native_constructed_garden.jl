#!/usr/bin/env julia
# Bounded gate on the exact production recipe, without changing numerical laws.
include("device_cellular_setup.jl")
include("native_constructed_garden.jl")
using Serialization,Statistics
length(ARGS)==2 || error("PROFILE FRESH_OUTPUT required")
profile,outdir=ARGS;ispath(outdir)&&error("fresh output required");mkpath(outdir)
cfg=TOML.parsefile(profile);dt=Float32(cfg["native"]["dt"])
w,p,funding=native_constructed_garden(backend,cfg)
function steps!(w,p,count,delta)
    for _ in 1:count
        Interior.device_cellular_step!(w,p.cells,delta;anabolism=p.habitat,check_status=false)
    end
    Interior.device_cellular_status!(p.cells)
    w.time=Float64(p.cells.attempted_steps)*Float64(delta)
end
warm_steps=parse(Int,get(ENV,"GARDEN_GATE_WARM_STEPS","100"))
sample_steps=parse(Int,get(ENV,"GARDEN_GATE_SAMPLE_STEPS","200"))
sample_count=parse(Int,get(ENV,"GARDEN_GATE_SAMPLES","5"))
1<=warm_steps<=1000 && 1<=sample_steps<=1000 && 1<=sample_count<=10 || error("invalid bounded gate controls")
steps!(w,p,warm_steps,dt)
samples=Float64[]
for _ in 1:sample_count
    start=time_ns();steps!(w,p,sample_steps,dt);push!(samples,(time_ns()-start)/1e9)
end
cp=Interior.device_habitat_checkpoint(w,p)
serialize(joinpath(outdir,"checkpoint.jls"),cp)
expected=(cells=cp.cells.snapshot,arrays=cp.arrays,feedstock=cp.feedstock,parameters=cp.parameters,feed_parameters=cp.feed_parameters)
final=Interior.device_habitat_audit(expected)
material_error=final.material-funding.initial.material-final.supplied+funding.initial.supplied+final.removed-funding.initial.removed
potential_error=final.potential-funding.initial.potential+final.work-funding.initial.work+final.heat-funding.initial.heat-final.supplied_potential+funding.initial.supplied_potential+final.removed_potential-funding.initial.removed_potential
rtol=cfg["run"]["ledger_relative_tolerance"]
abs(material_error)<=max(1e-4,rtol*(funding.initial.material+final.supplied)) || error("material closure failed")
abs(potential_error)<=max(1e-4,rtol*(funding.initial.potential+final.supplied_potential)) || error("potential closure failed")
final.finite || error("nonfinite final state")
wr=World(w.grid,backend,NamedTuple(),NamedTuple(),0.,w.seed)
restored=Interior.device_habitat_restore(wr,cp)
Interior.device_habitat_snapshot(restored)==expected || error("strict native restore changed arrays")
cp=nothing;expected=nothing;GC.gc()
for _ in 1:8
    steps!(w,p,1,dt);steps!(wr,restored,1,dt)
end
a=Interior.device_habitat_snapshot(p);b=Interior.device_habitat_snapshot(restored)
for key in (:active,:ids,:parents,:lineages,:generations,:bonds,:bond_ids)
    a.cells.state[key]==b.cells.state[key] || error("continuation identity/contact mismatch: $key")
end
a.cells.environment.metadata[1:6]==b.cells.environment.metadata[1:6] || error("continuation counter mismatch")
errors=Dict{String,Float64}()
for (label,x,y) in (("cells",a.cells.state,b.cells.state),("environment",a.cells.environment,b.cells.environment),("habitat",a.arrays,b.arrays),("feedstock",a.feedstock,b.feedstock))
    for key in keys(x)
        if eltype(x[key])<:AbstractFloat
            errors[label*"/"*string(key)]=maximum(abs,Float64.(x[key]).-Float64.(y[key]);init=0.)
        elseif key===:enzyme_bits
            errors[label*"/enzyme"]=maximum(abs,Float64.(reinterpret(Float32,x[key])).-Float64.(reinterpret(Float32,y[key]));init=0.)
        end
    end
end
maximum(values(errors))<=2e-5 || error("native continuation exceeded tolerance")
fill!(restored.cells.environment.status,Int32(8))
failed_before=Interior.device_habitat_snapshot(restored;check_status=false)
Interior.device_cellular_step!(wr,restored.cells,dt;anabolism=restored.habitat,check_status=false)
failed_after=Interior.device_habitat_snapshot(restored;check_status=false)
failed_after.cells.state==failed_before.cells.state && failed_after.cells.environment==failed_before.cells.environment &&
    failed_after.arrays==failed_before.arrays && failed_after.feedstock==failed_before.feedstock || error("latched failure mutated physical state")
rejected=false
try
    Interior.device_habitat_checkpoint(wr,restored)
catch err
    err isa ArgumentError || rethrow()
    global rejected=true
end
rejected || error("partial checkpoint was accepted")
result=(source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,runtime=device_cellular_runtime(),
    profile_sha256=bytes2hex(Interior.SHA.sha256(read(profile))),fixture_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"native_constructed_garden.jl")))),
    runner_sha256=bytes2hex(Interior.SHA.sha256(read(@__FILE__))),funding,final,material_error,potential_error,
    warm_steps,sample_steps,sample_seconds=samples,median_ms_per_step=median(samples)/sample_steps*1000,
    steps_per_second=sample_steps/median(samples),strict_native_restore_exact=true,continuation_steps=8,continuation_errors=errors,
    failure_latch_preserved=true,partial_checkpoint_rejected=true,checkpoint_bytes=filesize(joinpath(outdir,"checkpoint.jls")))
open(io->JSON3.pretty(io,result),joinpath(outdir,"gate.json"),"w")
println("NATIVE_GARDEN_GATE_COMPLETE steps=",p.cells.attempted_steps," population=",final.population," ms=",result.median_ms_per_step);flush(stdout)
