#!/usr/bin/env julia
# Fresh matched states only. Native whole steps, no historical resume or new law.
include("device_cellular_setup.jl")
include("native_constructed_garden.jl")
using Statistics
length(ARGS)==2 || error("PROFILE OUTPUT_JSON required")
function advance_batch(w,p,n,dt)
    for _ in 1:n;Interior.device_cellular_step!(w,p.cells,dt;anabolism=p.habitat,check_status=false);end
    Interior.device_cellular_status!(p.cells);w.time=Float64(p.cells.attempted_steps)*Float64(dt)
end
function residual(a,b)
    (material=b.material-a.material-b.supplied+a.supplied+b.removed-a.removed,
     potential=b.potential-a.potential+b.work-a.work+b.heat-a.heat-b.supplied_potential+a.supplied_potential+b.removed_potential-a.removed_potential)
end
cfg=TOML.parsefile(ARGS[1]);results=Dict{String,Any}()
small=deepcopy(cfg);small["grid"]["shape"]=[16,16,16];small["grid"]["extent"]=[1.,1.,1.]
small["cellular"]["founders"]=4;small["cellular"]["capacity"]=256
for pool in ("resources","feedstock"),key in ("source_rates","source_budgets","sink_rates","sink_budgets")
    small[pool][key]=small[pool][key]./64
end
# CPU/native parity in selected receptor composition.
cw,cp,_=native_constructed_garden(CPU(),small);gw,gp,_=native_constructed_garden(backend,small)
advance_batch(cw,cp,100,.02f0);advance_batch(gw,gp,100,.02f0)
cs=Interior.device_habitat_snapshot(cp);gs=Interior.device_habitat_snapshot(gp)
errors=Dict(string(key)=>maximum(abs,Float64.(cs.cells.state[key]).-Float64.(gs.cells.state[key])) for key in (:inputs,:hidden,:adaptation,:actions,:mass,:metabolites,:energy))
maximum(values(errors))<1e-3 || error("native sensory parity failed: $errors")
check=Interior.device_habitat_checkpoint(gw,gp);rw=World(gw.grid;backend,seed=small["seed"]);fill!(rw.fields.nutrient,0)
rp=Interior.device_habitat_restore(rw,check)
Interior.device_habitat_snapshot(rp)==gs || error("native receptor restart mismatch")
results["gate"]=(errors=errors,exact_restore=true)
cp=nothing;gp=nothing;rp=nothing;GC.gc()
for arm in ("off","totals","species")
    trial=deepcopy(small)
    arm=="off" && delete!(trial,"material_sensing")
    if arm=="species"
        trial["material_sensing"]=Dict("channels"=>["raw_a","raw_b","bound_total"],"scales"=>[2.,2.,.02])
    end
    w,p,funding=native_constructed_garden(backend,trial);trace=Any[]
    for sample in 0:160
        if sample>0;advance_batch(w,p,10,.02f0);end
        snap=Interior.device_habitat_snapshot(p);s=snap.cells.state;active=findall(==(1),s.active);h=Int(p.cells.geometry.hidden)
        push!(trace,(time=w.time,ids=s.ids[active],hidden=s.hidden[active,:],actions=s.actions[active,:],
            material_inputs=s.inputs[active,2h+16:2h+18],audit=Interior.device_habitat_audit(snap)))
    end
    final=last(trace).audit;r=residual(funding.initial,final)
    final.finite && abs(r.material)<1e-8 && abs(r.potential)<1e-8 || error("sensory trajectory accounting failed: $r")
    results[arm]=(config=trial,trace=trace,residual=r)
    open(io->JSON3.pretty(io,results),ARGS[2],"w");println(arm," done ",r);flush(stdout)
    p=nothing;GC.gc()
end
# Independent native64 complete-step performance, same production timestep.
for arm in ("off","totals")
    trial=deepcopy(cfg);arm=="off" && delete!(trial,"material_sensing")
    w,p,funding=native_constructed_garden(backend,trial);dt=Float32(trial["native"]["dt"])
    advance_batch(w,p,100,dt);times=Float64[]
    for _ in 1:3
        start=time_ns();advance_batch(w,p,200,dt);push!(times,(time_ns()-start)/1e9)
    end
    r=residual(funding.initial,Interior.device_habitat_audit(Interior.device_habitat_snapshot(p)))
    abs(r.material)<1e-8 && abs(r.potential)<1e-8 || error("timing closure failure")
    results["timing_"*arm]=(sample_seconds=times,sample_steps=200,median_ms_per_step=1000median(times)/200,residual=r)
    p=nothing;GC.gc()
end
results["runtime"]=device_cellular_runtime();results["source"]=Interior.DEVICE_CELLULAR_SOURCE_PIN
results["habitat_source"]=Interior.DEVICE_HABITAT_SOURCE_PIN
results["driver_sha256"]=bytes2hex(Interior.SHA.sha256(read(@__FILE__)))
open(io->JSON3.pretty(io,results),ARGS[2],"w");println("all material sensing gates complete")
