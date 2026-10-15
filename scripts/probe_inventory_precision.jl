#!/usr/bin/env julia
# Fresh bounded experiments only; never restore historical run checkpoints.
include("device_cellular_setup.jl")
include("native_constructed_garden.jl")
include("device_habitat_fixture.jl")
using Statistics
length(ARGS)==2 || error("PROFILE OUTPUT_JSON required")
results=Dict{String,Any}()
function residual(initial,final)
    (material=final.material-initial.material-final.supplied+initial.supplied+final.removed-initial.removed,
     potential=final.potential-initial.potential+final.work-initial.work+final.heat-initial.heat-final.supplied_potential+initial.supplied_potential+final.removed_potential-initial.removed_potential)
end
function advance_batch(w,p,n,dt)
    for _ in 1:n;Interior.device_cellular_step!(w,p.cells,dt;anabolism=p.habitat,check_status=false);end
    Interior.device_cellular_status!(p.cells)
    w.time=Float64(p.cells.attempted_steps)*Float64(dt)
end
function strict_gate()
    cw,cpu=constructive_device_fixture(CPU();energy_type=Float64,inventory_type=Float64)
    gw,gpu=constructive_device_fixture(backend;energy_type=Float64,inventory_type=Float64)
    initial=Interior.device_habitat_audit(Interior.device_habitat_snapshot(gpu))
    advance_batch(cw,cpu,1000,.01f0);advance_batch(gw,gpu,1000,.01f0)
    cs=Interior.device_habitat_snapshot(cpu);gs=Interior.device_habitat_snapshot(gpu)
    final=Interior.device_habitat_audit(gs);errors=Dict{String,Float64}()
    for (label,x,y) in (("cells",cs.cells.state,gs.cells.state),("environment",cs.cells.environment,gs.cells.environment),("habitat",cs.arrays,gs.arrays),("feedstock",cs.feedstock,gs.feedstock))
        for key in keys(x)
            if eltype(x[key])<:AbstractFloat
                errors[label*"/"*string(key)]=maximum(abs,Float64.(x[key]).-Float64.(y[key]);init=0.)
            elseif label=="cells" && key in (:active,:ids,:parents,:lineages,:generations)
                x[key]==y[key] || error("identity parity failed: $key")
            end
        end
    end
    r=residual(initial,final)
    maximum(values(errors))<1e-3 || error("numeric parity failed")
    abs(r.material)<1e-8 && abs(r.potential)<1e-8 || error("wide fixture accounting failed: $r")
    cp=Interior.device_habitat_checkpoint(gw,gpu)
    rw=World(gw.grid;backend,seed=911);fill!(rw.fields.nutrient,0)
    restored=Interior.device_habitat_restore(rw,cp)
    Interior.device_habitat_snapshot(restored)==gs || error("strict restoration failed")
    advance_batch(gw,gpu,7,.01f0);advance_batch(rw,restored,7,.01f0)
    a=Interior.device_habitat_snapshot(gpu);b=Interior.device_habitat_snapshot(restored)
    restart_error=maximum(maximum(abs,Float64.(a.cells.state[k]).-Float64.(b.cells.state[k]);init=0.) for k in (:mass,:metabolites,:energy))
    restart_error<1e-10 || error("continuation diverged")
    (;initial,final,residual=r,array_errors=errors,exact_restore=true,restart_error)
end
results["gate"]=strict_gate();open(io->JSON3.pretty(io,results),ARGS[2],"w");println("wide_gate_passed");flush(stdout)
for T in (Float32,Float64)
    cfg=TOML.parsefile(ARGS[1]);cfg["native"]["inventory_type"]=string(T)
    w,p,funding=native_constructed_garden(backend,cfg);dt=Float32(cfg["native"]["dt"])
    advance_batch(w,p,100,dt)
    times=Float64[]
    for _ in 1:3
        start=time_ns();advance_batch(w,p,200,dt);push!(times,(time_ns()-start)/1e9)
    end
    final=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
    closure=residual(funding.initial,final)
    final.finite || error("nonfinite native inventory")
    if T===Float64
        abs(closure.material)<1e-8 && abs(closure.potential)<1e-8 || error("native wide inventory closure failed: $closure")
    end
    bytes=sum(sizeof,values(p.cells.state))+sum(sizeof,values(p.cells.environment))+sum(sizeof,values(p.habitat.arrays))+sum(sizeof,values(p.habitat.feedstock.arrays))
    workspace=p.cells.ownership.stages.resource
    scratch_bytes=sizeof(workspace.grants)+sum(sizeof(x.high)+sizeof(x.low) for x in workspace.levels)
    results[string(T)]=(initial=funding.initial,final,residual=closure,sample_steps=200,sample_seconds=times,
        median_ms_per_step=1000median(times)/200,physical_and_snapshot_bytes=bytes,resource_workspace_bytes=scratch_bytes)
    open(io->JSON3.pretty(io,results),ARGS[2],"w");println(T," ",results[string(T)]);flush(stdout)
    p=nothing;w=nothing;GC.gc()
end
results["runtime"]=device_cellular_runtime()
results["source"]=Interior.DEVICE_CELLULAR_SOURCE_PIN
results["habitat_source"]=Interior.DEVICE_HABITAT_SOURCE_PIN
open(io->JSON3.pretty(io,results),ARGS[2],"w")
