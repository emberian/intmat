#!/usr/bin/env julia
include("device_cellular_setup.jl")
using KernelAbstractions
include(joinpath(@__DIR__,"../studies/native_lifetime_credit.jl"))
length(ARGS)==2||error("PROFILE FRESH_OUTPUT_JSON")
ispath(ARGS[2])&&error("fresh output required")
f=NativeLifetimeCredit.NativeActionCredit.fixture(ARGS[1],(.1,.5))
results=Any[]
for (critic,dt) in Iterators.product((false,true),(.02,.01))
    pairs=[]
    for be in (CPU(),backend)
        w=World(f.grid,be,NamedTuple(),NamedTuple(),0.,f.seed)
        p=Interior.device_habitat_restore(w,f.checkpoint)
        h=NativeLifetimeCredit.sidecar(p;dt,critic,critic_learning=critic)
        origin=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
        start=time();NativeLifetimeCredit.advance!(w,p,h,4.;dt);KernelAbstractions.synchronize(be)
        elapsed=time()-start;snap=Interior.device_habitat_snapshot(p);a=Interior.device_habitat_audit(snap)
        residual=(material=a.material-origin.material,potential=a.potential-origin.potential+a.work-origin.work+a.heat-origin.heat)
        maximum(abs,values(residual))<1e-8||error("coupled paid credit closure failed")
        push!(pairs,(;snapshot=snap,state=map(Array,h.state),elapsed,residual,decisions=h.steps÷h.options.decision_stride,reward_windows=h.steps÷h.options.reward_stride))
    end
    cpu,gpu=pairs
    errors=Dict(string(k)=>maximum(abs,Float64.(cpu.snapshot.cells.state[k]).-Float64.(gpu.snapshot.cells.state[k])) for k in (:energy,:metabolites,:hidden,:actions))
    sidecar=Dict(string(k)=>maximum(abs,Float64.(cpu.state[k]).-Float64.(gpu.state[k])) for k in (:learned,:trace,:drive,:baseline,:previous,:paid,:critic,:prediction))
    maximum(values(errors))<1e-4 && maximum(values(sidecar))<1e-4||error("paid learning CPU/AMD parity failed")
    cpu.decisions==gpu.decisions==20 && cpu.reward_windows==gpu.reward_windows==2||error("physical credit clocks drifted")
    push!(results,(;critic,dt,errors,sidecar,cpu_residual=cpu.residual,gpu_residual=gpu.residual,decisions=cpu.decisions,reward_windows=cpu.reward_windows,cpu_seconds=cpu.elapsed,gpu_seconds=gpu.elapsed))
    println(last(results));flush(stdout)
end
report=(;schema="native-lifetime-credit-gpu-v2",results,runtime=device_cellular_runtime(),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,
    study_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"../studies/native_lifetime_credit.jl")))),
    scope="bounded paid lifetime credit on shared coupled schedule; fixed physical decision cadence at two timesteps; timings include first compilation and are not a throughput claim")
open(io->JSON3.pretty(io,report),ARGS[2],"w")
