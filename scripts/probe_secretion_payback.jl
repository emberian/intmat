#!/usr/bin/env julia
include("device_cellular_setup.jl")
include(joinpath(@__DIR__,"../studies/native_secretion_payback.jl"))
length(ARGS)==2||error("PROFILE FRESH_OUTPUT_JSON")
ispath(ARGS[2])&&error("fresh output required")
results=Any[]
for raw in (0.,2.),dt in (.02,.01)
    f=NativeSecretionPayback.fixture(ARGS[1];raw)
    pairs=[]
    for be in (CPU(),backend)
        c=NativeSecretionPayback.branch(f;backend=be,investment=0.,dt,horizon=2.)
        p=NativeSecretionPayback.branch(f;backend=be,dt,horizon=2.)
        push!(pairs,(;control=c,pulse=p,comparison=NativeSecretionPayback.compare(c,p)))
    end
    cpu,gpu=pairs
    errors=Dict{String,Float64}()
    for arm in (:control,:pulse),field in (:usable,:uptake_a,:uptake_b,:work,:secreted_a,:secreted_b,:enzyme,:bound,:solubilized)
        c=getproperty(cpu,arm).trace;g=getproperty(gpu,arm).trace
        errors[string(arm)*"_"*string(field)]=maximum(maximum(abs,getproperty(x,field).-getproperty(y,field)) for (x,y) in zip(c,g))
    end
    maximum(values(errors))<1e-4||error("secretion CPU/AMD mismatch")
    cpu.comparison.private_payback==gpu.comparison.private_payback||error("private payback timing mismatch")
    residual=max(maximum(abs(x.potential) for pair in pairs for arm in (:control,:pulse) for x in getproperty(pair,arm).trace),
        maximum(abs(x.material) for pair in pairs for arm in (:control,:pulse) for x in getproperty(pair,arm).trace))
    push!(results,(;raw,dt,errors,residual,private_payback=cpu.comparison.private_payback,
        cpu_final=last(cpu.comparison.differences),gpu_final=last(gpu.comparison.differences)))
    println((;raw,dt,residual,maximum_error=maximum(values(errors)),payback=cpu.comparison.private_payback));flush(stdout)
end
report=(;schema="native-secretion-payback-gpu-v1",results,runtime=device_cellular_runtime(),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,
    habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,study_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"../studies/native_secretion_payback.jl")))),
    scope="matched paid pulse/control, raw substrate present/absent, two physics timesteps; bounded native parity and payback timing")
open(io->JSON3.pretty(io,report),ARGS[2],"w")
