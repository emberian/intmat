#!/usr/bin/env julia
include("device_cellular_setup.jl")
include(joinpath(@__DIR__,"../studies/construction_lag_map.jl"))
length(ARGS)==2||error("PROFILE FRESH_OUTPUT_DIRECTORY")
(ispath(ARGS[2]) || ispath(ARGS[2]*"-gate.json"))&&error("fresh output required")
M=ConstructionLagMap
M.register_construction_lag_world_batches!(M.P)
records=M.checkpoints(ARGS[1];seeds=(211,),starts=(7.8,))
cpu=M.cohort(records;backend=CPU(),captures=(0.,.2,.4,.6))
gpu=M.cohort(records;backend,captures=(0.,.2,.4,.6))
errors=Dict{String,Float64}()
for field in (:usable,:energy,:stores,:uptake_a,:uptake_b,:oxidized,:secretion_a,:secretion_b,:requested,:reserve_fraction,:enzyme,:bound,:solubilized,:work,:heat,:supplied)
    errors[string(field)]=maximum(maximum(abs,getproperty(c,field).-getproperty(g,field)) for (a,b) in zip(cpu,gpu) for (c,g) in zip(a.trace,b.trace))
end
maximum(values(errors))<1e-4 || error("causal pulse CPU/AMD mismatch")
all(all(t.resident_valid) for x in (cpu,gpu) for a in x for t in a.trace)||error("pulse parity resident changed")
residual=maximum(max(abs(t.material),abs(t.potential)) for x in (cpu,gpu) for a in x for t in a.trace)
println((gate="pulse CPU/AMD",errors,residual));flush(stdout)
gate=(;schema="native-construction-lag-gpu-v1",errors,residual,cpu,gpu,runtime=device_cellular_runtime(),
    source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
    scope="seven matched pulse/control worlds across source exchange at8, actual learned checkpoint7.8, pulse.2 and common future innovations; six full14-world map cohorts run only after this gate")
open(io->JSON3.pretty(io,gate),ARGS[2]*"-gate.json","w")
M.run(ARGS[1],ARGS[2];backend)
