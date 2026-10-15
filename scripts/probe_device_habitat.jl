#!/usr/bin/env julia
include("device_cellular_setup.jl")
isdefined(Interior,:DeviceCatalyticHabitat) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_habitat.jl"))
include("device_habitat_fixture.jl")
steps=parse(Int,get(ENV,"DEVICE_HABITAT_STEPS","20"))
0<steps<=100 || error("bounded probe requires1–100steps")
length(ARGS)==1 || error("outputJSONpathrequired")
cw,cpu=constructive_device_fixture(CPU());gw,gpu=constructive_device_fixture(backend)
initial=Interior.device_habitat_audit(Interior.device_habitat_snapshot(gpu))
started=time_ns()
for _ in 1:steps
    advance!(cw,cpu,.01);advance!(gw,gpu,.01)
end
seconds=(time_ns()-started)/1e9
cs=Interior.device_habitat_snapshot(cpu);gs=Interior.device_habitat_snapshot(gpu)
final=Interior.device_habitat_audit(gs)
array_errors=Dict{String,Float64}()
for (group,a,b) in (("cells",cs.cells.state,gs.cells.state),("environment",cs.cells.environment,gs.cells.environment),
        ("habitat",cs.arrays,gs.arrays),("feedstock",cs.feedstock,gs.feedstock))
    for key in keys(a)
        x=getproperty(a,key);y=getproperty(b,key)
        if eltype(x)<:AbstractFloat
            array_errors[group*"/"*string(key)]=isempty(x) ? 0. : maximum(abs,Float64.(x).-Float64.(y))
        elseif key===:enzyme_bits
            array_errors[group*"/enzyme"]=maximum(abs,Float64.(reinterpret(Float32,x)).-Float64.(reinterpret(Float32,y)))
        elseif group=="cells" && key in (:active,:ids,:parents,:lineages,:generations)
            x==y || error("CPU/device identity mismatch $key")
        end
    end
end
material_error=final.material-initial.material-final.supplied+initial.supplied+final.removed-initial.removed
potential_error=final.potential-initial.potential+final.work-initial.work+final.heat-initial.heat-
    final.supplied_potential+initial.supplied_potential+final.removed_potential-initial.removed_potential
maximum(values(array_errors))<1e-4 || error("CPU/device numeric parity exceeded tolerance")
abs(material_error)<1e-4 && abs(potential_error)<2e-4 || error("closed accounting exceeded tolerance")
checkpoint=Interior.device_habitat_checkpoint(gw,gpu)
rw=World(gw.grid;backend,seed=911);fill!(rw.fields.nutrient,0)
restored=Interior.device_habitat_restore(rw,checkpoint)
Interior.device_habitat_snapshot(restored)==gs || error("GPU checkpoint bitwise restore failed")
result=(runtime=device_cellular_runtime(),steps,seconds_including_compilation=seconds,
    timing_note="This is a bounded correctness/compile probe, not a throughput measurement.",
    initial,final,material_error,potential_error,array_errors,checkpoint_bitwise=true,
    source=Interior.DEVICE_HABITAT_SOURCE_PIN,base_source=Interior.DEVICE_CELLULAR_SOURCE_PIN)
open(io->JSON3.pretty(io,result),ARGS[1],"w");println(JSON3.write(result))
