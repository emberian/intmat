# Portability of the actual source-position wrapper and composed existing kernels.
include(joinpath(@__DIR__,"../studies/warm_thread.jl"))
function fixture(backend)
 w=World(Grid(18,18,18);backend,fields=(:nutrient,:deposit,:activity,:pressure),scratch_fields=())
 writers=(KnotWriter(ThermalWriter(w;period=.8,width=.13,heat=25,diffusion=.00012,cooling=8,threshold=.035,conversion=14,stop=.8),1),KnotWriter(ThermalWriter(w;period=.9,width=.13,heat=22,diffusion=.00012,cooling=8,threshold=.035,conversion=12,stop=.8),2))
 acoustic=AcousticField(w;speed=.37,material_feedback=14,source_amplitude=30,frequency=1.7,source_radius=.13,source_center=(-.34,-.1,-.15),activity_gain=0)
 for _ in 1:150;step!(w,(writers...,acoustic),.0032);end
 (;nutrient=Array(w.fields.nutrient),deposit=Array(w.fields.deposit),temperature1=Array(writers[1].writer.temperature),dose1=Array(writers[1].writer.dose),temperature2=Array(writers[2].writer.temperature),dose2=Array(writers[2].writer.dose),pressure=Array(acoustic.pressure))
end
cpu=fixture(CPU());metal=fixture(Metal.MetalBackend());errors=Dict(string(k)=>maximum(abs.(getfield(cpu,k).-getfield(metal,k))) for k in keys(cpu));@assert maximum(values(errors))<2e-5
result=(schema="warm-knot-source-wrapper-cpu-metal",shape=(18,18,18),steps=150,time=150*Float64(.0032f0),maximum_absolute_errors=errors,material_relative_error=(sum(Float64,metal.nutrient)+sum(Float64,metal.deposit)-18^3)/18^3,interpretation="Broader and faster sources exercise both exact knot-position functions, two shared-inventory production thermal kernels and mineral-dependent production acoustic kernel on CPU/Metal. A portability fixture, not artwork-resolution convergence or a benchmark.")
open(io->JSON3.pretty(io,result),ARGS[1],"w");println(JSON3.write(result))
