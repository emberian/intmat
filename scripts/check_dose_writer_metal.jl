# Correctness only: 12³, 25 steps; do not use concurrent timings as a benchmark.
using Interior,Metal
using Interior: JSON3
fixture=isempty(ARGS) ? "outputs/refinement-20261015/dose-fixture" : ARGS[1]
m=JSON3.read(read(joinpath(fixture,"fixture.json"),String),Dict{String,Any})
w=World(Grid(12,12,12),Metal.MetalBackend(),NamedTuple(),NamedTuple(),0.,91)
p=DoseThermalWriter(w;path=:helix,period=.08,stop=.045,radius=.3,height=.2,width=.12,heat=35,diffusion=.0003,cooling=4,threshold=.05,conversion=.37,pulse_rate=17)
errors=Any[]
for i in 1:25
    step!(w,p,.003)
    for frame in m["frames"]
        frame["step"]==i || continue
        for (name,array) in (("temperature",p.writer.temperature),("dose",p.writer.dose))
            expected=reshape(copy(reinterpret(Float32,read(joinpath(fixture,frame[name])))),12,12,12)
            difference=maximum(abs.(Array(array).-expected));@assert difference<=2f-5
            push!(errors,(;step=i,field=name,max_abs_error=difference))
        end
    end
end
n=similar(p.writer.dose);d=similar(n);writer_material!(n,d,p)
for (name,array) in (("nutrient",n),("deposit",d))
    expected=reshape(copy(reinterpret(Float32,read(joinpath(fixture,m["material_readout"][name])))),12,12,12)
    difference=maximum(abs.(Array(array).-expected));@assert difference<=2f-5
    push!(errors,(;step=25,field=name,max_abs_error=difference))
end
result=(;schema="dose-writer-metal-parity",shape=w.grid.size,steps=25,errors,timing_claim=false)
println(JSON3.write(result));open(io->JSON3.pretty(io,result),joinpath(fixture,"metal-parity.json"),"w")
