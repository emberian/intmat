using Interior
using Interior: JSON3
out=isempty(ARGS) ? "outputs/refinement-20261015/dose-fixture" : ARGS[1]
ispath(out) && error("fixture directory exists")
mkpath(out);w=World(Grid(12,12,12),CPU(),NamedTuple(),NamedTuple(),0.,91)
p=DoseThermalWriter(w;path=:helix,period=.08,stop=.045,radius=.3,height=.2,width=.12,heat=35,diffusion=.0003,cooling=4,threshold=.05,conversion=.37,pulse_rate=17);q=p.writer
dt=.003;frames=Any[];sources=Any[]
for i in 1:25
    t=w.time+Float64(Float32(dt))/2
    envelope=q.pulse_rate==0 ? 1f0 : Float32(.5*(1+cos(2pi*Float64(q.pulse_rate)*t)))
    power=q.start<=t<q.stop ? q.heat*envelope : 0f0
    push!(sources,(;step=i,dt=Float32(dt),position=beam_position(q,t),power))
    step!(w,p,dt)
    q=p.writer
    if i in (1,7,15,25)
        for (name,a) in (("temperature",q.temperature),("dose",q.dose));open(io->write(io,a),joinpath(out,"step-$i-$name.f32"),"w");end
        push!(frames,(;step=i,time=w.time,temperature="step-$i-temperature.f32",dose="step-$i-dose.f32"))
    end
end
n=similar(q.dose);d=similar(q.dose);writer_material!(n,d,p)
open(io->write(io,n),joinpath(out,"final-nutrient.f32"),"w");open(io->write(io,d),joinpath(out,"final-deposit.f32"),"w")
metadata=(;material_readout=(nutrient="final-nutrient.f32",deposit="final-deposit.f32"),schema="interior-dose-thermal-fixture",shape=w.grid.size,extent=w.grid.extent,spacing=w.grid.spacing,order="F",dtype="little-endian Float32",initial="zero temperature/dose",sources,frames,diffusion=q.diffusion,cooling=q.cooling,threshold=q.threshold,width=q.width,conversion=q.conversion,initial_nutrient=p.initial_nutrient,initial_deposit=p.initial_deposit,sources_sha=Interior.arena_sources())
open(io->JSON3.pretty(io,metadata),joinpath(out,"fixture.json"),"w");println(abspath(out))
