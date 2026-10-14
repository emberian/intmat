using Interior
using JSON3
using Statistics
out=isempty(ARGS) ? get(ENV,"INTERIOR_PHASE_OUTPUT",joinpath(@__DIR__,"..","outputs","phase")) : abspath(ARGS[1])
ispath(joinpath(out,"initial")) && error("Phase observations already exist; pass a fresh output directory or set INTERIOR_PHASE_OUTPUT")
mkpath(out)
w=World(Grid(32,32,32;extent=(32,32,32));seed=83)
p=PhaseSeparation(w;mean=-0.12,noise=0.24,kappa=0.7,
    nutrient_diffusion=0.04,interface_activity=0.8)
processes=(p,Deposition(rate=0.08,threshold=0.04,dissolution=0))
dt=0.003f0
metrics=[]
function observe(label,step)
    dir=joinpath(out,label)
    record_state(w,processes,dir;step)
    push!(metrics,Dict("label"=>label,"step"=>step,"time"=>w.time,
        "phase_mass"=>sum(Float64,p.phi),"energy"=>phase_energy(w,p),
        "composition_std"=>std(p.phi),"interface_mean"=>mean(p.interface),
        "deposit_mass"=>sum(Float64,w.fields.deposit)))
end
observe("initial",0)
for i in 1:18000
    step!(w,processes,dt)
    i==9000 && observe("demixed",i)
end
observe("coarsened",18000)
fields=phase_optical_fields(w,p)
export_volume(fields,out)
image=render_volume(fields;settings=RenderSettings(width=320,height=320,steps=120,cutaway=0.65))
# write_ppm is the renderer's dependency-free preview writer.
Interior.write_ppm(joinpath(out,"phase.ppm"),image)
open(joinpath(out,"metrics.json"),"w") do io; JSON3.pretty(io,metrics); end
println(JSON3.write(metrics))
