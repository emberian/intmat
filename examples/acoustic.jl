using Interior
using Interior: JSON3
out=isempty(ARGS) ? joinpath(@__DIR__,"..","outputs","acoustic") : abspath(ARGS[1]); mkpath(out)
w=World(Grid(48,48,48);seed=31)
p=AcousticField(w;speed=.36,frequency=.9,source_amplitude=22,source_radius=.065,
    source_center=(-.3,0,0),material_feedback=10,activity_gain=.28,intensity_scale=.00015,damping=.055)
processes=(p,Excitable(diffusion=.0015,threshold=.18,recovery_rate=.1,gain=2.0),
    Diffusion(:nutrient,.003;material_feedback=8),Deposition(rate=.45,threshold=.18,dissolution=.0002))
mass()=sum(Float64,Array(w.fields.nutrient))+sum(Float64,Array(w.fields.deposit))
dt=.01f0; initial_mass=mass()
frames=Any[]
elapsed=@elapsed for step in 1:600
    step!(w,processes,dt)
    if step%200==0
        all(a->all(isfinite,a),values(w.fields)) && all(isfinite,p.pressure) && all(isfinite,p.envelope) || error("Nonfinite acoustic specimen")
        push!(frames,(step=step,time=w.time,deposit=sum(w.fields.deposit),pressure_max=maximum(abs,w.fields.pressure)))
        println(last(frames))
    end
end
fields=optical_fields(w)
export_volume(fields,joinpath(out,"volume"))
for (name,a) in (("pressure",p.pressure),("intensity",p.envelope),("deposit",w.fields.deposit),("activity",w.fields.activity),("nutrient",w.fields.nutrient))
    open(joinpath(out,"$name.f32"),"w") do io; write(io,Array(a)); end
end
render_volume(fields,joinpath(out,"specimen.ppm");settings=RenderSettings(width=480,height=480,steps=128,absorption=1.4,emission=1.8,exposure=1.6))
open(joinpath(out,"run.json"),"w") do io
    JSON3.pretty(io,(title="Acoustic amber: waves writing mineral",shape=w.grid.size,dt=dt,steps=600,time=w.time,elapsed_seconds=elapsed,
        initial_mass=initial_mass,final_mass=mass(),frames=frames,
        parameters=(speed=p.speed,frequency=p.frequency,source_amplitude=p.source_amplitude,source_radius=p.source_radius,source_center=p.source_center,
            material_feedback=p.material_feedback,damping=p.damping,activity_gain=p.activity_gain,intensity_scale=p.intensity_scale,envelope_time=p.envelope_time),
        model="Damped scalar wave with mineral-dependent face coefficients; intensity stimulates excitable deposition; periodic domain, external transducer power"))
end
