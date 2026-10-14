using Interior, Printf
using Interior: JSON3
# Standalone experiment; fresh worlds receive identical beam paths.
out=if length(ARGS)>0
    abspath(ARGS[1])
else
    root=joinpath(@__DIR__,"..","outputs","writing");mkpath(root)
    mktempdir(root;prefix="run-",cleanup=false)
end
function specimen(living)
    name=living ? "living" : "scaffold"; directory=joinpath(out,name);mkpath(directory)
    w=World(Grid(48,48,48);seed=91)
    writer=ThermalWriter(w;path=:lissajous,width=.038,heat=22,threshold=.08,conversion=16,period=4,stop=4)
    neural=living ? NeuralField(w;seed=91,initial_noise=.004,nutrient_feedback=0,material_feedback=-1.4,diffusion=.001,gain=.65,response_time=.65) : nothing
    processes=living ? (writer,neural,Deposition(rate=.18,threshold=.35,dissolution=0)) : (writer,)
    initialmass=material_mass(w); captures=Any[]
    function capture(label,index)
        record_state(w,processes,joinpath(directory,"state-$label");step=index)
        fields=merge(w.fields,(temperature=writer.temperature,dose=writer.dose))
        neural!==nothing && (fields=merge(fields,(latent1=neural.latent1,latent2=neural.latent2)))
        for (key,value) in pairs(fields)
            open(joinpath(directory,"$label-$key.f32"),"w") do io;write(io,Array(value));end
        end
        push!(captures,(label=label,time=w.time,solid=sum(w.fields.deposit),max_temperature=maximum(writer.temperature),max_dose=maximum(writer.dose)))
    end
    capture("initial",0)
    seconds=@elapsed for i in 1:500
        step!(w,processes,.01f0)
        if i%100==0
            all(a->all(isfinite,a),values(w.fields)) || error("nonfinite physical state")
            all(isfinite,writer.temperature) && all(isfinite,writer.dose) || error("nonfinite writer state")
        end
        i==250 && capture("beam",i)
    end
    capture("final",500)
    f=optical_fields(w); deposit=Array(w.fields.deposit)
    # Amber carrier with written inclusions. Faint scaffold emission is an artistic
    # optical mapping; the dynamics remain encoded in the exported physical fields.
    fields=(density=.6f0.+1.2f0.*deposit,activity=living ? .22f0.*max.(Array(w.fields.activity),0f0).*sqrt.(deposit) : .1f0.*deposit,
            structure=deposit,mask=f.mask,deposit=deposit,nutrient=Array(w.fields.nutrient),temperature=Array(writer.temperature),dose=Array(writer.dose))
    export_volume(fields,joinpath(directory,"volume"))
    render_volume(fields,joinpath(directory,"specimen.ppm");settings=RenderSettings(width=320,height=320,steps=110,emission=5,absorption=.6,scattering=.25,exposure=2.5,orbit=.45,elevation=.25))
    metadata=(shape=collect(w.grid.size),extent=collect(w.grid.extent),dt=.01,steps=500,seed=w.seed,living= living,
              source="Prescribed 3D Lissajous Gaussian heating; phenomenological threshold writing",
              fields="Little-endian Float32, x fastest; initial, beam and final snapshots",
              writer=(path="lissajous",width=writer.width,heat=writer.heat,threshold=writer.threshold,conversion=writer.conversion,period=writer.period,stop=writer.stop),
              captures=captures,elapsed_seconds=seconds,initial_mass=initialmass,final_mass=material_mass(w),relative_mass_error=(material_mass(w)-initialmass)/initialmass)
    open(joinpath(directory,"run.json"),"w") do io;JSON3.pretty(io,metadata);end
    @printf("%s: solid %.3f, relative budget error %.3g; saved %s\n",name,sum(deposit),metadata.relative_mass_error,directory)
end
specimen(false)
specimen(true)
