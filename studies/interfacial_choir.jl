#= Nacre fault: recurrent mineral + pressure, sharing one changing material.
The perforated lamellae are deliberately seeded geometry. Subsequent images are saved
Julia states. No external waveform or surface animation is applied at render time.
=#
using Interior, Metal, SHA, TOML, Statistics
using Interior: JSON3, RecurrentMaterial
using Interior.KernelAbstractions: @kernel, @index

@kernel function pressure_to_memory!(state,pressure,dt,gain)
    i,j,k=@index(Global,NTuple)
    u=state[i,j,k,1];v=state[i,j,k,2]
    drive=tanh(3f0*pressure[i,j,k])
    state[i,j,k,1]=u+dt*gain*(drive-u)
    state[i,j,k,2]=v+dt*.2f0*gain*(-drive-v)
end
function initial(n;aperture=.85f0,frequency=1f0)
    deposit=zeros(Float32,n,n,n);food=similar(deposit); latent=zeros(Float32,n,n,n,8)
    host=zeros(Float32,n,n,n)
    for k=1:n,j=1:n,i=1:n
        x,y,z=ntuple(d->2f0*((i,j,k)[d]-.5f0)/n-1,3)
        theta=atan(z,x);r=sqrt(x*x+z*z)
        # Three folded perforated lamellae share actual reaction/transport state.
        # Unequal folds and displaced apertures expose their different depths.
        envelope=clamp((1f0-(x/.84f0)^8-(z/.86f0)^8)*8f0,0f0,1f0)
        support=0f0
        for layer=-1:1
            middle=.25f0*layer+.15f0*sin(4x+3z+layer*.65f0)+.08f0*cos(9z-2x)
            holes=clamp((aperture-cos(frequency*(9x+3z)+layer*.8f0)-cos(frequency*(8z-2x)-layer*.65f0))*3f0,0f0,1f0)
            # Broad oblique tears are deliberate initial boundaries.
            tear=clamp((.91f0-abs(x+.16f0*layer+.17f0*sin(4z)))*12f0,0f0,1f0)
            sheet=exp(-((y-middle)/.033f0)^2)*holes*envelope*tear
            support=max(support,sheet)
        end
        seed=.48f0+.14f0*sin(4theta+5y)+.08f0*cos(10r-3z)
        deposit[i,j,k]=support*seed
        food[i,j,k]=support*(1f0-seed)
        host[i,j,k]=support
        for pair=1:4
            phase=2theta+(7f0+pair)*r+3y+pair*.8f0+.3f0*sin(5z)
            latent[i,j,k,2pair-1]=.75f0*sin(phase)
            latent[i,j,k,2pair]=.75f0*cos(phase)
        end
    end
    deposit,food,latent,host
end
function run(config,out; initializer=initial,
             source_files=("src/recurrent_material.jl","src/acoustics.jl","studies/interfacial_choir.jl"),
             initialization_text="Three analytically seeded, folded perforated lamellae with unequal depth and displaced apertures; coherent latent phases. Geometry is seeded, not claimed to emerge from noise.")
    c=TOML.parsefile(config); n=c["resolution"];dt=Float32(c["dt"]);steps=c["steps"];capture=c["capture_every"]
    ispath(out) && error("Refusing existing run: $out")
    mkpath(out);cp(config,joinpath(out,"recipe.toml"));cp(@__FILE__,joinpath(out,"study-source.jl"))
    # Pin the actual loaded implementations beside this new run.
    hashes=Dict{String,String}()
    for f in source_files
        data=read(joinpath(dirname(@__DIR__),f));hashes[f]=bytes2hex(sha256(data));write(joinpath(out,basename(f)),data)
    end
    w=World(Grid(n,n,n;extent=(64,64,64));backend=MetalBackend(),seed=813)
    p=RecurrentMaterial(w;coupling=.23,feedback=1.4,writing=1.1,erosion=.48,contrast=.8)
    sound=AcousticField(w;speed=3.4,material_feedback=14,damping=.035,source_amplitude=2.8,frequency=.19,source_radius=2.5,source_center=(-17.,-6.,-2.),activity_gain=0)
    d,f,l,host=initializer(n;aperture=Float32(get(c,"aperture",.85)),frequency=Float32(get(c,"perforation_frequency",1.)));copyto!(w.fields.deposit,d);copyto!(w.fields.nutrient,f);copyto!(p.state,l)
    Interior.seed_acoustic!(w,sound;amplitude=1.7,radius=3.5,center=(-13.,-4.,-1.),dt)
    d=nothing;f=nothing;l=nothing;GC.gc()
    open(io->write(io,host),joinpath(out,"initial-support.f32"),"w")
    initialmass=sum(Float64,Array(w.fields.deposit))+sum(Float64,Array(w.fields.nutrient))
    history=Any[]; manifest=Any[]
    function snapshot(s)
        dir=joinpath(out,"frames",lpad(s,5,'0'));mkpath(dir)
        solid=Array(w.fields.deposit);activity=Array(w.fields.activity);pressure=Array(sound.pressure)
        fields=(solid=solid,activity=activity,pressure=pressure)
        for (name,a) in pairs(fields);open(io->write(io,a),joinpath(dir,string(name)*".f32"),"w");end
        open(io->JSON3.pretty(io,(shape=w.grid.size,extent=w.grid.extent,fields=collect(keys(fields)),dtype="float32",order="F",time=w.time,step=s)),joinpath(dir,"volume.json"),"w")
        mass=sum(Float64,solid)+sum(Float64,Array(w.fields.nutrient))
        metric=(step=s,time=w.time,material_error=mass-initialmass,solid_mean=mean(solid),solid_max=maximum(solid),nutrient_min=minimum(w.fields.nutrient),pressure_extrema=extrema(pressure),finite=all(isfinite,solid)&&all(isfinite,pressure)&&all(isfinite,Array(p.state)))
        metric.finite || error("nonfinite")
        push!(history,metric);push!(manifest,(step=s,time=w.time,directory=relpath(dir,out)))
        println(JSON3.write(metric));flush(stdout)
    end
    snapshot(0)
    elapsed=@elapsed for s=1:steps
        # Acoustic propagation reads the current mineral. Recurrence then writes
        # the shared mineral, and bounded pressure forcing changes its next state.
        advance!(w,sound,dt);advance!(w,p,dt)
        pressure_to_memory!(w.backend,256)(p.state,sound.pressure,dt,Float32(c["pressure_coupling"]);ndrange=w.grid.size)
        w.time+=Float64(dt)
        (s%capture==0 || s==steps) && snapshot(s)
    end
    Metal.synchronize()
    open(io->write(io,Array(p.state)),joinpath(out,"final-latent.f32"),"w")
    open(io->write(io,Array(w.fields.nutrient)),joinpath(out,"final-nutrient.f32"),"w")
    for (name,array) in (("final-previous-pressure",sound.previous),("final-pressure-envelope",sound.envelope))
        open(io->write(io,Array(array)),joinpath(out,name*".f32"),"w")
    end
    report=(initial_material_mass=initialmass,id=c["id"],title=c["title"],description=c["description"],family="interfaces",shape=w.grid.size,extent=w.grid.extent,time=w.time,steps=steps,dt=dt,elapsed_seconds=elapsed,source_sha256=hashes,params=c,history=history,frames=manifest,simulation="Eight fixed recurrent channels, conservative nutrient/solid exchange, scalar pressure with mineral-dependent face speed; bounded pressure-to-latent forcing. Driven invented material, not calibrated biology or acoustics.",initialization=initialization_text,optics="Isosurfaces of saved solid and pressure; chosen opal/copper/glass materials. No solved optical feedback, physical luminescence or material stress birefringence.")
    open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w")
    open(io->JSON3.pretty(io,(shape=w.grid.size,fields=["solid","activity","pressure"])),joinpath(out,"volume.json"),"w")
    println("READY ",out)
end
function verify(output)
    n=24;dt=.025f0;d,f,l,_=initial(n)
    function fixture(backend;pressure_coupling=.75f0,feedback=14f0)
        w=World(Grid(n,n,n;extent=(64,64,64));backend,seed=813)
        p=RecurrentMaterial(w;coupling=.23,feedback=1.4,writing=1.1,erosion=.48,contrast=.8)
        sound=AcousticField(w;speed=3.4,material_feedback=feedback,damping=.035,source_amplitude=2.8,frequency=.19,source_radius=2.5,source_center=(-17.,-6.,-2.),activity_gain=0)
        copyto!(w.fields.deposit,d);copyto!(w.fields.nutrient,f);copyto!(p.state,l)
        Interior.seed_acoustic!(w,sound;amplitude=1.7,radius=3.5,center=(-13.,-4.,-1.),dt)
        before=sum(Float64,d)+sum(Float64,f)
        for _=1:160
            advance!(w,sound,dt);advance!(w,p,dt)
            pressure_to_memory!(w.backend,256)(p.state,sound.pressure,dt,pressure_coupling;ndrange=w.grid.size)
            w.time+=Float64(dt)
        end
        (solid=Array(w.fields.deposit),nutrient=Array(w.fields.nutrient),latent=Array(p.state),pressure=Array(sound.pressure),initialmass=before)
    end
    cpu=fixture(Interior.CPU());gpu=fixture(MetalBackend());uncoupled=fixture(Interior.CPU();pressure_coupling=0f0);unretuned=fixture(Interior.CPU();feedback=0f0)
    errors=Dict(string(key)=>maximum(abs.(getfield(cpu,key).-getfield(gpu,key))) for key in (:solid,:nutrient,:latent,:pressure))
    maximum(values(errors))<5e-5 || error("CPU/Metal parity failed: $errors")
    balance=sum(Float64,gpu.solid)+sum(Float64,gpu.nutrient)-gpu.initialmass
    abs(balance)/gpu.initialmass<1e-6 || error("material balance")
    minimum(gpu.nutrient)>=0 || error("negative nutrient")
    effect=maximum(abs.(cpu.solid.-uncoupled.solid));retuning=maximum(abs.(cpu.pressure.-unretuned.pressure))
    effect>1e-4 && retuning>1e-3 || error("couplings have no measured effect")
    result=(shape=(n,n,n),steps=160,dt=dt,cpu_metal_maximum_errors=errors,material_relative_error=balance/gpu.initialmass,minimum_nutrient=minimum(gpu.nutrient),pressure_to_mineral_max_effect=effect,mineral_to_pressure_max_effect=retuning,interpretation="Matched ablations establish both feedback directions in this fixture; they do not establish emergent intelligence, physical acoustic calibration, or an ecological conclusion.")
    open(io->JSON3.pretty(io,result),output,"w");println(JSON3.write(result))
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    ARGS[1]=="--verify" ? verify(ARGS[2]) : run(ARGS[1],ARGS[2])
end
