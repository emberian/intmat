using Interior, Printf, TOML, SHA, Statistics
using Interior: JSON3
# Prescribed toolpaths are imposed heat sources, not emergent trajectories.
const ROOT=normpath(joinpath(@__DIR__,".."))
const RECIPES=joinpath(ROOT,"experiments","collection","writing")
function bake(recipe; n=96, pilot=false)
    cfg=TOML.parsefile(recipe); slug=cfg["id"]
    out=joinpath(ROOT,"outputs","collection-20261015",(pilot ? "pilot-write-" : "write-")*slug); mkpath(out)
    w=World(Grid(n,n,n);seed=cfg["seed"])
    writers=map(cfg["writers"]) do spec
        kwargs=Dict(Symbol(k)=>(k=="path" ? Symbol(v) : k=="center" ? Tuple(v) : v) for (k,v) in spec)
        ThermalWriter(w;kwargs...)
    end
    processes=Tuple(writers)
    if haskey(cfg,"neural")
        neural=NeuralField(w;seed=cfg["seed"],initial_noise=0,nutrient_feedback=0,pressure_feedback=0,
            material_feedback=cfg["neural"]["material_feedback"],gain=cfg["neural"]["gain"],response_time=.55,
            diffusion=cfg["neural"]["diffusion"],recurrence=(1.2,0,0,.5),bias=(-.12,0))
        processes=(processes...,neural,Deposition(rate=Float32(cfg["neural"]["rate"]),threshold=.15f0,dissolution=0f0))
    end
    mass0=material_mass(w); dt=.02f0; steps=round(Int,cfg["duration"]/dt)
    if !pilot; record_state(w,processes,joinpath(out,"state-initial");step=0); end
    elapsed=@elapsed for i in 1:steps; step!(w,processes,dt); end
    all(a->all(isfinite,a),values(w.fields)) || error("Nonfinite world")
    deposit=Array(w.fields.deposit); dose=sum(Array(p.dose) for p in writers)
    metrics=(mass_error=(material_mass(w)-mass0)/mass0,deposit_sum=sum(deposit),deposit_max=maximum(deposit),
       fraction_over_01=count(>(.1),deposit)/length(deposit),fraction_over_025=count(>(.25),deposit)/length(deposit),
       dose_max=maximum(dose),elapsed_seconds=elapsed)
    @show slug n metrics
    if !pilot
        record_state(w,processes,joinpath(out,"state-final");step=steps)
        mask=optical_fields(w).mask
        fields=(density=.3f0.+.45f0.*deposit,activity=.025f0.*deposit,structure=deposit,mask=mask,deposit=deposit,
                nutrient=Array(w.fields.nutrient),dose=dose,temperature=sum(Array(p.temperature) for p in writers))
        export_volume(fields,joinpath(out,"volume"))
        cp(recipe,joinpath(out,"recipe.toml");force=true)
        hashes=Dict(relpath(f,ROOT)=>bytes2hex(sha256(read(f))) for f in [@__FILE__,joinpath(ROOT,"src","writing.jl"),joinpath(ROOT,"src","neural.jl"),recipe])
        metadata=(id="write-"*slug,family="writing",title=cfg["title"],description=cfg["description"],shape=collect(w.grid.size),extent=collect(w.grid.extent),
           time=w.time,dt=dt,steps=steps,parameters=cfg,metrics=metrics,source_sha256=hashes,
           model="Prescribed Gaussian heat source, diffusive cooling and conservative threshold nutrient-to-solid conversion. Optional fixed recurrent neural field; no trained network.",
           optical_mapping="Amber or ice carrier is artistic; inclusion surface is actual deposit=0.18 isosurface.")
        open(joinpath(out,"run.json"),"w") do io; JSON3.pretty(io,metadata); end
        args=["--palette",cfg["palette"],"--polished-carrier","--carrier-shape","ellipsoid","--inclusion-surface","deposit","--inclusion-isovalue","0.18","--emission","0.0","--absorption","0.16","--scattering","0.04","--surface-roughness","0.09","--orbit",string(cfg["orbit"]),"--elevation","0.35"]
        open(joinpath(out,"render-job.json"),"w") do io; JSON3.pretty(io,(id="write-"*slug,family="writing",title=cfg["title"],description=cfg["description"],volume="volume",renderer="volume",args=args,size=1600,samples=96)); end
    end
end
pilot="--pilot" in ARGS
recipes=sort(filter(f->endswith(f,".toml"),readdir(RECIPES;join=true)))
selected=filter(a->!startswith(a,"--"),ARGS)
for recipe in recipes
    isempty(selected)||any(s->occursin(s,basename(recipe)),selected)||continue
    bake(recipe;n=pilot ? 40 : 96,pilot=pilot)
    GC.gc()
end
