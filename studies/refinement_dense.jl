# True spatial refinement: fixed physical extent, no interpolation/upscaling.
using Interior, TOML, Statistics
using Interior: JSON3
"""Writer-specific carrier: only nutrient/deposit; no generic field/scratch promise."""
function compact_writing_world(grid;seed=91,backend=CPU())
    nutrient=Interior.KernelAbstractions.zeros(backend,Float32,grid.size...); fill!(nutrient,1f0)
    deposit=Interior.KernelAbstractions.zeros(backend,Float32,grid.size...)
    World(grid,backend,(;nutrient,deposit),NamedTuple(),0.0,Int(seed))
end
function dense_refinement(config; out=nothing,pilot_steps=0)
    out!==nothing && ispath(out) && error("output already exists: $out")
    c=TOML.parsefile(config); n=Int(c["resolution"]); extent=Tuple(Float32.(c["extent"])); mode=c["mode"]
    compact=mode=="writing" && get(c,"compact_world",false)
    arrays=mode=="writing" ? (compact ? 5 : 19) : 20
    println(JSON3.write((;stage="memory-plan",shape=(n,n,n),extent,cell_count=n^3,resident_arrays=arrays,resident_bytes=arrays*4*n^3,export_peak_additional_bytes=8*n^3)))
    arrays*4*n^3<=parse(Int,get(ENV,"DENSE_MEMORY_BUDGET_BYTES","8000000000")) || error("resident plan exceeds configured memory budget")
    backend=if get(c,"backend","cpu")=="metal"
        Metal.functional() || error("Metal unavailable")
        Metal.MetalBackend()
    else; CPU(); end
    w=(compact ? compact_writing_world : World)(Grid(n,n,n;extent);seed=get(c,"seed",91),backend)
    processes=if mode=="writing"
        (ThermalWriter(w;path=:lissajous,width=get(c,"width",.018),heat=get(c,"heat",30),threshold=get(c,"threshold",.05),conversion=get(c,"conversion",18),period=get(c,"period",4),stop=get(c,"period",4)),)
    elseif mode=="wave"
        (AcousticField(w;speed=.36,frequency=.9,source_amplitude=22,source_radius=.035,source_center=(-.3,0,0),material_feedback=10,activity_gain=.28,intensity_scale=.00015,damping=.055),Deposition(rate=.45,threshold=.18,dissolution=.0002))
    else; error("unknown mode"); end
    dt=Float64(c["dt"]); limit=minimum(stability_limit(w,p) for p in processes)
    dt<=limit || error("configured dt $dt exceeds physical stability limit $limit")
    steps=pilot_steps>0 ? pilot_steps : ceil(Int,c["duration"]/dt)
    step!(w,processes,dt) # warm compiled path, counts as physical step
    timings=Float64[]
    for i in 2:steps
        push!(timings,@elapsed step!(w,processes,dt))
        i%max(1,steps÷10)==0 && (println(JSON3.write((;stage="advance",step=i,time=w.time,seconds_per_step=median(timings)))); flush(stdout))
    end
    report=(;schema="interior-true-spatial-refinement",config=abspath(config),recipe=c,backend=string(typeof(backend)),sources=Interior.arena_sources(),study_sha256=bytes2hex(Interior.SHA.sha256(read(@__FILE__))),field_format="little-endian Float32 x-fastest native grid",mode,shape=w.grid.size,extent=w.grid.extent,spacing=w.grid.spacing,dt,stability_limit=limit,steps,time=w.time,median_step_seconds=isempty(timings) ? nothing : median(timings),evolution_seconds=sum(timings),resident_array_bytes=arrays*4*n^3,deposit_max=maximum(w.fields.deposit),deposit_sum=sum(w.fields.deposit),finite=all(isfinite,w.fields.deposit) && all(p->all(a->all(isfinite,a), p isa ThermalWriter ? (p.temperature,p.dose) : p isa AcousticField ? (p.pressure,p.envelope) : ()),processes))
    if out!==nothing
        ispath(out) && error("output already exists: $out")
        mkpath(out)
        # Stream one field at a time: bounded host staging, full native resolution.
        fields=mode=="writing" ? (("deposit",w.fields.deposit),("temperature",processes[1].temperature),("dose",processes[1].dose)) : (("deposit",w.fields.deposit),("pressure",processes[1].pressure),("intensity",processes[1].envelope))
        for (name,field) in fields
            open(joinpath(out,name*".f32"),"w") do io; write(io,Array(field)); end
        end
        open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w")
    end
    println(JSON3.write(report)); report
end
if abspath(PROGRAM_FILE)==(@__FILE__) || (isdefined(Main,:INTERIOR_RUNNING_STUDY) && Main.INTERIOR_RUNNING_STUDY=="refinement_dense")
    get(TOML.parsefile(ARGS[1]),"backend","cpu")=="metal" && (@eval using Metal)
    Base.invokelatest(dense_refinement,ARGS[1];out=length(ARGS)>=2 && ARGS[2]!="-" ? abspath(ARGS[2]) : nothing,pilot_steps=parse(Int,get(ENV,"DENSE_PILOT_STEPS","0")))
end
