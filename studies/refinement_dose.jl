# Three resident arrays, shared ThermalWriter physics, bounded slab exports.
using Interior,TOML,Statistics
using Interior: JSON3
isdefined(Interior,:DoseThermalWriter) || Base.include(Interior,joinpath(@__DIR__,"../src/dose_writer.jl"))
function dose_refinement(config,out;pilot_steps=0)
    c=TOML.parsefile(config); n=Int(c["resolution"]);bytes=12*n^3
    bytes<=get(c,"resident_budget_bytes",22000000000) || error("three-array resident memory plan exceeds recipe budget")
    out!==nothing && ispath(out) && error("output already exists")
    backend=get(c,"backend","cpu")=="metal" ? Metal.MetalBackend() : CPU()
    println(JSON3.write((;stage="memory-plan",shape=(n,n,n),resident_arrays=3,resident_bytes=bytes,host_slab_bytes=3*4*n*n*4)));flush(stdout)
    w=World(Grid(n,n,n;extent=Tuple(c["extent"])),backend,NamedTuple(),NamedTuple(),0.,Int(get(c,"seed",91)))
    kwargs=(; (Symbol(k)=>v for (k,v) in c if k in ("width","heat","diffusion","cooling","threshold","conversion","period","stop","pulse_rate"))...)
    p=Interior.DoseThermalWriter(w;kwargs...);q=p.writer; dt=c["dt"];dt<=stability_limit(w,p) || error("unsafe timestep")
    if out!==nothing
        mkpath(dirname(out));mkdir(out)
        open(io->TOML.print(io,c),joinpath(out,"recipe.toml"),"w")
        cp(@__FILE__,joinpath(out,"study-source.jl"))
        for name in ("writing.jl","dose_writer.jl")
            cp(joinpath(@__DIR__,"../src",name),joinpath(out,"source-"*name))
        end
        open(io->JSON3.pretty(io,(;pid=getpid(),status="running",recipe=c,resident_array_bytes=bytes)),joinpath(out,"running.json"),"w")
    end
    steps=pilot_steps>0 ? pilot_steps : ceil(Int,c["duration"]/dt);times=Float64[]
    for i in 1:steps
        elapsed=@elapsed step!(w,p,dt);i>1 && push!(times,elapsed)
        i%max(1,steps÷10)==0 && (println(JSON3.write((;step=i,time=w.time,median_seconds=isempty(times) ? elapsed : median(times))));flush(stdout))
    end
    q=p.writer # reacquire the current ping-pong handle after advancing
    finite=all(isfinite,q.temperature)&&all(isfinite,q.dose);finite||error("nonfinite writer state")
    deposit_sum=0.;deposit_max=0.
    if out!==nothing
        open(joinpath(out,"dose.f32"),"w") do dose_io
            open(joinpath(out,"temperature.f32"),"w") do temp_io
                open(joinpath(out,"deposit.f32"),"w") do deposit_io
                    for z in 1:4:n
                        range=z:min(n,z+3)
                        dose=Array(view(q.dose,:,:,range));temperature=Array(view(q.temperature,:,:,range))
                        deposit=p.initial_deposit .+ p.initial_nutrient .* (-expm1.(-q.conversion.*dose))
                        deposit_sum+=sum(Float64,deposit);deposit_max=max(deposit_max,maximum(deposit))
                        write(dose_io,dose);write(temp_io,temperature);write(deposit_io,deposit)
                    end
                end
            end
        end
    end
    report=(;schema="interior-dose-owned-refinement",mode="writing",recipe=c,config=abspath(config),shape=w.grid.size,extent=w.grid.extent,spacing=w.grid.spacing,dt,steps,time=w.time,median_step_seconds=isempty(times) ? nothing : median(times),evolution_seconds=sum(times),resident_array_bytes=bytes,finite,deposit_sum,deposit_max,backend=string(typeof(backend)),material_ownership="private fixed uniform inventory; World material unchanged",sources=Interior.arena_sources(),study_sha256=bytes2hex(Interior.SHA.sha256(read(@__FILE__))))
    out!==nothing && open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w")
    out!==nothing && open(io->JSON3.pretty(io,(;pid=getpid(),status="complete",steps,time=w.time)),joinpath(out,"running.json"),"w")
    println(JSON3.write(report));report
end
if abspath(PROGRAM_FILE)==(@__FILE__) || (isdefined(Main,:INTERIOR_RUNNING_STUDY)&&Main.INTERIOR_RUNNING_STUDY=="refinement_dose")
    get(TOML.parsefile(ARGS[1]),"backend","cpu")=="metal" && (@eval using Metal)
    Base.invokelatest(dose_refinement,ARGS[1],length(ARGS)>1&&ARGS[2]!="-" ? abspath(ARGS[2]) : nothing;pilot_steps=parse(Int,get(ENV,"DENSE_PILOT_STEPS","0")))
end
