# Small 2D active-turnover studies. Advance laws live in the composed processes.
using Interior,TOML,SHA,Statistics
using Interior: JSON3,KernelAbstractions
include("living_tissue_turnover.jl")
using .LivingTissueTurnover
include("living_tissue_observations.jl")
using .TissueObservations
const BASE=LivingTissueTurnover.LivingTissueStudy

function run_turnover(recipe,output;backend_kind="cpu")
    backend=BASE.tissue_backend(backend_kind)
    Base.invokelatest(run_turnover_backend,recipe,output,backend,backend_kind)
end
function run_turnover_backend(recipe,output,backend,backend_kind;
        make_base=turnover_base_parameters,make_phase=turnover_parameters,
        make_process=TurnoverTissue,measure=turnover_statistics,
        save_checkpoint=turnover_checkpoint,extra_sources=(),law="living-tissue-turnover-v1")
    c=TOML.parsefile(recipe);dims=Tuple(c["dimensions"]);dt=Float32(c["dt"])
    all(d->0<d<=256,dims) || error("this runner is bounded to sides <=256")
    backend_kind=="cpu" && (dims[3]==1 ? maximum(dims)>128 : maximum(dims)>16) &&
        error("large CPU volume requires a separately coordinated runner; use a reserved GPU backend")
    steps=c["steps"];every=c["capture_every"]
    0<steps<=20000 && every>0 || error("invalid bounded study schedule")
    ispath(output)&&error("fresh output required")
    capture_steps=sort!(collect(union(Set(0:every:steps),Set((steps,)))))
    resident=big(prod(dims))*64
    captures=big(length(capture_steps))*prod(dims)*32+3big(prod(dims))*32+10_000_000
    resident<=get(c,"resident_budget_mb",64)*1024^2 || error("resident memory budget exceeded")
    capture_budget=get(c,"capture_budget_bytes",get(c,"capture_budget_mb",256)*1024^2)
    captures<=capture_budget || error("capture budget exceeded")
    base=make_base(;NamedTuple(Symbol(k)=>v for (k,v) in get(c,"base",Dict()))...)
    phase=make_phase(;NamedTuple(Symbol(k)=>v for (k,v) in get(c,"phase",Dict()))...)
    w=tissue_world(Grid(dims...;extent=Tuple(c["extent"]));seed=c["seed"],backend)
    noise_model=get(c,"noise_model","independent")
    p=make_process(w;base_parameters=base,parameters=phase,reserve=c["reserve"],noise_model)
    Interior.check_dt(w,p,dt)
    initial=measure(w,p);history=Any[];probes=Any[]
    solidity(state)=dims[3]==1 ? state[:,:,1,2].+state[:,:,1,3] : state[:,:,:,2].+state[:,:,:,3]
    observation=TissueTurnoverObservation(solidity(Array(p.tissue.state)))
    mkpath(output);mkdir(joinpath(output,"frames"));cp(recipe,joinpath(output,"recipe.toml"))
    sources=Dict{String,String}()
    for name in ("living_tissue_turnover.jl","living_tissue.jl","living_tissue_runner.jl","living_tissue_observations.jl","living_tissue_turnover_run.jl",extra_sources...)
        source=read(joinpath(@__DIR__,name));write(joinpath(output,name),source);sources["studies/$name"]=bytes2hex(sha256(source))
    end
    for (name,source) in BASE.SHARED_SOURCES;write(joinpath(output,name),source);sources["src/$name"]=BASE.SHARED_SHA[name];end
    function capture(step)
        KernelAbstractions.synchronize(backend)
        state=Array(p.tissue.state);reserve=Array(p.reserve);s=measure(w,p;state,reserve)
        density=solidity(state);events=observe_turnover!(observation,density)
        residual=s.material-initial.material
        s.finite && s.physical_min>=0 && s.neural_maxabs<=1.00001 && abs(residual)<5e-5initial.material || error("turnover state/balance failed")
        (!hasproperty(s,:capacity_excess) || s.capacity_excess<=16*eps(Float32)*s.capacity) || error("solid occupancy capacity exceeded")
        push!(history,merge((step=step,material_residual=residual,
            topology=dims[3]==1 ? tissue_topology(state) : nothing,
            minimum_solid=minimum(density),maximum_solid=maximum(density)),events,s))
        dir=joinpath(output,"frames",lpad(step,7,'0'));BASE.write_tissue_frame(dir,state,w.grid,step,w.time;native=true)
        open(io->write(io,reserve),joinpath(dir,"reserve.f32"),"w")
        volume=JSON3.read(read(joinpath(dir,"volume.json"),String),Dict{String,Any})
        push!(volume["fields"],"reserve")
        open(io->JSON3.pretty(io,volume),joinpath(dir,"volume.json"),"w")
        open(io->println(io,JSON3.write(last(history))),joinpath(output,"metrics.jsonl"),"a")
        (step==0 || step==steps || length(history)%10==0) && (println(JSON3.write(last(history)));flush(stdout))
        push!(probes,(step=step,time=w.time,values=collect(state[cld(dims[1],2),cld(dims[2],2),cld(dims[3],2),:])))
    end
    capture(0)
    # Compilation is separate and does not spend the specimen's finite reserve.
    warm=tissue_world(Grid(3,3,dims[3]==1 ? 1 : 3;extent=(3,3,3));backend);warm_p=make_process(warm;base_parameters=base,parameters=phase)
    step!(warm,warm_p,.001);KernelAbstractions.synchronize(backend)
    elapsed=@elapsed for k in 1:steps
        step!(w,p,dt)
        k%64==0 && KernelAbstractions.synchronize(backend)
        (k%every==0 || k==steps)&&capture(k)
    end
    save_checkpoint(w,p,joinpath(output,"final.jls"));record_state(w,p,joinpath(output,"final-state");step=steps)
    report=(id=c["id"],title=c["title"],description=c["description"],dimensions=dims,extent=w.grid.extent,backend=backend_kind,
        dt,steps,law,parameters=(base=base,phase=phase),resident_array_bytes=Int(resident),estimated_capture_bytes=Int(captures),capture_budget_bytes=capture_budget,elapsed_seconds=elapsed,voxel_steps_per_second=prod(dims)*steps/elapsed,
        source_sha256=sources,frame_steps=capture_steps,history,probes,probe_site=(cld(dims[1],2),cld(dims[2],2),cld(dims[3],2)),
        section=(shape=(min(160,dims[1]),min(160,dims[2]),3),dtype="uint8",array_order="F",channels=("biomass","matrix","excitation"),ranges=((0,1.2),(0,.6),(-1,1))),
        topology_observation="2D only: periodic B+M>.15 components. All dimensions: sampled hysteretic crossings below.10/above.24; a reformation requires a prior observed opening. These are concentration observations, not organisms or mechanical cavities.",
        initialization="Nearly uniform tissue and precursor with 2% $(noise_model) noise; excitation 2% noise, no drawn anatomy or time-dependent forcing. Bandlimited initialization uses exact cell averages of one shared continuous random field.",
        budget="Finite immobile high-potential reserve R→N renews local precursor toward target. N/B/M/R potential1; W potential0. All material conserved; only paid N→W dissipates the declared potential. Abstract signal channels are not chemical species.")
    open(io->JSON3.pretty(io,report),joinpath(output,"run.json"),"w")
    maps=Dict{String,Any}()
    for (name,array) in (("reopenings",Int32.(observation.reopenings)),("reformations",Int32.(observation.reformations)),
        ("minimum_solid",observation.minimum),("maximum_solid",observation.maximum))
        suffix=eltype(array)==Int32 ? ".i32" : ".f32"
        filename=name*suffix;open(io->write(io,array),joinpath(output,filename),"w")
        maps[name]=(file=filename,dtype=eltype(array)==Int32 ? "int32" : "float32")
    end
    open(io->JSON3.pretty(io,(shape=size(observation.occupied),array_order="F",low=.10,high=.24,fields=maps)),joinpath(output,"turnover-maps.json"),"w")
    println("READY ",output)
end
if abspath(PROGRAM_FILE)==@__FILE__;length(ARGS) in (2,3)||error("recipe fresh-output [cpu|amd|metal]");run_turnover(ARGS[1],ARGS[2];backend_kind=length(ARGS)==3 ? ARGS[3] : "cpu");end
