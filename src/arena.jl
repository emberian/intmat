"""Externally funded nutrient exchange. Rates are inventory per cell per simulation-time;
outflow removes a fraction `1-exp(-sink_rate*dt)` of dissolved nutrient only."""
mutable struct Chemostat
    source_rate::Float64
    sink_rate::Float64
    supplied::Float64
    removed::Float64
end
Chemostat(;source_rate=0.01,sink_rate=0.01)=Chemostat(source_rate,sink_rate,0,0)
function validate(p::Chemostat)
    all(x->isfinite(x)&&x>=0,(p.source_rate,p.sink_rate,p.supplied,p.removed)) || throw(ArgumentError("invalid chemostat rates or budgets"))
end
function advance!(w::World,p::Chemostat,dt)
    validate(p)
    isfinite(Float32(p.source_rate*dt)) || throw(ArgumentError("chemostat source overflows Float32"))
    all(x->isfinite(x)&&x>=0 && isfinite(x+Float32(p.source_rate*dt)),w.fields.nutrient) || throw(ArgumentError("invalid nutrient or chemostat overflow"))
    for i in eachindex(w.fields.nutrient)
        old=w.fields.nutrient[i]
        incoming=Float32(p.source_rate*dt)
        after=old+incoming
        new=Float32(Float64(after)*exp(-p.sink_rate*dt))
        # Ledger follows actual Float32 inventory changes, including rounding.
        p.supplied+=Float64(after)-Float64(old)
        p.removed+=Float64(after)-Float64(new)
        w.fields.nutrient[i]=new
    end
    p
end
record_process_fields(p::Chemostat)=(source_rate=p.source_rate,sink_rate=p.sink_rate,supplied=p.supplied,removed=p.removed)
mutable struct Arena
    world::World
    processes::Tuple
    dt::Float64
    steps::Int
    paused::Bool
    stopped::Bool
    initial_mass::Float64
    last_command::String
    config::Dict{String,Any}
end
# Follow references once: wrappers borrow pools/media, they do not own copies.
function arena_inventory(w,processes)
    total=sum(Float64,w.fields.nutrient)+sum(Float64,w.fields.deposit)
    seen=IdDict{Any,Bool}()
    function visit(p)
        haskey(seen,p) && return 0.0
        seen[p]=true
        mass=p isa MetabolicEcology ? 0.0 : owned_material_mass(p)
        for name in (:controller,:resources,:boundary)
            hasproperty(p,name) || continue
            child=getproperty(p,name)
            child===nothing || (mass+=visit(child))
        end
        mass
    end
    total+sum(visit(p) for p in processes;init=0.0)
end
function arena_energy(processes)
    seen=IdDict{Any,Bool}()
    function visit(p)
        haskey(seen,p) && return 0.0
        seen[p]=true
        if p isa ResourcePool
            return resource_potential(p)
        elseif p isa MetabolicEcology
            c=p.controller
            return sum(Float64,c.energy;init=0.0)+c.heat+p.chemical_heat+visit(p.resources)
        elseif p isa CellularEcology
            return sum(p.energy)+p.work_spent+p.chemical_heat+
                sum(p.metabolites[s,i]*p.resources.potentials[s] for s in 1:2,i in eachindex(p.ids);init=0.)+visit(p.resources)
        end
        0.0
    end
    sum(visit(p) for p in processes;init=0.0)
end
function Arena(w::World,processes::Tuple;dt=0.05,config=Dict{String,Any}())
    validate_composition(processes)
    foreach(p->check_dt(w,p,dt),processes)
    cfg=Dict{String,Any}(config)
    cfg["_initial_energy"]=arena_energy(processes)
    Arena(w,processes,Float64(dt),0,false,false,arena_inventory(w,processes),"",cfg)
end
function arena_step!(a::Arena)
    a.paused || a.stopped || (step!(a.world,a.processes,a.dt); a.steps+=1)
    a
end
# Capture the implementation that this module loaded, rather than mutable disk
# contents at checkpoint time. Julia invalidates included source dependencies on
# the next module load, while a running arena keeps its original fingerprint.
const ARENA_LOADED_SOURCES = Dict(basename(p)=>bytes2hex(SHA.sha256(read(p)))
    for p in sort(readdir(@__DIR__;join=true)) if endswith(p,".jl"))
arena_sources() = copy(ARENA_LOADED_SOURCES)
function arena_atomic(f,path)
    mkpath(dirname(abspath(path)))
    temp,io=mktemp(dirname(abspath(path)))
    try
        f(io); flush(io); close(io); mv(temp,path;force=true)
    finally
        isopen(io) && close(io)
        isfile(temp) && rm(temp)
    end
    path
end
"""Atomic trusted-local restart checkpoint. Serializes complete state and aliases;
restore requires the exact Julia version and implementation hashes by default."""
function checkpoint!(a::Arena,path::AbstractString)
    envelope=(schema="interior-arena-checkpoint",version=1,julia=string(VERSION),sources=arena_sources(),arena=a)
    arena_atomic(io->Serialization.serialize(io,envelope),path)
end
function restore_arena(path::AbstractString;allow_source_change=false)
    e=open(Serialization.deserialize,path)
    e.schema=="interior-arena-checkpoint" && e.version==1 || throw(ArgumentError("unsupported checkpoint schema"))
    e.julia==string(VERSION) || throw(ArgumentError("checkpoint Julia version differs"))
    allow_source_change || e.sources==arena_sources() || throw(ArgumentError("checkpoint implementation differs; explicitly opt into migration"))
    e.arena isa Arena || throw(ArgumentError("checkpoint contains no Arena"))
    e.arena
end
function arena_snapshot(a::Arena;max_side=96,max_agents=3000)
    max_side>=1 && max_agents>=0 || throw(ArgumentError("invalid snapshot bounds"))
    w=a.world; dims=w.grid.size; k=cld(dims[3],2)
    xs=collect(1:max(1,cld(dims[1],max_side)):dims[1]); ys=collect(1:max(1,cld(dims[2],max_side)):dims[2])
    fields=Dict{String,Any}()
    addfield(name,arr)= fields[string(name)]=Dict{String,Any}("shape"=>[length(xs),length(ys)],"values"=>vec(Array(arr[xs,ys,k:k])))
    for name in propertynames(w.fields); addfield(name,getproperty(w.fields,name)); end
    agents=Any[]; metrics=Dict{String,Any}("material"=>arena_inventory(w,a.processes),"nutrient_total"=>sum(Float64,w.fields.nutrient))
    supplied=0.0; removed=0.0
    for p in a.processes
        if p isa CellularEcology
            cellular_arena_snapshot!(fields,agents,metrics,p,w,xs,ys,k,max_agents)
            continue
        end
        if p isa Chemostat; supplied+=p.supplied; removed+=p.removed; end
        if p isa ResourcePool
            for (name,arr) in zip(p.names,p.fields)
                addfield(name,arr); fields[string(name)]["units"]="material / cell"
            end
            stats=resource_statistics(p)
            metrics["resources"]=Dict(string(k)=>v for (k,v) in pairs(stats))
            supplied+=sum(p.supplied); removed+=sum(p.removed)
        end
        if p isa BoundaryMaintenance
            addfield(:film,p.film); fields["film"]["units"]="material / cell"
            metrics["boundary_material"]=boundary_material_mass(p)
            for name in (:contributed_material,:spent_work,:recycled_material)
                metrics["boundary_"*string(name)]=getproperty(p,name)
            end
            p=p.medium
        end
        metabolic=p isa MetabolicEcology ? p : nothing
        if metabolic!==nothing
            for (key,value) in pairs(ecology_statistics(metabolic)); metrics[string(key)]=value; end
            p=metabolic.controller
        end
        for name in (:permeability,:signal,:temperature,:phi,:interface)
            hasproperty(p,name) && getproperty(p,name) isa AbstractArray{<:Real,3} && addfield(name,getproperty(p,name))
        end
        if hasproperty(p,:ids) && hasproperty(p,:positions)
            for i in 1:min(length(p.ids),max_agents-length(agents))
                cell=p.positions[:,i]
                ag=Dict{String,Any}("id"=>p.ids[i],"cell_position"=>cell,"position"=>[cell[d]*w.grid.spacing[d]-w.grid.extent[d]/2 for d in 1:3])
                for (key,name) in (("lineage",:lineages),("generation",:generations),("energy",:energy),("storage",:storage))
                    hasproperty(p,name) && (ag[key]=getproperty(p,name)[i])
                end
                for (key,name) in (("hidden",:hidden),("actions",:last_actions))
                    hasproperty(p,name) && (ag[key]=getproperty(p,name)[:,i])
                end
                hasproperty(p,:weights) && (ag["genome_signature"]=bytes2hex(SHA.sha256(reinterpret(UInt8,vec(p.weights[:,:,i]))))[1:12])
                hasproperty(p,:plastic_marks) && (ag["plastic_norm"]=sqrt(sum(abs2,p.plastic_marks[:,:,i])))
                if metabolic!==nothing
                    ag["metabolism"]=Dict("allocation"=>metabolic.allocation[:,i],
                        "assimilation"=>metabolic.assimilation[i],"efficiency"=>metabolic.efficiency[i],"cooperation"=>metabolic.cooperation[i])
                    hasproperty(metabolic,:byproduct_fraction) && (ag["metabolism"]["byproduct_fraction"]=metabolic.byproduct_fraction[i])
                    for name in (:contributed_material,:contributed_work,:recent_contributions)
                        hasproperty(metabolic,name) && (ag[string(name)]=get(getproperty(metabolic,name),p.ids[i],0.0))
                    end
                    if metabolic.boundary!==nothing
                        ag["boundary_exposure"]=local_boundary_benefit(metabolic.boundary,Tuple(_neural_cell(view(p.positions,:,i),w.grid.size)))
                    end
                end
                push!(agents,ag)
            end
            if p isa Union{Ecology,NeuralEcology}
                for (key,value) in pairs(ecology_statistics(p)); metrics[string(key)]=value; end
            end
            metrics["population"]=length(p.ids)
            for name in (:births,:deaths,:blocked_births,:heat,:energy_input,:signal_work,:remodeling_work,:movement_work)
                hasproperty(p,name) && (metrics[string(name)]=getproperty(p,name))
            end
            hasproperty(p,:lineages) && (metrics["surviving_lineages"]=length(unique(p.lineages)))
            hasproperty(p,:generations) && (metrics["max_generation"]=maximum(p.generations;init=0))
            if hasproperty(p,:last_actions)
                actions=p.last_actions
                metrics["action_diversity"]=isempty(actions) ? 0.0 : sum(sum(abs2, Float64.(actions).-sum(Float64.(actions);dims=2)./size(actions,2);dims=2)./size(actions,2))
            end
        end
    end
    metrics["supplied_material"]=supplied; metrics["removed_material"]=removed
    metrics["material_error"]=metrics["material"]-a.initial_mass-supplied+removed
    if haskey(metrics,"resources")
        r=metrics["resources"]
        metrics["chemical_energy_accounted"]=arena_energy(a.processes)
        metrics["chemical_energy_error"]=metrics["chemical_energy_accounted"]-get(a.config,"_initial_energy",0.0)-r["supplied_potential"]+r["removed_potential"]
    end
    dep=Array(w.fields.deposit); metrics["deposit_occupancy"]=count(>(0.05f0),dep)/length(dep)
    metrics["deposit_range"]=[minimum(dep),maximum(dep)]
    Dict("schema"=>"interior-arena-live","schema_version"=>1,"generated_at_unix"=>time(),"step"=>a.steps,"time"=>w.time,"paused"=>a.paused,"stopped"=>a.stopped,
         "grid"=>Dict("shape"=>collect(dims),"extent"=>collect(w.grid.extent),"slice_z"=>k),"fields"=>fields,"agents"=>agents,"metrics"=>metrics)
end
function arena_command!(a::Arena,c)
    id=string(c["id"]); id==a.last_command && return nothing
    a.last_command=id
    action=string(c["action"])
    if action=="pause"; a.paused=true
    elseif action=="resume"; a.paused=false
    elseif action=="stop"; a.stopped=true
    elseif action=="checkpoint"; return :checkpoint
    elseif action=="set"
        values=c["values"]
        all(k->string(k) in ("dt","source_rate","sink_rate","snapshot_seconds","checkpoint_seconds"),keys(values)) || throw(ArgumentError("unsupported live setting"))
        all(v->v isa Real && isfinite(v) && v>=0,Base.values(values)) || throw(ArgumentError("settings must be finite nonnegative numbers"))
        haskey(values,"dt") && foreach(p->check_dt(a.world,p,values["dt"]),a.processes)
        for name in ("snapshot_seconds","checkpoint_seconds"); haskey(values,name) && values[name]<=0 && throw(ArgumentError("cadence must be positive")); end
        if any(p->p isa ResourcePool,a.processes) && any(k->haskey(values,k),("source_rate","sink_rate"))
            throw(ArgumentError("metabolic sources require per-species rates and finite budgets; use the profile"))
        end
        for p in a.processes
            if p isa Chemostat
                candidate=Chemostat(source_rate=get(values,"source_rate",p.source_rate),sink_rate=get(values,"sink_rate",p.sink_rate))
                check_dt(a.world,candidate,get(values,"dt",a.dt))
            end
        end
        for (name,value) in pairs(values)
            if name=="dt"; a.dt=Float64(value)
            elseif name in ("source_rate","sink_rate")
                for p in a.processes; p isa Chemostat && setproperty!(p,Symbol(name),Float64(value)); end
            else; a.config[name]=Float64(value)
            end
        end
    else; throw(ArgumentError("unknown action $action"))
    end
    nothing
end
# Atomic mkdir excludes competing writers. PID liveness uses POSIX kill(pid,0);
# an EPERM result still means alive. Lock contents are only our own protocol files.
function arena_pid_alive(pid::Integer)
    pid>0 || return true
    ccall(:kill,Cint,(Cint,Cint),pid,0)==0 || Base.Libc.errno()!=3
end
function arena_lock(dir;resume=false)
    mkpath(dir)
    marker=joinpath(dir,"arena-owner.json"); lock=joinpath(dir,".arena-writer")
    if isdir(lock)
        pidfile=joinpath(lock,"pid.json")
        isfile(pidfile) || throw(ArgumentError("arena lock has no PID; inspect it manually"))
        info=JSON3.read(read(pidfile,String),Dict{String,Any})
        get(info,"schema","")=="interior-arena-lock" || throw(ArgumentError("unrecognized arena lock"))
        arena_pid_alive(Int(info["pid"])) && throw(ArgumentError("arena already has a live writer PID $(info["pid"])"))
        readdir(lock)==["pid.json"] || throw(ArgumentError("arena lock contains unknown files; inspect it manually"))
        # Atomically claim stale directory before removing only protocol-owned contents.
        stale=lock*".stale-"*string(getpid())
        ispath(stale) && throw(ArgumentError("stale recovery destination exists"))
        mv(lock,stale)
        rm(joinpath(stale,"pid.json")); rm(stale)
    end
    try
        mkdir(lock)
    catch
        throw(ArgumentError("arena writer lock is already held"))
    end
    try
        arena_atomic(io->JSON3.write(io,Dict("schema"=>"interior-arena-lock","pid"=>getpid())),joinpath(lock,"pid.json"))
        if isfile(marker)
            info=JSON3.read(read(marker,String),Dict{String,Any})
            get(info,"schema","")=="interior-arena-output" || throw(ArgumentError("unrecognized arena ownership marker"))
            resume || throw(ArgumentError("existing arena output requires explicit resume=true / CLI --resume"))
        else
            isempty(filter(!=(".arena-writer"),readdir(dir))) || throw(ArgumentError("refusing nonempty unowned arena output directory"))
            arena_atomic(io->JSON3.write(io,Dict("schema"=>"interior-arena-output","version"=>1)),marker)
        end
    catch
        arena_unlock(lock); rethrow()
    end
    lock
end
function arena_unlock(lock)
    pidfile=joinpath(lock,"pid.json")
    if isfile(pidfile)
        info=JSON3.read(read(pidfile,String),Dict{String,Any})
        Int(info["pid"])==getpid() || return
        rm(pidfile)
    end
    isdir(lock) && isempty(readdir(lock)) && rm(lock)
    nothing
end
"""Run until a requested number of additional steps, stop command, or interrupt.
Owns latest.json, control.json and a fixed ring of checkpoint-N.jls files in `dir`.
Physics timestep and wall-clock publication cadence are independent."""
function run_arena!(a::Arena,dir;steps=nothing,snapshot_seconds=1.0,checkpoint_seconds=60.0,checkpoint_slots=3,batch=5,sleep_seconds=0.001,resume=false)
    all(isfinite,(snapshot_seconds,checkpoint_seconds,sleep_seconds)) && sleep_seconds>=0 && checkpoint_slots>=1 && batch>=1 && snapshot_seconds>0 && checkpoint_seconds>0 || throw(ArgumentError("invalid runner cadence"))
    steps===nothing || steps>=0 || throw(ArgumentError("steps must be nonnegative"))
    lock=arena_lock(dir;resume)
    a.config["snapshot_seconds"]=snapshot_seconds; a.config["checkpoint_seconds"]=checkpoint_seconds
    target=steps===nothing ? typemax(Int) : a.steps+steps
    lastview=-Inf; lastsave=time(); slot=0
    save()=(slot=mod(slot,checkpoint_slots)+1;checkpoint!(a,joinpath(dir,"checkpoint-$slot.jls"));arena_atomic(io->JSON3.write(io,Dict("file"=>"checkpoint-$slot.jls","step"=>a.steps)),joinpath(dir,"checkpoint.json")))
    try
        while !a.stopped && a.steps<target
            command=joinpath(dir,"command.json")
            if isfile(command)
                try
                    c=JSON3.read(read(command,String),Dict{String,Any})
                    if string(c["id"])!=a.last_command
                        result=arena_command!(a,c); result===:checkpoint && save()
                        arena_atomic(io->JSON3.write(io,Dict("id"=>a.last_command,"ok"=>true,"message"=>"applied")),joinpath(dir,"control.json"))
                    end
                catch e
                    arena_atomic(io->JSON3.write(io,Dict("id"=>a.last_command,"ok"=>false,"message"=>sprint(showerror,e))),joinpath(dir,"control.json"))
                end
            end
            for _ in 1:min(batch,target-a.steps); (a.paused || a.stopped) && break; arena_step!(a); end
            now=time()
            if now-lastview>=a.config["snapshot_seconds"]
                arena_atomic(io->JSON3.write(io,arena_snapshot(a)),joinpath(dir,"latest.json"));lastview=now
            end
            if now-lastsave>=a.config["checkpoint_seconds"]; save();lastsave=now; end
            sleep(a.paused ? 0.1 : sleep_seconds)
        end
    catch e
        e isa InterruptException || rethrow()
    finally
        try
            save()
            arena_atomic(io->JSON3.write(io,arena_snapshot(a)),joinpath(dir,"latest.json"))
        finally
            arena_unlock(lock)
        end
    end
    a
end
"""Build a CPU arena from a TOML-compatible configuration. Physical extent is an
independent parameter, allowing grid refinement without changing domain physics."""
function build_arena(config::AbstractDict=Dict{String,Any}())
    cfg=Dict{String,Any}(config)
    get(cfg,"model","neural")=="cellular" && return build_cellular_arena(cfg)
    grid=get(cfg,"grid",Dict()); dims=Tuple(Int.(get(grid,"shape",[96,96,1])))
    length(dims)==3 || throw(ArgumentError("grid.shape needs three dimensions"))
    extent=Tuple(get(grid,"extent",Float32.(dims)))
    model=get(cfg,"model","neural")
    model in ("neural","metabolic") || throw(ArgumentError("unknown arena model $model"))
    seed=Int(get(cfg,"seed",17)); w=World(Grid(dims...;extent);seed)
    phasecfg=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in get(cfg,"phase",Dict()))
    defaults=(mineral_pinning=0,nutrient_diffusion=0,interface_activity=0,activity_relaxation=0)
    phase=ImplicitPhaseSeparation(w;merge(Dict(pairs(defaults)),phasecfg)...)
    rng=MersenneTwister(seed+31)
    # Low modes create many independently interacting interfaces immediately.
    maxwave=Int(get(cfg,"initial_modes",max(2,round(Int,dims[1]/8))))
    modecount=Int(get(cfg,"initial_mode_count",12))
    initial_amplitude=Float32(get(cfg,"initial_amplitude",0.12f0))
    initial_mean=Float32(get(phasecfg,:mean,0f0))
    maxwave>=0 && modecount>0 && isfinite(initial_amplitude) && initial_amplitude>=0 && isfinite(initial_mean) || throw(ArgumentError("invalid arena initialization settings"))
    fill!(phase.phi,0)
    for _ in 1:modecount
        wave=ntuple(d->dims[d]==1 ? 0 : rand(rng,-maxwave:maxwave),3)
        amplitude=initial_amplitude*randn(rng,Float32); offset=2pi*rand(rng)
        for I in CartesianIndices(phase.phi)
            phase.phi[I]+=amplitude*cos(2pi*sum(wave[d]*(I[d]-1)/dims[d] for d in 1:3)+offset)
        end
    end
    phase.phi .-= sum(phase.phi)/length(phase.phi)
    initial_mean==0 || (phase.phi .+= initial_mean)
    phase_observe!(w,phase)
    mediumcfg=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in get(cfg,"medium",Dict()))
    medium=ActiveMedium(w,phase;mediumcfg...)
    ecfg=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in get(cfg,"ecology",Dict()))
    n=Int(pop!(ecfg,:population,1000))
    positions=rand(rng,Float32,3,n).*reshape(Float32.(collect(dims)),3,1)
    dims[3]==1 && (positions[3,:].=0)
    if model=="metabolic"
        rcfg=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in get(cfg,"resources",Dict()))
        names=Symbol.(pop!(rcfg,:names,["substrate_a","substrate_b","byproduct"]))
        spatial=Bool(pop!(rcfg,:spatial_niches,true))
        # Complementary periodic masks establish different local opportunities.
        # They are initial/boundary conditions, not an ecological conclusion.
        masks=[ones(Float32,dims) for _ in names]
        if spatial
            for I in CartesianIndices(masks[1])
                value=Float32(.5+.45*cos(2pi*(I[1]-1)/dims[1]))
                masks[1][I]=value
                length(masks)>1 && (masks[2][I]=1-value)
            end
        end
        resources=ResourcePool(w;names,medium,source_masks=masks,rcfg...)
        spatial && foreach(i->resources.fields[i].*=masks[i],eachindex(names))
        bcfg=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in get(cfg,"boundary",Dict()))
        boundary=BoundaryMaintenance(w,medium;bcfg...)
        mcfg=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in get(cfg,"metabolism",Dict()))
        if haskey(mcfg,:allocation)
            mcfg[:allocation]=repeat(reshape(Float32.(mcfg[:allocation]),:,1),1,n)
        end
        ecology=MetabolicEcology(positions,resources;boundary,medium,seed,mcfg...,ecfg...)
        return Arena(w,(resources,ecology,boundary,phase);dt=get(cfg,"dt",0.1),config=cfg)
    end
    ecology=NeuralEcology(positions;medium,seed,ecfg...)
    sourcecfg=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in get(cfg,"chemostat",Dict()))
    source=Chemostat(;sourcecfg...)
    Arena(w,(source,ecology,medium,phase);dt=get(cfg,"dt",0.1),config=cfg)
end
function check_dt(w::World,p::Chemostat,dt)
    isfinite(dt) && dt>0 && isfinite(Float32(dt)) && Float32(dt)>0 || throw(ArgumentError("invalid timestep"))
    validate(p)
    incoming=Float32(p.source_rate*Float32(dt))
    isfinite(incoming) && all(x->isfinite(x)&&x>=0&&isfinite(x+incoming),w.fields.nutrient) || throw(ArgumentError("invalid nutrient or source overflow"))
    nothing
end
