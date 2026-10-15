#!/usr/bin/env julia
# Fresh study only; no live arena access or edits to core physics.
using Interior, Random, Statistics, JSON3, Serialization, TOML, SHA
include(joinpath(@__DIR__, "..", "scripts", "metabolic_experiments.jl"))

const TRAITS=(:assimilation,:efficiency,:cooperation,:byproduct_fraction)
function genotype(p)
    (;weights=copy(p.controller.weights),allocation=copy(p.allocation),
      (k=>copy(getfield(p,k)) for k in TRAITS)...)
end
function install!(p,g,indices,slots)
    p.controller.weights[:,:,slots].=g.weights[:,:,indices]
    p.allocation[:,slots].=g.allocation[:,indices]
    for k in TRAITS; getfield(p,k)[slots].=getfield(g,k)[indices]; end
end
function fresh(cfg,seed,environment,n;mutation=true,plasticity=.03)
    a,phase,m,p,b=metabolic_experiment(;seed,condition=environment,founders=n,
        capacity=cfg["capacity"],uptake_rate=cfg["uptake_rate"],lifespan=cfg["lifespan"])
    a.dt=Float32(cfg["dt"])
    p.resources.source_remaining .= [cfg["source_budget"],cfg["source_budget"],0.]
    p.controller.plasticity_rate=Float32(plasticity)
    if !mutation; p.controller.mutation=0; p.trait_mutation=0; end
    a,phase,m,p,b
end
function observation(a,phase,m,p,b;group_births=nothing,group_deaths=nothing)
    c=p.controller
    q=metabolic_measure(a,phase,m,p,b)
    merge(q,(max_generation=isempty(c.ids) ? 0 : maximum(c.generations),
             source_remaining=copy(p.resources.source_remaining),
             chemical_potential=resource_potential(p.resources),
             living_groups=[count(==(i),c.lineages) for i in 1:2],
             group_births=group_births,group_deaths=group_deaths))
end
# Explicit total energy ledger uses exact initial potential + initial controller energy.
function measured(state,initial_energy;kwargs...)
    a,phase,m,p,b=state; q=observation(state...;kwargs...); c=p.controller
    accounted=resource_potential(p.resources)+sum(Float64,c.energy;init=0.)+c.heat+p.chemical_heat
    incoming=sum(p.resources.supplied.*Float64.(p.resources.potentials))
    outgoing=sum(p.resources.removed.*Float64.(p.resources.potentials))
    merge(q,(chemical_energy_error=accounted-initial_energy-incoming+outgoing,))
end
function evolve!(state,cfg,horizon,dir;groups=false)
    a,phase,m,p,b=state; c=p.controller; mkpath(dir)
    initial_energy=resource_potential(p.resources)+sum(Float64,c.energy;init=0.)
    births=zeros(Int,2); deaths=zeros(Int,2); samples=Any[]
    steps=round(Int,horizon/a.dt)
    open(joinpath(dir,"history.jsonl"),"w") do io
        for step in 0:steps
            if step%cfg["capture_every"]==0 || step==steps
                q=measured(state,initial_energy;group_births=copy(births),group_deaths=copy(deaths))
                push!(samples,q); JSON3.write(io,q); write(io,'\n'); flush(io)
                open(joinpath(dir,"lineages.jsonl"),"a") do li
                    JSON3.write(li,(time=a.world.time,ids=c.ids,parents=c.parents,lineages=c.lineages,generations=c.generations));write(li,'\n')
                end
            end
            step==steps && break
            oldnext=c.next_id
            oldgroups=groups ? Dict(zip(c.ids,c.lineages)) : nothing
            arena_step!(a)
            if groups
                live=Set(c.ids)
                for (id,g) in oldgroups; id in live || (deaths[g]+=1); end
                for i in eachindex(c.ids); c.ids[i]>=oldnext && (births[c.lineages[i]]+=1); end
            end
        end
    end
    all(s->s.finite,samples) || error("nonfinite study state")
    minimum(p.resources.source_remaining[1:2])>0 || error("funding exhausted")
    serialize(joinpath(dir,"terminal.jls"),(sources=Interior.arena_sources(),state=state))
    (samples=samples,group_births=births,group_deaths=deaths)
end
function assay(cfg,founder,descendant,seed,env,rep,plasticity,dir)
    n=cfg["assay_per_group"]; rng=MersenneTwister(seed+100000rep)
    state=fresh(cfg,seed+100000rep,env,2n;mutation=false,plasticity); p=state[4];c=p.controller
    c.capacity=cfg["assay_capacity"]
    # Interleave labels after a shared random permutation to avoid ID-priority effects.
    slots=randperm(rng,2n); ds=slots[1:n]; fs=slots[n+1:end]
    install!(p,descendant,rand(rng,axes(descendant.weights,3),n),ds)
    install!(p,founder,rand(rng,axes(founder.weights,3),n),fs)
    c.lineages[ds].=1;c.lineages[fs].=2
    @assert all(iszero,c.hidden) && all(iszero,c.plastic_marks) && all(iszero,c.age)
    @assert all(==(.4f0),c.storage) && all(==(1f0),c.energy)
    run=evolve!(state,cfg,cfg["assay_time"],dir;groups=true)
    counts=[count(==(i),c.lineages) for i in 1:2]
    (environment=env,replicate=rep,plasticity=plasticity,living=counts,
     births=run.group_births,deaths=run.group_deaths,
     log_frequency_change=log((counts[1]+.5)/(counts[2]+.5)),
     both_extinct=sum(counts)==0,descendant_extinct=counts[1]==0,founder_extinct=counts[2]==0,
     final=last(run.samples))
end
function main()
    cfgpath=length(ARGS)>1 ? ARGS[2] : joinpath(@__DIR__,"..","experiments","alife","adaptation.toml")
    cfg=TOML.parsefile(cfgpath);out=abspath(first(ARGS));ispath(out)&&!isempty(readdir(out))&&error("Use fresh output directory")
    mkpath(out); cp(cfgpath,joinpath(out,"protocol.toml")); cp(@__FILE__,joinpath(out,"study-source.jl"))
    # Max source mass per substrate = rate * grid cells * horizon, since masks sum to grid size.
    @assert cfg["source_budget"]>2*.015*24*24*max(cfg["selection_time"],cfg["assay_time"])
    results=Any[]
    for seed in cfg["seeds"], env in cfg["environments"], mutation in (true,false)
        label="seed-$seed-$env-$(mutation ? "mutation" : "no_mutation")";dir=joinpath(out,label)
        state=fresh(cfg,seed,env,cfg["founders"];mutation); founder=genotype(state[4]);mkpath(dir)
        serialize(joinpath(dir,"founders.jls"),founder)
        selected=evolve!(state,cfg,cfg["selection_time"],joinpath(dir,"selection"))
        descendant=genotype(state[4]); assays=Any[]
        if size(descendant.weights,3)>0
            for garden in cfg["environments"],rep in 1:cfg["position_replicates"],plasticity in cfg["assay_plasticity"]
                push!(assays,assay(cfg,founder,descendant,seed,garden,rep,plasticity,
                    joinpath(dir,"assay-$garden-$rep-plasticity-$plasticity")))
            end
        end
        r=(seed=seed,selection_environment=env,mutation=mutation,
           selection_final=last(selected.samples),selection_extinct=size(descendant.weights,3)==0,assays=assays)
        open(joinpath(dir,"result.json"),"w") do io;JSON3.write(io,r);end
        push!(results,r);open(joinpath(out,"results.json"),"w") do io;JSON3.write(io,results);end
        println((label=label,population=size(descendant.weights,3),assays=length(assays)));flush(stdout)
    end
    # Life-only null: identical founder pool sampled independently into both groups.
    nulls=Any[]
    for seed in cfg["seeds"],env in cfg["environments"],rep in 1:cfg["position_replicates"]
        state=fresh(cfg,seed,env,cfg["founders"]);g=genotype(state[4])
        push!(nulls,merge((seed=seed,),assay(cfg,g,g,seed,env,rep,0.,joinpath(out,"life-only-$seed-$env-$rep"))))
    end
    open(joinpath(out,"life-only.json"),"w") do io;JSON3.write(io,nulls);end
    open(joinpath(out,"provenance.json"),"w") do io;JSON3.write(io,(sources=Interior.arena_sources(),julia=string(VERSION),script_sha256=bytes2hex(sha256(read(@__FILE__)))));end
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
