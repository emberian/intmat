#!/usr/bin/env julia
# Finite matched mechanism/evolution pilots. Never overwrites existing run files.
using Interior, Random, Statistics, JSON3, SHA, Serialization
include(joinpath(@__DIR__, "arena_experiments.jl"))

function metabolic_experiment(;seed=11,condition="niches",dims=(24,24,1),capacity=180,founders=24,uptake_rate=.7f0,lifespan=80f0)
    w=World(Grid(dims...;extent=Float32.(dims));seed)
    phase=ImplicitPhaseSeparation(w;mean=-.12f0,noise=0f0,mineral_pinning=0f0,nutrient_diffusion=0f0,activity_relaxation=0f0)
    copyto!(phase.phi,mode_initial_state(dims;seed)); Interior.phase_observe!(w,phase)
    medium=ActiveMedium(w,phase)
    boundary=BoundaryMaintenance(w,medium;decay=.02,work_per_mass=.6)
    masks=[ones(Float32,dims),ones(Float32,dims),zeros(Float32,dims)]
    if condition!="uniform"
        for I in CartesianIndices(masks[1])
            masks[1][I]=I[1]<=dims[1]÷2 ? 2f0 : 0f0
            masks[2][I]=2f0-masks[1][I]
        end
    end
    resources=ResourcePool(w;names=(:substrate_a,:substrate_b,:byproduct),potentials=(4,4,1),initial=(.5,.5,0),diffusion=(.03,.03,.03),source_rates=(.015,.015,0),source_budgets=(2000,2000,0),source_masks=masks,medium)
    rng=MersenneTwister(seed+1000)
    positions=rand(rng,Float32,3,founders).*reshape(Float32.(collect(dims)),3,1)
    dims[3]==1 && (positions[3,:].=0f0)
    allocation=fill(1f0/3,3,founders)
    # Same heterogeneous founders in each matched arm; not a claim of evolved specialization.
    for i in 1:founders; allocation[:,i].=isodd(i) ? Float32[.75,.15,.1] : Float32[.15,.75,.1]; end
    p=MetabolicEcology(positions,resources;boundary,allocation,medium,seed,capacity,lifespan,
        uptake_rate,trait_mutation=condition=="no_mutation" ? 0f0 : .04f0,
        byproduct_fraction=condition=="no_crossfeed" ? 0f0 : .2f0,
        cooperation=condition=="withhold" ? 0f0 : .35f0,boundary_rate=condition=="withhold" ? 0f0 : .04f0,
        mutation=condition=="no_mutation" ? 0f0 : .08f0,secretion_rate=0f0)
    arena=Arena(w,(resources,p,boundary,phase);dt=.2)
    arena,phase,medium,p,boundary
end

function metabolic_measure(a,phase,medium,p,b)
    c=p.controller; n=length(c.ids)
    avg(v)=isempty(v) ? 0.0 : mean(v)
    variance(v)=length(v)<2 ? 0.0 : var(v)
    cells=[Interior._neural_cell(c.positions[:,i],a.world.grid.size) for i in eachindex(c.ids)]
    left=[i for i in eachindex(c.ids) if cells[i][1]<=a.world.grid.size[1]÷2]
    right=setdiff(collect(eachindex(c.ids)),left)
    contributors=[i for i in eachindex(c.ids) if get(p.recent_contributions,c.ids[i],0.)>0]
    withholding=setdiff(collect(eachindex(c.ids)),contributors)
    exposure=[local_boundary_benefit(b,Tuple(cell)) for cell in cells]
    generalism=[-sum(x>0 ? x*log(x) : 0. for x in p.allocation[:,i])/log(size(p.allocation,1)) for i in 1:n]
    finite=all(arr->all(isfinite,arr),vcat(p.resources.fields,[p.allocation,p.efficiency,p.cooperation,b.film,c.energy,c.storage]))
    (time=a.world.time,population=n,births=c.births,deaths=c.deaths,finite=finite,
        material_error=Interior.arena_inventory(a.world,a.processes)-a.initial_mass-sum(p.resources.supplied)+sum(p.resources.removed),
        founder_effective_count=founder_effective_count(c.lineages),controller_variance=controller_diversity(c.weights),
        allocation_mean=[avg(p.allocation[s,:]) for s in axes(p.allocation,1)],allocation_variance=[variance(p.allocation[s,:]) for s in axes(p.allocation,1)],
        generalism_mean=avg(generalism),efficiency_mean=avg(p.efficiency),efficiency_variance=variance(p.efficiency),
        assimilation_mean=avg(p.assimilation),cooperation_mean=avg(p.cooperation),cooperation_variance=variance(p.cooperation),
        byproduct_fraction_mean=avg(p.byproduct_fraction),byproduct_fraction_variance=variance(p.byproduct_fraction),uptake_by_resource=copy(p.consumed),exported_by_resource=copy(p.resources.returned),resource_mass=[sum(Float64,f) for f in p.resources.fields],supplied=copy(p.resources.supplied),
        patch_population=[length(left),length(right)],patch_a_allocation=[avg(p.allocation[1,left]),avg(p.allocation[1,right])],
        inhabited_cells=length(unique(cells)),positive_phase_residents=count(cell->phase.phi[cell]>0,cells),
        boundary_exposure_mean=avg(exposure),contributor_population=length(contributors),withholding_population=length(withholding),contributor_exposure=avg(exposure[contributors]),withholding_exposure=avg(exposure[withholding]),maintenance_material=b.contributed_material,maintenance_work=b.spent_work,film_mass=boundary_material_mass(b),
        chemical_heat=p.chemical_heat,energy_input=c.energy_input,stored_energy=sum(Float64,c.energy;init=0.),controller_heat=c.heat,energy_error=sum(Float64,c.energy;init=0.)+c.heat-c.energy_input,movement_work=c.movement_work,signal_work=c.signal_work,remodeling_work=c.remodeling_work,reproduction_work=c.births*c.birth_energy,morphology=morphology_metrics(phase.phi))
end
function metabolic_sources()
    Interior.arena_sources()
end
function run_metabolic(;seed,condition,steps=400,capture_every=40,checkpoint=nothing,kwargs...)
    sources=metabolic_sources()
    a,phase,m,p,b=metabolic_experiment(;seed,condition,kwargs...)
    samples=[metabolic_measure(a,phase,m,p,b)]; first_step=1; prior_wall=0.; migrated_sources=nothing
    if checkpoint!==nothing && isfile(checkpoint)
        state=deserialize(checkpoint)
        mismatches=[k for (k,v) in state.sources if get(sources,k,nothing)!=v]
        if !isempty(mismatches)
            get(ENV,"METABOLIC_ALLOW_ARENA_METADATA_MIGRATION","0")=="1" && mismatches==["arena.jl"] || error("Checkpoint source differs; use a fresh output directory")
            migrated_sources=state.sources
        end
        state.seed==seed && state.condition==condition || error("Checkpoint comparison differs")
        a,phase,m,p,b=state.state; samples=state.samples; first_step=state.step+1; prior_wall=state.wall_seconds
    end
    started=time()
    for step in first_step:steps
        condition=="no_crossfeed" && fill!(p.byproduct_fraction,0f0)
        arena_step!(a)
        if step%capture_every==0 || step==steps
            push!(samples,metabolic_measure(a,phase,m,p,b))
            if checkpoint!==nothing
                serialize(checkpoint*".tmp",(seed=seed,condition=condition,sources=sources,state=(a,phase,m,p,b),samples=samples,step=step,wall_seconds=prior_wall+time()-started))
                mv(checkpoint*".tmp",checkpoint;force=true)
            end
            println((seed=seed,condition=condition,step=step,population=length(p.controller.ids)));flush(stdout)
        end
    end
    (seed=seed,condition=condition,steps=steps,dt=a.dt,wall_seconds=prior_wall+time()-started,model="finite-multisubstrate-funded-film-v1",settings=(capacity=p.controller.capacity,lifespan=p.controller.lifespan,uptake_rate=p.uptake_rate,founders=first(samples).population),sources=sources,migrated_checkpoint_sources=migrated_sources,source_unchanged=sources==metabolic_sources(),summary=summarize_windows(samples),samples=samples)
end
function metabolic_main()
    out=isempty(ARGS) ? "outputs/metabolic-$(Dates.format(Dates.now(), "yyyymmdd-HHMMSS"))" : first(ARGS)
    mkpath(out)
    steps=parse(Int,get(ENV,"METABOLIC_STEPS","400"))
    seeds=parse.(Int,split(get(ENV,"METABOLIC_SEEDS","11,29"),','))
    conditions=split(get(ENV,"METABOLIC_CONDITIONS","niches,uniform,withhold,no_mutation"),',')
    results=Any[]
    for seed in seeds,condition in conditions
        path=joinpath(out,"seed-$seed-$condition.json")
        if isfile(path); push!(results,JSON3.read(read(path,String))); continue; end
        r=run_metabolic(;seed,condition,steps,checkpoint=path*".checkpoint",uptake_rate=parse(Float32,get(ENV,"METABOLIC_UPTAKE","0.7")),lifespan=parse(Float32,get(ENV,"METABOLIC_LIFESPAN","80")))
        open(path,"w") do io; JSON3.write(io,r); end
        push!(results,r)
        open(joinpath(out,"summary.json"),"w") do io; JSON3.write(io,results); end
        println((seed=seed,condition=condition,wall_seconds=r.wall_seconds,summary=r.summary));flush(stdout)
    end
end
using Dates
abspath(PROGRAM_FILE)==(@__FILE__) && metabolic_main()
