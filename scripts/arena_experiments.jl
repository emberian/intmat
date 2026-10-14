#!/usr/bin/env julia
# Matched causal comparisons; this file can also be included for its diagnostics.
using Interior, Random, JSON3, Statistics

"""Small seeded Fourier perturbations around a conserved mean; no privileged blob."""
function mode_initial_state(dims; seed=11, mean=-0.12f0, amplitude=0.08f0, modes=24, max_wave=minimum(filter(>(1),collect(dims)))<=16 ? 2 : 4)
    rng=MersenneTwister(seed); phi=fill(Float32(mean),dims)
    for _ in 1:modes
        wave=ntuple(d->dims[d]==1 ? 0 : rand(rng,-max_wave:max_wave),3)
        all(iszero,wave) && continue
        angle=2pi*rand(rng); strength=amplitude*randn(rng)/sqrt(modes)
        for I in CartesianIndices(phi)
            phi[I]+=strength*cos(2pi*sum(wave[d]*(I[d]-1)/dims[d] for d in 1:3)+angle)
        end
    end
    phi .+= Float32(mean)-sum(Float64,phi)/length(phi)
    phi
end

"""Periodic nearest-neighbor phase components, including all three active axes."""
function morphology_metrics(phi; threshold=0f0)
    mask=phi .> threshold; visited=falses(size(phi)); sizes=Int[]
    for start in CartesianIndices(mask)
        (!mask[start] || visited[start]) && continue
        pending=[start]; visited[start]=true; n=0
        while !isempty(pending)
            I=pop!(pending); n+=1
            for d in 1:3, direction in (-1,1)
                size(phi,d)==1 && continue
                J=CartesianIndex(ntuple(k->k==d ? mod1(I[k]+direction,size(phi,k)) : I[k],3))
                if mask[J] && !visited[J]; visited[J]=true; push!(pending,J); end
            end
        end
        push!(sizes,n)
    end
    occupied=sum(sizes;init=0)
    (occupancy=occupied/length(phi), components=length(sizes),
     largest_component_share=occupied==0 ? 0.0 : maximum(sizes)/occupied,
     component_size_cv=length(sizes)<2 ? 0.0 : std(sizes)/mean(sizes),
     phi_mean=sum(Float64,phi)/length(phi),phi_variance=var(vec(phi)))
end

function founder_effective_count(lineages)
    isempty(lineages) && return 0.0
    counts=Dict{Int,Int}()
    for id in lineages; counts[id]=get(counts,id,0)+1; end
    exp(-sum((n/length(lineages))*log(n/length(lineages)) for n in values(counts)))
end
function controller_diversity(weights)
    size(weights,3)<2 && return 0.0
    sum(var(weights;dims=3))/prod(size(weights)[1:2])
end

function summarize_windows(samples)
    mid=cld(length(samples),2); a=samples[mid]; b=last(samples); duration=b.time-a.time
    (late_duration=duration,late_births=b.births-a.births,late_deaths=b.deaths-a.deaths,
     late_population_change=b.population-a.population,
     late_birth_rate=duration==0 ? 0.0 : (b.births-a.births)/duration,
     late_death_rate=duration==0 ? 0.0 : (b.deaths-a.deaths)/duration,
     final=b)
end

function experiment_arena(;seed=11,dims=(32,32,1),condition="full",capacity=240,founders=32,
        source_rate=.015,sink_rate=.01,lifespan=35f0,uptake_rate=.5f0,dt=.2)
    w=World(Grid(dims...;extent=Float32.(dims));seed)
    phase=ImplicitPhaseSeparation(w;mean=-.12f0,noise=0f0,mineral_pinning=0f0,
        nutrient_diffusion=0f0,activity_relaxation=0f0)
    copyto!(phase.phi,mode_initial_state(dims;seed)); Interior.phase_observe!(w,phase)
    medium=ActiveMedium(w,phase)
    rng=MersenneTwister(seed+1000)
    positions=rand(rng,Float32,3,founders).*reshape(Float32.(collect(dims)),3,1)
    dims[3]==1 && (positions[3,:].=0f0)
    feedback=condition!="no_feedback"
    cells=NeuralEcology(positions;medium,seed,capacity,lifespan,uptake_rate,
        plasticity_rate=condition=="no_plasticity" ? 0f0 : .03f0,
        mutation=condition=="no_mutation" ? 0f0 : .08f0,
        secretion_rate=feedback ? .025f0 : 0f0,
        signal_rate=feedback ? .1f0 : 0f0,remodeling_rate=feedback ? .03f0 : 0f0)
    flux=Chemostat(;source_rate,sink_rate)
    arena=Arena(w,(flux,cells,medium,phase);dt)
    arena,phase,medium,cells,flux
end
function arena_measure(a,phase,m,c,flux)
    s=ecology_statistics(c); morphology=morphology_metrics(phase.phi)
    finite=all(x->all(isfinite,x),(phase.phi,m.signal,m.permeability,c.weights,c.plastic_marks,a.world.fields.nutrient))
    mass=Interior.arena_inventory(a.world,a.processes)
    merge((time=a.world.time,finite=finite,material_error=mass-a.initial_mass-flux.supplied+flux.removed,
        supplied=flux.supplied,removed=flux.removed,founder_effective_count=founder_effective_count(c.lineages),
        genotype_variance=controller_diversity(c.weights),plastic_rms=isempty(c.plastic_marks) ? 0.0 : sqrt(sum(abs2,c.plastic_marks)/length(c.plastic_marks)),
        signal_variance=var(vec(m.signal)),signal_max=maximum(m.signal),
        inhabited_cells=length(unique(Interior._neural_cell(c.positions[:,i],a.world.grid.size) for i in eachindex(c.ids)))),s,morphology)
end
function run_comparison(;seed,dims,condition,steps=800,capture_every=40,kwargs...)
    a,phase,m,c,flux=experiment_arena(;seed,dims,condition,kwargs...)
    samples=[arena_measure(a,phase,m,c,flux)]; started=time()
    for k in 1:steps
        arena_step!(a)
        (k%capture_every==0 || k==steps) && push!(samples,arena_measure(a,phase,m,c,flux))
    end
    (seed=seed,dims=dims,condition=condition,steps=steps,dt=a.dt,wall_seconds=time()-started,
     settings=(capacity=c.capacity,founders=first(samples).population,source_rate=flux.source_rate,sink_rate=flux.sink_rate,lifespan=c.lifespan,uptake_rate=c.uptake_rate),
     summary=summarize_windows(samples),samples=samples)
end
function main()
    out=length(ARGS)>0 ? ARGS[1] : "outputs/arena-experiments"; mkpath(out)
    steps=parse(Int,get(ENV,"ARENA_EXPERIMENT_STEPS","800")); results=Any[]
    for seed in (11,29), condition in ("full","no_feedback","no_plasticity","no_mutation")
        r=run_comparison(;seed,dims=(32,32,1),condition,steps)
        push!(results,r)
        open(joinpath(out,"seed-$seed-$condition.json"),"w") do io; JSON3.write(io,r); end
        println((seed=seed,condition=condition,summary=r.summary,wall_seconds=r.wall_seconds)); flush(stdout)
    end
    r=run_comparison(;seed=11,dims=(16,16,16),condition="full",steps,capacity=160,founders=24)
    push!(results,r)
    open(joinpath(out,"small-3d.json"),"w") do io; JSON3.write(io,r); end
    open(joinpath(out,"summary.json"),"w") do io; JSON3.write(io,results); end
    println((small3d=r.summary,wall_seconds=r.wall_seconds))
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
