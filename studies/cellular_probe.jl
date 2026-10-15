#!/usr/bin/env julia
# A bounded machinery probe with matched founders. No live arena is read or changed.
using Interior, Random, Statistics, JSON3, Serialization, SHA, Dates

function probe_config(;seed=113,founders=24,shape=(24,24,1),extent=(12.,12.,1.))
    volume=prod(Float64(extent[d])/shape[d] for d in 1:3 if shape[d]>1)
    Dict{String,Any}("model"=>"cellular","seed"=>seed,"dt"=>.1,
        "grid"=>Dict("shape"=>collect(shape),"extent"=>collect(extent)),
        "cellular"=>Dict("population"=>founders,"capacity"=>512),
        # ResourcePool stores voxel mass, not concentration: normalize both the
        # initial fields and per-voxel source rates by physical voxel volume.
        # Config values below are concentrations; the arena builder scales them.
        "resources"=>Dict("names"=>["precursor_a","precursor_b","waste"],
            "potentials"=>[4.,3.,0.],"initial"=>[1.5,1.5,0.],
            "diffusion"=>[.06,.06,.02],"source_rates"=>[.025,.025,0.],
            "source_budgets"=>[120.,120.,0.],"spatial_niches"=>true),
        "probe_cell_volume"=>volume)
end

function fresh_probe(cfg)
    arena=build_arena(cfg)
    w=arena.world;processes=arena.processes
    p=only(filter(p->p isa CellularEcology,processes));resources=p.resources
    backend="build_arena"
    # The study controls smooth complementary physical stripes explicitly,
    # independent of the product constructor's environmental defaults.
    for I in CartesianIndices(resources.fields[1])
        coordinate=(I[1]-.5)*w.grid.spacing[1]
        a=.5+.45cos(2pi*coordinate/w.grid.extent[1])
        resources.source_masks[1][I]=a;resources.source_masks[2][I]=1-a
        resources.source_masks[3][I]=0
        resources.fields[1][I]=Float32(cfg["resources"]["initial"][1]*cfg["probe_cell_volume"]*a)
        resources.fields[2][I]=Float32(cfg["resources"]["initial"][2]*cfg["probe_cell_volume"]*(1-a))
    end
    if arena!==nothing
        arena.initial_mass=Interior.arena_inventory(w,processes)
        arena.config["_initial_energy"]=Interior.arena_energy(processes)
    end
    (;world=w,processes,cells=p,resources,arena,backend,config=cfg)
end

function component_sizes(ids,edges;indices=false)
    n=length(ids);index=Dict(id=>i for (i,id) in enumerate(ids));parent=collect(1:n)
    function root(i)
        while parent[i]!=i;parent[i]=parent[parent[i]];i=parent[i];end
        i
    end
    for (a,b) in edges
        i=indices ? a : get(index,a,0);j=indices ? b : get(index,b,0)
        (i==0||j==0) && continue
        ri=root(i);rj=root(j);parent[ri]=rj
    end
    counts=Dict{Int,Int}()
    for i in 1:n;r=root(i);counts[r]=get(counts,r,0)+1;end
    sort!(collect(values(counts));rev=true)
end

function toroidal_extent(w,x)
    n=size(x,2)
    ntuple(3) do d
        (n<2||w.grid.size[d]==1) && return 0.
        coordinates=sort(mod.(x[d,:],Float64(w.grid.extent[d])))
        gaps=diff(vcat(coordinates,coordinates[1]+Float64(w.grid.extent[d])))
        Float64(w.grid.extent[d])-maximum(gaps)
    end
end

function total_chemical_energy(state)
    p=state.cells
    resource_potential(state.resources)+sum(p.energy)+p.work_spent+p.chemical_heat+
        sum(p.metabolites[s,i]*state.resources.potentials[s] for s in 1:2,i in eachindex(p.ids);init=0.)
end
probe_material(state)=resource_mass(state.resources)+sum(state.cells.biomass)+sum(state.cells.metabolites)+state.cells.inert_residual

function probe_metrics(state,baseline;previous=nothing)
    w=state.world;p=state.cells;r=state.resources;n=length(p.ids)
    contacts=cellular_pairs(w.grid,p.positions,cellular_radii(w,p))
    connected=component_sizes(p.ids,p.bonds.endpoints)
    touching=component_sizes(p.ids,contacts;indices=true)
    action_spread=n>1 ? vec(std(Float64.(p.last_actions);dims=2,corrected=false)) : zeros(11)
    hidden_spread=n>1 ? vec(std(Float64.(p.hidden);dims=2,corrected=false)) : zeros(size(p.hidden,1))
    temporal=Float64[]
    if previous!==nothing
        for (i,id) in enumerate(p.ids)
            haskey(previous.actions,id) && push!(temporal,sqrt(sum(abs2,p.last_actions[:,i].-previous.actions[id])))
        end
    end
    period=previous===nothing ? 0. : w.time-previous.time
    mechanics=p.mechanical
    merge(cellular_statistics(p),(time=w.time,contacts=length(contacts),
        bond_component_sizes=connected,contact_component_sizes=touching,
        largest_bond_component=isempty(connected) ? 0 : first(connected),
        toroidal_extent=toroidal_extent(w,p.positions),
        action_population_std=action_spread,hidden_population_std=hidden_spread,
        action_temporal_rms_change=isempty(temporal) ? 0. : sqrt(mean(abs2,temporal)),
        junction_flux=period>0 ? (p.junction_transferred-previous.junction_transferred)/period : 0.,
        supplied=copy(r.supplied),source_remaining=copy(r.source_remaining),
        mass_error=probe_material(state)-baseline.mass-sum(r.supplied)+sum(r.removed),
        chemical_energy_error=total_chemical_energy(state)-baseline.energy-sum(r.supplied.*r.potentials)+sum(r.removed.*r.potentials),
        mechanical_potential=mechanics===nothing ? 0. : mechanics.potential_after,
        mechanical_energy_residual=mechanics===nothing ? 0. : mechanics.energy_residual,
        mechanics_substeps=mechanics===nothing ? 0 : mechanics.substeps,
        finite=all(isfinite,p.positions)&&all(isfinite,p.biomass)&&all(isfinite,p.energy)&&all(a->all(isfinite,a),r.fields)))
end

function evolve_probe!(state,arm,directory;steps=200,capture_every=25)
    mkpath(directory);p=state.cells
    arm=="no_junction" && (p.parameters=merge(p.parameters,(junction_rate=0.,)))
    if arm=="no_adhesion"
        p.parameters=merge(p.parameters,(max_degree=0,))
        # Clamp the adhesion regulatory output off. Preserve all other founder
        # weights; this is an intervention, not an independently sampled lineage.
        row=size(p.hidden,1)+8;p.weights[row,:,:].=0;p.weights[row,end,:].=-6
    end
    baseline=(mass=probe_material(state),energy=total_chemical_energy(state))
    previous=nothing;history=Any[];wall_start=time()
    open(joinpath(directory,"metrics.jsonl"),"w") do io
        for step in 0:steps
            if step%capture_every==0||step==steps
                metric=merge(probe_metrics(state,baseline;previous),(step=step,arm=arm))
                push!(history,metric);JSON3.write(io,metric);write(io,'\n');flush(io)
                # Complete checkpoint includes world, chemical fields and budgets,
                # lineage, genotype, regulatory state, bonds and mechanics ledger.
                serialize(joinpath(directory,"snapshot-$(lpad(step,6,'0')).jls"),
                    (state=state,metrics=metric,baseline=baseline,sources=Interior.arena_sources()))
                previous=(time=state.world.time,junction_transferred=p.junction_transferred,
                    actions=Dict(id=>copy(p.last_actions[:,i]) for (i,id) in enumerate(p.ids)))
                metric.finite || error("nonfinite state in $arm at step $step")
            end
            step==steps && break
            arena_step!(state.arena)
            arm=="no_adhesion"&&!isempty(p.bonds.endpoints) && error("no-adhesion arm acquired bonds")
        end
    end
    (;arm,backend=state.backend,elapsed_seconds=time()-wall_start,initial=first(history),final=last(history),
        max_abs_mass_error=maximum(abs(row.mass_error) for row in history),
        max_abs_chemical_energy_error=maximum(abs(row.chemical_energy_error) for row in history),
        parameters=p.parameters,snapshots=length(history))
end

function contact_benchmark(;sizes=(100,300,1000),repeats=3)
    results=Any[]
    for n in sizes
        side=ceil(Int,sqrt(n));extent=Float64(side)
        grid=Grid(side,side;extent=(extent,extent,1.))
        x=zeros(3,n)
        for i in 1:n;x[1,i]=mod(i-1,side)+.5;x[2,i]=div(i-1,side)+.5;end
        ids=collect(1:n);radii=fill(.51,n);bonds=CellBonds()
        rng=MersenneTwister(9);motility=.03randn(rng,3,n);motility[3,:].=0.
        cellular_mechanics!(grid,ids,copy(x),radii,deepcopy(bonds),.1;motility)
        cellular_pairs(grid,x,radii)
        pairs=Any[];mechanics=Any[]
        for rep in 1:repeats
            pair_timing=@timed cellular_pairs(grid,x,radii)
            mechanics_timing=@timed cellular_mechanics!(grid,ids,copy(x),radii,deepcopy(bonds),.1;motility)
            push!(pairs,(seconds=pair_timing.time,bytes=pair_timing.bytes,pairs=length(pair_timing.value)))
            push!(mechanics,(seconds=mechanics_timing.time,bytes=mechanics_timing.bytes,
                substeps=mechanics_timing.value.substeps,bonds=mechanics_timing.value.bonds))
        end
        push!(results,(population=n,pairs=pairs,mechanics=mechanics,
            median_pair_seconds=median(r.seconds for r in pairs),
            median_mechanics_seconds=median(r.seconds for r in mechanics)))
    end
    results
end

function main(args=ARGS)
    root=length(args)>=1 ? abspath(args[1]) : joinpath(@__DIR__,"..","outputs","cellular-probe-"*Dates.format(now(),"yyyymmdd-HHMMSS"))
    ispath(root) && error("output exists; choose a fresh directory: $root")
    steps=length(args)>=2 ? parse(Int,args[2]) : 200
    steps>=1 || error("positive step count required")
    mkpath(root);cfg=probe_config()
    JSON3.write(joinpath(root,"configuration.json"),cfg)
    states=[fresh_probe(deepcopy(cfg)) for _ in 1:3]
    # Verify matching founder state before applying any intervention.
    all(s->s.cells.positions==states[1].cells.positions&&s.cells.weights==states[1].cells.weights,states) || error("founders differ")
    summaries=Any[]
    for (arm,state) in zip(("full","no_junction","no_adhesion"),states)
        println("Running ",arm," for ",steps," steps");flush(stdout)
        summary=evolve_probe!(state,arm,joinpath(root,arm);steps,capture_every=max(1,min(25,steps÷4)))
        push!(summaries,summary);JSON3.write(joinpath(root,"summary.json"),summaries)
    end
    println("Benchmarking contact mechanics at 100, 300, and 1000 cells");flush(stdout)
    benchmark=contact_benchmark();JSON3.write(joinpath(root,"contact-benchmark.json"),benchmark)
    metadata=(purpose="Bounded machinery and ablation probe, not an emergence or evolutionary assay",
        interpretation="Short runs diagnose mechanisms and numerical behavior. An absent behavior here does not exclude its emergence under other developmental histories, parameters, environments or longer selection.",
        model="Synthetic physiology; chemical ledger and coarse mechanical potential are separate",
        source_normalization="ResourcePool voxel masses and source rates scale with physical voxel volume; finite budgets are total mass",
        dt=cfg["dt"],steps=steps,arms=summaries,benchmark=benchmark,
        source_hashes=Interior.arena_sources(),julia=string(VERSION))
    JSON3.write(joinpath(root,"report.json"),metadata)
    println(root)
end

(abspath(PROGRAM_FILE)==(@__FILE__) || (isdefined(Main,:INTERIOR_RUNNING_STUDY) && Main.INTERIOR_RUNNING_STUDY=="cellular_probe")) && main()
