# Matched initial genotypes with closed vs finite-budget resource-flux arms.
# Initialization only is in the fixture; every numerical stage is Interior's.
using Interior, JSON3, Serialization
isdefined(Main,:ConstructiveBatchCells)||include(joinpath(@__DIR__,"..","studies","constructive_world_batches.jl"))
isdefined(Main,:constructive_batch_fixture)||include(joinpath(@__DIR__,"..","studies","constructive_batch_fixture.jl"))
function study_constructive_world_batches(out;steps=240,dt=.01,every=40,seeds=(11,29))
    ispath(out)&&throw(ArgumentError("fresh output directory required"))
    steps>0&&every>0||throw(ArgumentError("positive steps/sample interval required"))
    labels=[(seed=seed,arm=arm,budget_scale=arm=="closed" ? 0. : 1.) for seed in seeds for arm in ("closed","finite_flux")]
    serial=[constructive_batch_fixture(x.seed;flux=true,budget_scale=x.budget_scale,capacity=32,mechanics=true) for x in labels]
    worlds=[World(w.grid,CPU(),NamedTuple(),NamedTuple(),0.,w.seed) for (w,_) in serial]
    batch=constructive_world_batch(last.(serial));mkdir(out)
    initial=[Interior.device_habitat_audit(constructive_world_batch_snapshot(batch,i)) for i in eachindex(labels)]
    observations=Any[];maximum_material=0.;maximum_potential=0.;matched=true;restarted=false
    function observe(step)
        rows=map(eachindex(labels)) do i
            snapshot=constructive_world_batch_snapshot(batch,i)
            reference=Interior.device_habitat_snapshot(serial[i][2])
            constructive_snapshot_equal(snapshot,reference)||error("complete constructive parity failed at step=$step world=$i")
            worlds[i].time==serial[i][1].time||error("World clock parity failed")
            a=initial[i];b=Interior.device_habitat_audit(snapshot)
            material=b.material-a.material-b.supplied+a.supplied+b.removed-a.removed
            potential=b.potential-a.potential+b.work-a.work+b.heat-a.heat-b.supplied_potential+a.supplied_potential+b.removed_potential-a.removed_potential
            maximum_material=max(maximum_material,abs(material));maximum_potential=max(maximum_potential,abs(potential))
            active=findall(==(1),snapshot.cells.state.active)
            lineages=length(unique(snapshot.cells.state.lineages[active]))
            (label=labels[i],step=step,time=worlds[i].time,population=b.population,lineages=lineages,
                births=Int(snapshot.cells.environment.metadata[3]),deaths=Int(snapshot.cells.environment.metadata[4]),
                blocked_births=Int(snapshot.cells.environment.metadata[5]),work=b.work,heat=b.heat,supplied=b.supplied,
                secreted=b.secreted,solubilized=b.solubilized,bound=b.bound,enzyme=b.enzyme,
                material_error=material,potential_work_heat_error=potential,finite=b.finite)
        end
        push!(observations,rows)
    end
    observe(0)
    elapsed=@elapsed for step in 1:steps
        for (w,p) in serial;step!(w,p,dt);end
        constructive_world_batch_step!(worlds,batch,dt)
        if step==div(steps,2)
            payload=constructive_world_batch_checkpoint(worlds,batch)
            open(joinpath(out,"midpoint.jls"),"w") do io;serialize(io,payload);end
            payload=open(deserialize,joinpath(out,"midpoint.jls"))
            restored_worlds=[World(w.grid,CPU(),NamedTuple(),NamedTuple(),0.,w.seed) for w in worlds]
            batch=constructive_world_batch_restore(restored_worlds,payload);worlds=restored_worlds;restarted=true
        end
        (step%every==0||step==steps)&&observe(step)
    end
    maximum_material<1e-3&&maximum_potential<2e-3||error("constructive budget drift gate failed")
    result=(model="existing DeviceConstructiveEcology",backend="CPU-KA",source=CONSTRUCTIVE_BATCH_SOURCE,
        geometry=(shape=(6,6,1),extent=(6,6,1),exchange_radius=1.,capacity_per_world=32,mechanical_substeps=16),
        steps=steps,dt=dt,worlds=length(labels),all_complete_snapshots_exact=matched,checkpoint_restart=restarted,
        maximum_material_error=maximum_material,maximum_potential_work_heat_error=maximum_potential,
        elapsed_with_serial_reference_and_audits=elapsed,observations=observations,
        interpretation="bounded matched resource-budget study; shared mutation streams within seed-matched arms; no adaptation or long-run diversity claim")
    open(joinpath(out,"result.json"),"w") do io;JSON3.pretty(io,result);println(io);end
    result
end
if abspath(PROGRAM_FILE)==@__FILE__
    out=isempty(ARGS) ? "outputs/constructive-world-batches-study" : ARGS[1]
    result=study_constructive_world_batches(out;steps=parse(Int,get(ENV,"CONSTRUCTIVE_BATCH_STEPS","240")))
    println((worlds=result.worlds,steps=result.steps,exact=result.all_complete_snapshots_exact,restarted=result.checkpoint_restart,
        material=result.maximum_material_error,potential=result.maximum_potential_work_heat_error))
end
