using Interior, Test, KernelAbstractions
include(joinpath(@__DIR__,"..","studies","learning_world_batches.jl"))
include(joinpath(@__DIR__,"..","studies","constructive_batch_fixture.jl"))
@testset "Paid learners share launches, retain independent worlds" begin
    sense=(channels=(:raw_a,:bound_total,:catalyst_balance),scales=(.3,.02,.004))
    fixtures=[constructive_batch_fixture(seed;flux=true,strategy=:materialized,energy_type=Float64,inventory_type=Float64,material_sensing=sense) for seed in (19,31)]
    worlds=first.(fixtures);ps=last.(fixtures)
    learners=[NativeLifetimeCredit.sidecar(p;dt=.01,decision_time=.02,reward_time=.04,seed=100+i,rate=i==1 ? .003 : -.003,critic=i==2) for (i,p) in enumerate(ps)]
    batch=learning_world_batch(ps,learners)
    initial_audits=[Interior.device_habitat_audit(x.physical) for x in learning_world_batch_snapshot(batch)]
    clocks=[World(w.grid;seed=w.seed) for w in worlds]
    for w in clocks;fill!(w.fields.nutrient,0);end
    for step in 1:10
        for (w,p,h) in zip(worlds,ps,learners)
            Interior.device_cellular_step!(w,p.cells,.01;anabolism=h)
            w.time+=Float64(Float32(.01))
        end
        learning_world_batch_step!(clocks,batch,.01)
        snapshots=learning_world_batch_snapshot(batch)
        for i in eachindex(ps)
            @test constructive_snapshot_equal(Interior.device_habitat_snapshot(ps[i]),snapshots[i].physical)
            @test isequal(learners[i].state,snapshots[i].state)
            @test learners[i].options==snapshots[i].options
            @test learners[i].steps==snapshots[i].steps
            @test worlds[i].time==clocks[i].time
        end
    end
    for (before,x) in zip(initial_audits,learning_world_batch_snapshot(batch))
        after=Interior.device_habitat_audit(x.physical)
        material=after.material-before.material-after.supplied+before.supplied+after.removed-before.removed
        potential=after.potential-before.potential+after.work-before.work+after.heat-before.heat-after.supplied_potential+before.supplied_potential+after.removed_potential-before.removed_potential
        @test after.finite
        @test abs(material)<1e-8
        @test abs(potential)<1e-8
    end
    @test sum(learners[1].state.paid)>0
    @test sum(learners[2].state.paid)>0
    # Both clocks are between reward windows; restoring only learned weights
    # would lose the accumulated trace, held drive and previous-value account.
    checkpoint=learning_world_batch_checkpoint(clocks,batch)
    restored_worlds=[World(w.grid;seed=w.seed) for w in clocks]
    for w in restored_worlds;fill!(w.fields.nutrient,0);end
    restored=learning_world_batch_restore(restored_worlds,checkpoint)
    @test learning_world_batch_snapshot(restored)==learning_world_batch_snapshot(batch)
    @test getproperty.(checkpoint.journal,:identity)==["world-1","world-2"]
    for _ in 1:3
        learning_world_batch_step!(clocks,batch,.01)
        learning_world_batch_step!(restored_worlds,restored,.01)
        @test learning_world_batch_snapshot(restored)==learning_world_batch_snapshot(batch)
        @test getproperty.(restored_worlds,:time)==getproperty.(clocks,:time)
    end
    learning_world_batch_advance!(clocks,batch,5,.01)
    for _ in 1:5;learning_world_batch_step!(restored_worlds,restored,.01);end
    @test learning_world_batch_snapshot(restored)==learning_world_batch_snapshot(batch)
    @test getproperty.(restored_worlds,:time)==getproperty.(clocks,:time)
    wrong=merge(checkpoint,(source=merge(checkpoint.source,(adapter="wrong",)),))
    @test_throws ArgumentError learning_world_batch_restore(restored_worlds,wrong)
    journal=copy(checkpoint.journal);journal[1]=merge(journal[1],(time=19.,))
    @test_throws ArgumentError learning_world_batch_restore(restored_worlds,merge(checkpoint,(;journal)))
    # Failure is local. No failed physical state or private learner can leak into
    # a selected healthy member's strict restart payload.
    frozen=deepcopy(learning_world_batch_snapshot(batch)[2].state)
    batch.physical.cells.environment.status.storage[1,2]=8
    before=getproperty.(clocks,:time)
    learning_world_batch_step!(clocks,batch,.01;check_status=false)
    @test clocks[1].time>before[1]
    @test clocks[2].time==before[2]
    @test learning_world_batch_snapshot(batch;check_status=false)[2].state==frozen
    @test_throws ArgumentError learning_world_batch_checkpoint(clocks,batch)
    healthy=learning_world_batch_checkpoint(clocks,batch;members=[1])
    w=World(clocks[1].grid;seed=19);fill!(w.fields.nutrient,0)
    selected=learning_world_batch_restore([w],healthy)
    @test only(learning_world_batch_snapshot(selected))==first(learning_world_batch_snapshot(batch;check_status=false))

end
