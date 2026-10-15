using Interior, Test, KernelAbstractions, Serialization
include(joinpath(@__DIR__,"..","studies","constructive_world_batches.jl"))
include(joinpath(@__DIR__,"..","studies","constructive_batch_fixture.jl"))
@testset "Constructive ensemble whole-law parity and budgets" begin
    for flux in (false,true),strategy in (:axis,:materialized)
        fixtures=[constructive_batch_fixture(seed;flux,strategy) for seed in (11,29)]
        worlds=first.(fixtures);serial=last.(fixtures);batch=constructive_world_batch(serial)
        c=batch.cells;h=batch.habitat
        @test h.cells===c&&h.argument.status===c.environment.status&&h.feedstock.argument.status===c.environment.status
        @test h.argument.enzyme_bits===h.arrays.enzyme_bits
        @test Array(c.seeds)==UInt32[11,29]&&Array(c.habitat_seeds)==UInt32[1011,1029]
        initial=[Interior.device_habitat_audit(constructive_world_batch_snapshot(batch,i)) for i in 1:2]
        for step in 1:12
            for (w,p) in fixtures;advance!(w,p,.01);end
            advance!(worlds[1],batch,.01)
            for i in 1:2
                @test constructive_snapshot_equal(Interior.device_habitat_snapshot(serial[i]),constructive_world_batch_snapshot(batch,i))
            end
        end
        for i in 1:2
            x=constructive_world_batch_snapshot(batch,i);a=initial[i];b=Interior.device_habitat_audit(x)
            @test b.finite&&x.cells.environment.metadata[3]>0
            @test all(>(0),b.secreted)&&all(>(0),b.solubilized)&&all(>(0),b.bound)
            @test abs(b.material-a.material-b.supplied+a.supplied+b.removed-a.removed)<1e-4
            @test abs(b.potential-a.potential+b.work-a.work+b.heat-a.heat-b.supplied_potential+a.supplied_potential+b.removed_potential-a.removed_potential)<2e-4
        end
    end
end
@testset "Constructive batch strict checkpoint, restart and alias ownership" begin
    fixtures=[constructive_batch_fixture(seed;flux=true) for seed in (7,19)]
    worlds=first.(fixtures);batch=constructive_world_batch(last.(fixtures))
    for _ in 1:4;advance!(worlds[1],batch,.01);end
    worlds[1].time=.1;worlds[2].time=.2
    checkpoint=constructive_world_batch_checkpoint(worlds,batch)
    io=IOBuffer();serialize(io,checkpoint);seekstart(io);payload=deserialize(io)
    restored_worlds=[World(w.grid;seed=99) for w in worlds]
    for w in restored_worlds;fill!(w.fields.nutrient,0);end
    restored=constructive_world_batch_restore(restored_worlds,payload)
    @test [w.time for w in restored_worlds]==[.1,.2]
    @test [w.seed for w in restored_worlds]==[7,19]
    @test restored.habitat.cells===restored.cells
    @test restored.habitat.argument.status===restored.cells.environment.status
    @test restored.habitat.feedstock.argument.status===restored.cells.environment.status
    @test restored.habitat.argument.enzyme_bits===restored.habitat.arrays.enzyme_bits
    @test restored.cells.state.position.storage!==batch.cells.state.position.storage
    for _ in 1:4
        advance!(worlds[1],batch,.01);advance!(restored_worlds[1],restored,.01)
        @test all(constructive_snapshot_equal(constructive_world_batch_snapshot(batch,i),constructive_world_batch_snapshot(restored,i)) for i in 1:2)
    end
    @test_throws ArgumentError constructive_world_batch_restore(restored_worlds,merge(payload,(adapter_source=nothing,)))
    @test_throws ArgumentError constructive_world_batch_restore(restored_worlds,merge(payload,(member_indices=[1,1],)))
    badworld=World(worlds[1].grid)
    @test_throws ArgumentError constructive_world_batch_restore([badworld,restored_worlds[2]],payload)
    corrupted=deepcopy(payload);corrupted.members[1].arrays.bound[1,1]=-1
    @test_throws ArgumentError constructive_world_batch_restore(restored_worlds,corrupted)
end
@testset "Constructive failure isolation and checkpoint refusal" begin
    fixtures=[constructive_batch_fixture(seed) for seed in (31,43)]
    worlds=first.(fixtures);serial=last.(fixtures);serial[1].cells.environment.status[1]=8
    batch=constructive_world_batch(serial)
    for (w,p) in fixtures
        Interior.device_cellular_step!(w,p.cells,.01;check_status=false,anabolism=p.habitat)
    end
    Interior.device_cellular_step!(worlds[1],batch.cells,.01;check_status=false,anabolism=batch.habitat)
    for i in 1:2
        @test constructive_snapshot_equal(Interior.device_habitat_snapshot(serial[i];check_status=false),constructive_world_batch_snapshot(batch,i;check_status=false))
    end
    healthy=constructive_world_batch_snapshot(batch,2;check_status=false)
    @test healthy.cells.environment.status==[0]&&healthy.cells.environment.metadata[6]==1
    @test_throws ArgumentError constructive_world_batch_checkpoint(worlds,batch)
    healthy_checkpoint=constructive_world_batch_checkpoint(worlds,batch;members=[2])
    @test healthy_checkpoint.member_indices==[2]
    target=World(worlds[2].grid;seed=1);fill!(target.fields.nutrient,0)
    resumed=constructive_world_batch_restore([target],healthy_checkpoint)
    @test constructive_snapshot_equal(constructive_world_batch_snapshot(resumed,1),healthy)
    constructive_world_batch_step!(worlds,batch,.01;check_status=false)
    constructive_world_batch_step!([target],resumed,.01)
    @test worlds[1].time==0.0 && worlds[2].time==Float64(.01f0)
    @test constructive_snapshot_equal(constructive_world_batch_snapshot(batch,2;check_status=false),constructive_world_batch_snapshot(resumed,1))
end

@testset "constructive batch World and backend ownership guards" begin
    fixtures=[constructive_batch_fixture(i) for i in (71,72)]
    worlds=first.(fixtures);batch=constructive_world_batch(last.(fixtures))
    @test_throws ArgumentError constructive_world_batch(Any[])
    @test_throws ArgumentError constructive_world_batch(last.(fixtures);backend=:unsupported)
    @test_throws ArgumentError constructive_world_batch_step!([worlds[1],worlds[1]],batch,.01)
    worlds[2].time=.1
    constructive_world_batch_step!(worlds,batch,.01)
    @test worlds[2].time == .1 + Float64(.01f0)
    worlds[2].time=0.
    worlds[2].fields.nutrient[1]=1.
    @test_throws ArgumentError constructive_world_batch_checkpoint(worlds,batch)
end

@testset "World-isolated multilevel resource reduction" begin
    fixtures=[constructive_batch_fixture(seed;flux=true,shape=(6,6,6),extent=(6.,6.,6.)) for seed in (81,82)]
    worlds=first.(fixtures);batch=constructive_world_batch(last.(fixtures))
    @test length(Interior.device_resource_workspace(batch.cells).levels)==2
    @test batch_resource_validate(batch.cells)
    for step in 1:2
        for (w,p) in fixtures;Interior.device_cellular_step!(w,p.cells,.01;anabolism=p.habitat);end
        constructive_world_batch_step!(worlds,batch,.01)
        for i in 1:2
            @test constructive_snapshot_equal(Interior.device_habitat_snapshot(fixtures[i][2]),constructive_world_batch_snapshot(batch,i))
        end
    end
end
