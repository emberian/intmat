using Interior, Test, KernelAbstractions
isdefined(@__MODULE__,:ConstructiveBatchCells) || include(joinpath(@__DIR__,"..","studies","constructive_world_batches.jl"))
isdefined(@__MODULE__,:constructive_batch_fixture) || include(joinpath(@__DIR__,"..","studies","constructive_batch_fixture.jl"))
@testset "Material sensing follows independent batch members" begin
    sense=(channels=(:raw_a,:bound_total,:catalyst_balance),scales=(.3,.02,.004))
    fixtures=[constructive_batch_fixture(seed;flux=true,energy_type=Float64,inventory_type=Float64,material_sensing=sense) for seed in (19,31)]
    # Deliberately unequal habitats detect accidental failure to slice the new
    # receptor descriptor along the world axis.
    fixtures[2][2].habitat.feedstock.arrays.fields[:,1].*=2
    batch=constructive_world_batch(last.(fixtures))
    for _ in 1:4
        for (w,p) in fixtures;advance!(w,p,.01);end
        advance!(fixtures[1][1],batch,.01)
        for i in 1:2
            @test constructive_snapshot_equal(Interior.device_habitat_snapshot(fixtures[i][2]),constructive_world_batch_snapshot(batch,i))
        end
    end
    cp=constructive_world_batch_checkpoint(first.(fixtures),batch)
    worlds=[World(w.grid;seed=w.seed) for (w,_) in fixtures]
    for w in worlds;fill!(w.fields.nutrient,0);end
    restored=constructive_world_batch_restore(worlds,cp)
    for i in 1:2
        @test constructive_snapshot_equal(constructive_world_batch_snapshot(restored,i),constructive_world_batch_snapshot(batch,i))
    end
end
