using Interior, Test, KernelAbstractions
isdefined(Interior,:DeviceCatalyticHabitat) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_habitat.jl"))
isdefined(Main,:constructive_device_fixture) || include(joinpath(@__DIR__,"..","scripts","device_habitat_fixture.jl"))
@testset "Constructive habitat shares dense and axis support law" begin
    wd,dense=constructive_device_fixture()
    wa,axis=constructive_device_fixture(;support_strategy=:axis)
    for _ in 1:20
        step!(wd,(dense,),.01);step!(wa,(axis,),.01)
    end
    xd=Interior.device_habitat_snapshot(dense);xa=Interior.device_habitat_snapshot(axis)
    for name in keys(xd.arrays)
        @test isequal(xd.arrays[name],xa.arrays[name])
    end
    for name in keys(xd.cells.environment)
        @test isequal(xd.cells.environment[name],xa.cells.environment[name])
    end
    for name in keys(xd.cells.state)
        name in (:sites,:support,:support_count) && continue
        @test isequal(xd.cells.state[name],xa.cells.state[name])
    end
    @test Interior.device_habitat_audit(xd)==Interior.device_habitat_audit(xa)
    checkpoint=Interior.device_habitat_checkpoint(wa,axis)
    restored_world=World(wa.grid;seed=911);fill!(restored_world.fields.nutrient,0)
    restored=Interior.device_habitat_restore(restored_world,checkpoint)
    @test Interior.device_habitat_snapshot(restored)==xa
end
