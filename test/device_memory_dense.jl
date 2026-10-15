using Interior,Test
include(joinpath(@__DIR__,"../studies/device_memory_dense.jl"))
include(joinpath(@__DIR__,"../scripts/device_habitat_fixture.jl"))
const DMD=DeviceMemoryDense
@testset "direct dense specialized storage" begin
    source,p=constructive_device_fixture(CPU();support_strategy=:axis)
    for _ in 1:3;step!(source,p,.01);end
    checkpoint=Interior.device_habitat_checkpoint(source,p)
    target=World(source.grid,CPU(),NamedTuple(),NamedTuple(),0.,source.seed)
    q=DMD.restore(target,checkpoint)
    @test Base.mightalias(q.cells.environment.scratch,q.cells.environment.scale)
    @test Base.mightalias(q.cells.environment.scratch,q.habitat.feedstock.arrays.scratch.data)
    for _ in 1:8
        step!(source,p,.01);step!(target,q,.01)
        @test isequal(DMD.M.persistent_snapshot(p),DMD.M.persistent_snapshot(q))
    end
    final=Interior.device_habitat_checkpoint(target,q)
    world=World(source.grid,CPU(),NamedTuple(),NamedTuple(),0.,source.seed)
    restored=Interior.device_habitat_restore(world,final)
    for _ in 1:4
        step!(target,q,.01);step!(world,restored,.01)
        @test isequal(DMD.M.persistent_snapshot(q),DMD.M.persistent_snapshot(restored))
    end
    @test DMD.retained_bytes(q,a->a isa Array)<DMD.retained_bytes(p,a->a isa Array)
    bad=q.habitat.feedstock.arrays.fields
    bad[1,3]=1f0
    @test q.cells.environment.status[1]==4096
    @test_throws ArgumentError Interior.device_habitat_checkpoint(target,q)
end
