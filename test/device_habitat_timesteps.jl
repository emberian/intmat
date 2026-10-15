using Test, Interior
include(joinpath(@__DIR__,"..","studies","device_habitat_timesteps.jl"))
@testset "Fixed physical duration habitat timesteps" begin
    result=habitat_timestep_study(CPU();size=8,founders=2,capacity=32,reference_steps=16)
    rows=result.rows
    @test length(rows)==5
    @test all(r->r.duration==first(rows).duration,rows)
    @test all(r->r.dt<=r.declared_limit&&r.audit.finite,rows)
    @test first(rows).difference.free.maxabs==0
    @test first(rows).difference.position_max==0
    @test all(r->r.difference.identities_equal,rows)
    @test all(r->abs(r.material_error)<1e-4,rows)
    @test all(r->abs(r.potential_error)<1e-4,rows)
    @test all(r->r.simulated_time_per_wall_second==r.duration/r.seconds,rows)
end
