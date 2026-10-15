using Test, Interior
include(joinpath(@__DIR__,"..","studies","multirate_material_ecology.jl"))
using .MultirateMaterialEcology
@testset "Conservative multirate material clock" begin
    r=multirate_comparison(size=8,founders=2,capacity=32,fast_steps=8,ratios=(1,2,4,8))
    @test first(r.rows).all_arrays_exact
    @test all(x->x.audit.finite&&x.identities_equal,r.rows)
    @test all(x->x.outer_steps*x.ratio==x.fast_steps,r.rows)
    @test all(x->abs(x.material_error)<1e-4,r.rows)
    @test all(x->abs(x.potential_error)<1e-4,r.rows)
    @test all(x->x.fresh_catalyst_release_bound>=0,r.rows)
    @test_throws ErrorException multirate_fixture(ratio=0)
    @test_throws ErrorException multirate_fixture(size=8,ratio=1000,fast_dt=.001)
end
@test_throws ErrorException multirate_comparison(fast_steps=0)
