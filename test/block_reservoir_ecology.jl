using Test, Interior
include(joinpath(@__DIR__,"..","studies","block_reservoir_ecology.jl"))
using .BlockReservoirEcology
@testset "Conservative block reservoir transfers" begin
    blocks=[[1,2],[3,4]]
    fine=[1. 0.;3. 0.;0. 2.;0. 6.]
    coarse=zeros(Float32,2,2);restrict_mass!(coarse,fine,blocks)
    @test coarse==Float32[4 0;0 8]
    before=copy(coarse);after=Float32[1 2;3 0]
    reconcile_mass!(fine,before,after,blocks)
    @test minimum(fine)>=0
    @test vec(sum(fine;dims=1))≈[4.,2.]
    @test fine[:,1]≈[.25,.75,1.5,1.5]
    @test fine[:,2]≈[1.,1.,0.,0.]
    # Rounded-down representatives do not erase fine residual inventory.
    f=reshape([1.0 + 2.0^-30,0.],2,1);c=zeros(Float32,1,1)
    restrict_mass!(c,f,[[1,2]]);old=copy(c)
    reconcile_mass!(f,old,c,[[1,2]])
    @test sum(f)==1.0 + 2.0^-30
end
@testset "Complete coarse/fine closed ecology" begin
    b=block_fixture(fine=8,coarse=4,founders=2,capacity=32)
    initial=block_audit(b)
    for _ in 1:3;block_step!(b,.001);end
    final=block_audit(b)
    @test final.finite
    @test final.minimum_fine>=0
    @test abs(final.material-initial.material)<1e-4
    @test abs(final.potential+final.work+final.heat-initial.potential)<1e-4
end
@testset "Implicit block geometry and fine body porosity" begin
    for (fine,coarse) in ((8,4),(12,3),(7,7))
        blocks=BlockReservoirEcology.NestedBlocks(fine,coarse)
        @test sort(vcat([collect(x) for x in blocks]...))==collect(1:fine^3)
        @test sizeof(blocks)==3sizeof(Int)
    end
    b=block_fixture(fine=8,coarse=4,founders=2,capacity=32)
    w,p=BlockReservoirEcology.device_habitat_pilot_fixture(Interior.CPU();size=8,founders=2,capacity=32,support_strategy=:materialized)
    @test maximum(abs.(b.free.permeability.-p.cells.environment.permeability))<1e-6
end
