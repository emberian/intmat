using Test, Random, Statistics
include(joinpath(@__DIR__,"..","studies","controller_capability_atlas.jl"))
using .ControllerCapabilityAtlas
@testset "Held-out controller capability measurement" begin
    rng=MersenneTwister(19)
    a=randn(rng,400,3);b=randn(rng,300,3)
    ya=2a[:,1].-a[:,2].+4;yb=2b[:,1].-b[:,2].+4
    @test ControllerCapabilityAtlas.readout(a,b,ya,yb)>.999
    # A training-constant feature must not acquire test-only information.
    @test ControllerCapabilityAtlas.readout(zeros(400,1),reshape(yb,:,1),ya,yb)<.02
    @test_throws ErrorException ControllerCapabilityAtlas.readout(a,b,ya,ones(300))
    @test ControllerCapabilityAtlas.sequence(seed=4)==ControllerCapabilityAtlas.sequence(seed=4)
    cfg=(gain=1.,adaptation=.3,plasticity=.002)
    u=ControllerCapabilityAtlas.sequence(samples=20,burn=8)
    a=ControllerCapabilityAtlas.record(cfg,u;burn=8,n=2)
    b=ControllerCapabilityAtlas.record(cfg,u;burn=8,n=2)
    @test a==b
    @test size(a.features)==(20,16,2)
    @test a.inputs==u[9:end]
    @test all(isfinite,a.differences)
    @test maximum(a.differences[1,:])>0
    @test_throws ErrorException ControllerCapabilityAtlas.record(cfg,u;burn=8,dt=.03,n=2)
end
