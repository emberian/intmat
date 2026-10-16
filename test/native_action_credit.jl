using Test, KernelAbstractions, Random, Statistics
include(joinpath(@__DIR__,"../studies/native_action_credit.jl"))
@testset "action-score credit mechanism" begin
    # Independent Gaussian quadrature by antithetic Monte Carlo: the score
    # estimator must agree with a pathwise derivative for a nonlinear payoff.
    rng=MersenneTwister(7301);z=randn(rng,100_000);noise=.3vcat(z,-z)
    x=.7;theta=.2;y=tanh.(theta*x.+noise)
    score=mean(y.*noise.*x./.3^2)
    exact=mean((1 .-y.^2).*x)
    @test isapprox(score,exact;rtol=.015)
    # The pre-tanh score has no extra tanh derivative, and payoff centering by
    # an action-independent baseline leaves its mean unchanged.
    @test isapprox(mean((y.-.4).*noise.*x./.3^2),score;atol=1e-12)
    learned=zeros(Float32,2,4,2);marks=fill(.13f0,2,4,5)
    perturb=Float32[.2 -.1;-.3 .4];rows=(1,3);h=2
    NativeActionCredit.install!(CPU(),64)(marks,learned,perturb,h,rows;ndrange=size(learned))
    @test all(marks[:,:,1:2].==.13f0)
    @test all(marks[:,:,4].==.13f0)
    @test marks[:,4,3]==perturb[:,1]
    @test marks[:,4,5]==perturb[:,2]
    @test all(marks[:,1:3,3].==0)
    trace=zero(learned);inputs=Float32[.1 .2 .3 1;.3 .2 .1 1]
    NativeActionCredit.eligible!(CPU(),64)(trace,inputs,perturb,4f0;ndrange=size(trace))
    @test trace[:,:,1]≈inputs.*perturb[:,1].*4
    @test trace[:,:,2]≈inputs.*perturb[:,2].*4
    NativeActionCredit.reinforce!(CPU(),64)(learned,trace,Float32[1,-1],0f0,.3f0;ndrange=size(learned))
    @test all(iszero,learned)
    NativeActionCredit.reinforce!(CPU(),64)(learned,trace,Float32[1,-1],.1f0,.3f0;ndrange=size(learned))
    @test learned[:,:,1]≈ 0.1 .*trace[:,:,1].*Float32[1,-1]
end
