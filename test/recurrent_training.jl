using Test, Random, LinearAlgebra
include("../studies/recurrent_training.jl")
const RT=RecurrentTraining
@testset "actor and prediction both train recurrent representation" begin
    rng=MersenneTwister(9);D=4;H=5;A=3;P=2;B=2;T=4
    theta=RT.initialize(;D,H,A,P)
    x=randn(rng,Float32,D,B,T);h0=randn(rng,Float32,H,B)
    mask=fill(true,B,T);mask[1,3]=0
    valid=fill(true,B,T);pred=randn(rng,Float32,P,B,T)
    base=(;x,h0,memory_continuation=mask,z=randn(rng,Float32,A,B,T),
      noise=randn(rng,Float32,A,B,T),reward=randn(rng,Float32,B,T),return_continuation=copy(mask),
      valid,bootstrap=zeros(Float32,B),prediction=pred,predvalid=valid,count=B*T)
    t=RT.trainer(theta,B,T;prediction_scale=ones(Float32,P))
    before=deepcopy(theta);result=RT.train_step!(theta,t,base)
    @test isfinite(result.loss)
    @test norm(theta.Wh-before.Wh)>0
    @test norm(theta.Wx-before.Wx)>0
    # Isolate prediction and actor dY routes through exact generic backward.
    RT.R.forward!(t.cache,theta,x,h0,mask)
    for rows in (1:A,A+2:A+1+P)
        dy=zeros(Float32,A+1+P,B,T);dy[rows,:,:].=1
        RT.R.backward!(t.cache,theta,x,h0,mask,dy)
        @test norm(t.cache.g.Wh)>0
        @test norm(t.cache.g.Wx)>0
    end
    @test isempty(intersect([s.seed for s in RT.epoch_specs().training],[s.seed for s in RT.epoch_specs().heldout]))
    @test_throws ArgumentError RT.trainer(theta,B,T;sigma=NaN,prediction_scale=ones(P))
    @test_throws ArgumentError RT.trainer(theta,B,T;discount=1.01,prediction_scale=ones(P))
    @test_throws ArgumentError RT.trainer(theta,B,T;prediction_weight=-1,prediction_scale=ones(P))
    @test_throws ArgumentError RT.trainer(theta,B,T;prediction_scale=fill(1e-100,P))
    @test_throws ArgumentError RT.trainer(theta,B,T;lower=zeros(3),upper=zeros(3),prediction_scale=ones(P))
    @test_throws DimensionMismatch RT.trainer(theta,B,T;prediction_scale=ones(P+1))
end
@testset "composed frozen-target objective gradient into recurrent weights" begin
    rng=MersenneTwister(92);D=2;H=3;A=3;P=2;B=2;T=3
    theta=RT.R.initialize(rng,D,H,A+1+P;T=Float64)
    x=randn(rng,D,B,T);h0=randn(rng,H,B);mask=fill(true,B,T);mask[1,2]=0
    z=randn(rng,A,B,T);noise=randn(rng,A,B,T);targets=randn(rng,B,T)
    advantages=randn(rng,B,T);pred=randn(rng,P,B,T);valid=fill(true,B,T)
    cache=RT.R.allocate(theta,B,T)
    function objective()
        RT.R.forward!(cache,theta,x,h0,mask)
        RT.O.loss_gradient(cache.y,z,noise,targets,advantages,pred,valid,valid)
    end
    r=objective();RT.R.backward!(cache,theta,x,h0,mask,r.gradient)
    expected=copy(cache.g.Wh)
    for i in eachindex(theta.Wh)
        old=theta.Wh[i];eps=1e-6
        theta.Wh[i]=old+eps;plus=objective().loss
        theta.Wh[i]=old-eps;minus=objective().loss
        theta.Wh[i]=old
        @test expected[i]≈(plus-minus)/(2eps) atol=2e-8 rtol=2e-6
    end
end
