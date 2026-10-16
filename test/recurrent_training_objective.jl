using Test, Random
include("../studies/recurrent_training_objective.jl")
const RTO=RecurrentTrainingObjective
@testset "external objective and bounded-action derivatives" begin
    rng=MersenneTwister(713);A=3;P=2;N=2;T=3
    y=randn(rng,A+1+P,N,T);z=randn(rng,A,N,T);noise=randn(rng,A,N,T)
    target=randn(rng,N,T);adv=randn(rng,N,T);pred=randn(rng,P,N,T)
    valid=fill(true,N,T);valid[2,2]=false;pv=copy(valid);pv[1,3]=false
    f(v)=RTO.loss_gradient(v,z,noise,target,adv,pred,valid,pv;lower=[0.,0.,.02],upper=[.3,1.,.25],slope=[.5,.5,1.])
    r=f(y);eps=1e-6
    for i in eachindex(y)
        yp=copy(y);ym=copy(y);yp[i]+=eps;ym[i]-=eps
        @test r.gradient[i]≈(f(yp).loss-f(ym).loss)/(2eps) atol=2e-8 rtol=2e-6
    end
    @test all(iszero,r.gradient[:,2,2])
    @test all(iszero,r.gradient[A+2:end,1,3])
    @test isfinite(RTO.log_sech2(1000.))
    @test RTO.bounded_action(0.,.02,.25)≈.135
    out=zeros(2,3);reward=ones(2,3);cont=ones(2,3);cont[1,2]=0
    RTO.returns!(out,reward,cont,trues(2,3),[10.,10.];discount=.5)
    @test out[1,:]≈[1.5,1.,6.]
    @test out[2,:]≈[3.,4.,6.]
    # Slot reuse cannot transmit the new resident's rewards through a terminal.
    reward[1,3]=1000.;RTO.returns!(out,reward,cont,trues(2,3),[10.,10.];discount=.5)
    @test out[1,1:2]≈[1.5,1.]
end
@testset "native objective agrees with CPU oracle" begin
    rng=MersenneTwister(29);A=3;P=12;N=3;T=8
    y=randn(rng,Float32,A+1+P,N,T);z=randn(rng,Float32,A,N,T)
    noise=randn(rng,Float32,A,N,T);reward=randn(rng,Float32,N,T)
    cont=ones(Float32,N,T);cont[1,4]=0;valid=fill(true,N,T);valid[2,1]=false
    pred=randn(rng,Float32,P,N,T);pv=copy(valid);pv[1,4]=false
    bootstrap=randn(rng,Float32,N);lower=Float32[0,0,.02];upper=Float32[.3,1,.25]
    slope=fill(.5f0,A);scales=ones(Float32,P);w=RTO.workspace(y)
    RTO.evaluate!(w,y,z,noise,reward,cont,valid,bootstrap,pred,pv;count=sum(valid),lower,upper,slope,prediction_scale=scales)
    RTO.KernelAbstractions.synchronize(RTO.KernelAbstractions.get_backend(y))
    expected=zeros(Float32,N,T);RTO.returns!(expected,reward,cont,valid,bootstrap;discount=.995f0)
    @test w.targets≈expected
    adv=expected.-y[A+1,:,:]
    oracle=RTO.loss_gradient(y,z,noise,expected,adv,pred,valid,pv;lower,upper,slope,prediction_scale=scales)
    @test w.gradient≈oracle.gradient atol=2e-6 rtol=2e-6
    @test sum(w.parts[1,:,:])≈oracle.actor atol=2e-6 rtol=2e-6
    @test sum(w.parts[2,:,:])≈oracle.value atol=2e-6 rtol=2e-6
    @test sum(w.parts[3,:,:])≈oracle.prediction atol=2e-6 rtol=2e-6
    @test sum(w.parts[4,:,:])≈oracle.entropy atol=2e-6 rtol=2e-6
end
