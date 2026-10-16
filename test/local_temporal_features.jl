using Test, KernelAbstractions
include(joinpath(@__DIR__,"../studies/local_temporal_features.jl"))
const TF=LocalTemporalFeatures

@testset "temporal coordinates use history, not feature gain" begin
    q=TF.coefficients(.2)
    z=TF.state(zeros(Float32,1,1),2)
    active=Int32[1,1];ids=Int32[11,21]
    inputs=Float32[1 .5 .5; -1 .5 .5]
    sample!()=TF.observe_kernel!(CPU(),64)(z,active,ids,inputs,1,q.fast,q.slow;ndrange=2)
    sample!()
    @test all(TF.contrast(z,i,k)==0 for i in 1:2,k in 1:3)
    # Same present inputs, mirrored histories: exactly matched feature norms.
    inputs[:,1].=0;sample!()
    a=[TF.contrast(z,1,k) for k in 1:3]
    b=[TF.contrast(z,2,k) for k in 1:3]
    @test a == -b
    @test sum(abs2,a)==sum(abs2,b)>0
    @test a[1]<0<b[1]
    # One frozen readout discriminates history at equal present/gain/cost.
    weights=Float32[.2,0,0]
    action_a=tanh(sum(weights.*a));action_b=tanh(sum(weights.*b))
    @test action_a == -action_b
    @test action_a != action_b
    @test inputs[1,:]==inputs[2,:]
    # Current-only copies at the same present cannot distinguish the twins.
    @test sum(weights.*inputs[1,:])==sum(weights.*inputs[2,:])
    # Reset only one history; other row is preserved and first sample is neutral.
    other=copy(z.fast[2,:]);ids[1]=12;sample!()
    @test all(TF.contrast(z,1,k)==0 for k in 1:3)
    @test z.ids[1]==12
    @test z.fast[2,:]!=other
    @test all(abs(TF.contrast(z,i,k))<=1 for i in 1:2,k in 1:3)
end

@testset "physical clock, independent state and replay" begin
    q=TF.coefficients(.2);half=TF.coefficients(.1)
    z=TF.state(zeros(Float32,1),1);w=deepcopy(z)
    TF.observe!(z,1,Int32(1),0f0,0f0,0f0,q.fast,q.slow)
    TF.observe!(w,1,Int32(1),0f0,0f0,0f0,half.fast,half.slow)
    for step in 1:10
        TF.observe!(z,1,Int32(1),1f0,-1f0,.5f0,q.fast,q.slow)
        for sub in 1:2
            TF.observe!(w,1,Int32(1),1f0,-1f0,.5f0,half.fast,half.slow)
        end
    end
    @test z.fast≈w.fast atol=2e-7
    @test z.slow≈w.slow atol=2e-7
    restored=deepcopy(z)
    for x in Float32[-.3,.7,-.1]
        TF.observe!(z,1,Int32(1),x,x,x,q.fast,q.slow)
        TF.observe!(restored,1,Int32(1),x,x,x,q.fast,q.slow)
    end
    @test z==restored
    restored.fast[1,1]=9f0
    @test z.fast[1,1]!=restored.fast[1,1]
    @test_throws ArgumentError TF.coefficients(0)
    @test_throws ArgumentError TF.coefficients(.2;fast=2.,slow=.4)
    @test_throws ArgumentError TF.coefficients(Inf)
end

@testset "explicit descendant contract and inactive slots" begin
    q=TF.coefficients(.2);z=TF.state(zeros(Float32,1),3)
    ids=Int32[3,4,5];active=Int32[1,1,0];parents=Int32[0,1,0]
    TF.observe!(z,1,ids[1],1f0,-1f0,.5f0,q.fast,q.slow)
    TF.observe!(z,1,ids[1],0f0,0f0,0f0,q.fast,q.slow)
    before=deepcopy(z)
    TF.inherit_kernel!(CPU(),64)(z,active,ids,parents,true;ndrange=3)
    @test z.fast[1,:]==before.fast[1,:]
    @test z.slow[1,:]==before.slow[1,:]
    @test z.fast[2,:]==z.fast[1,:]
    @test z.slow[2,:]==z.slow[1,:]
    @test z.ids==Int32[3,4,0]
    TF.observe!(z,2,ids[2],.8f0,.8f0,.8f0,q.fast,q.slow)
    @test z.fast[1,:]==before.fast[1,:]
    @test z.fast[1,:]!=z.fast[2,:]
    TF.inherit_kernel!(CPU(),64)(z,active,ids,parents,false;ndrange=3)
    @test z.ids[2]==0
    TF.observe!(z,2,ids[2],.8f0,.8f0,.8f0,q.fast,q.slow)
    @test all(TF.contrast(z,2,k)==0 for k in 1:3)
    # Reused/stale parent state must not leak into a newborn.
    ids[1]=9
    TF.inherit_kernel!(CPU(),64)(z,active,ids,parents,true;ndrange=3)
    @test z.ids[2]==0
    TF.observe_kernel!(CPU(),64)(z,active,ids,ones(Float32,3,3),1,q.fast,q.slow;ndrange=3)
    @test all(iszero,z.fast[3,:]) && z.ids[3]==0
end
