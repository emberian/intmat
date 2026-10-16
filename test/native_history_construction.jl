using Test, KernelAbstractions
include(joinpath(@__DIR__,"../studies/native_history_construction.jl"))
const HC=NativeHistoryConstruction

@testset "occluded cue gives norm-matched causal history control" begin
    n=2;inputs=zeros(Float32,n,22);inputs[:,18].=Float32[.4,-.4]
    s=(active=Int32[1,1],ids=Int32[11,12],inputs=inputs,energy=ones(n),metabolites=zeros(Float32,n,2))
    q=(stride=10,cue_end=20,decision=50,hidden=1,af=HC.T.coefficients(.2).fast,as=HC.T.coefficients(.2).slow,
       sigma=.5f0,efficiency=.65f0,potentials=(4f0,4f0))
    function runmode(mode;swap=false)
        z=(learned=zeros(Float32,n,2,1),critic=zeros(Float32,n,8),trace=zeros(Float32,n,2,1),drive=zeros(Float32,n,1),
            action_feature=zeros(Float32,n),signed_feature=zeros(Float32,n),features=zeros(Float32,n,8),
            prediction=zeros(n),previous=zeros(n),ids=zeros(Int32,n))
        z.learned[:,2,1].=.3f0
        m=HC.T.state(inputs,n)
        for step in 0:10:50
            HC.observe_decide!(CPU(),64)(s,Int32[0],z,m,q,step,zeros(Float32,n),mode,NaN32,swap;ndrange=n)
        end
        z
    end
    signed=runmode(0);unsigned=runmode(1);swapped=runmode(0;swap=true)
    @test signed.action_feature[1]==-signed.action_feature[2]
    @test unsigned.action_feature==abs.(signed.action_feature)
    @test unsigned.action_feature.^2==signed.action_feature.^2
    @test signed.drive[1]==-signed.drive[2]
    @test unsigned.drive[1]==unsigned.drive[2]
    @test swapped.drive==-signed.drive
    @test inputs[:,18]==Float32[.4,-.4] # no mutation of native receptors
    @test all(iszero,signed.trace) # deterministic probe has no score innovation
end

@testset "allocation changes composition at fixed total request" begin
    s=(active=Int32[1,1],mass=Float32[1.6,1.6],actions=zeros(Float32,2,11))
    a=(expression=zeros(Float32,2,3),requests=zeros(Float64,2,2))
    z=(drive=reshape(Float32[-1,1],2,1),)
    q=(secretion_rate=.12f0,allocation_exponent=1f0)
    HC.actuate!(CPU(),64)(a,s,Int32[0],z,q,.02f0,true;ndrange=2)
    @test sum(a.requests[1,:])≈sum(a.requests[2,:]) atol=1e-9
    @test a.requests[1,1]<a.requests[1,2]
    @test a.requests[2,1]>a.requests[2,2]
    @test s.actions[:,11]==fill(-.6f0,2)
    HC.actuate!(CPU(),64)(a,s,Int32[0],z,q,.02f0,false;ndrange=2)
    @test all(iszero,a.requests)
end
