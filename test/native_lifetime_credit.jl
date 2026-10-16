using Test, KernelAbstractions
include(joinpath(@__DIR__,"../studies/native_lifetime_credit.jl"))
@testset "paid lifetime credit boundaries" begin
    s=(active=Int32[1],ids=Int32[17],energy=[1.],metabolites=zeros(1,2),
       weights=reshape(Float32[.2,.1],1,2,1),marks=zeros(Float32,1,2,1),inputs=reshape(Float32[.3,1],1,2),
       ledger=zeros(Float32,1,8),ledger_low=zeros(Float32,1,8))
    z=(ids=Int32[0],ready=Int32[0],learned=zeros(Float32,1,2,1),trace=zeros(Float32,1,2,1),drive=zeros(Float32,1,1),
       baseline=zeros(1),previous=zeros(1),paid=zeros(1),last_reward=zeros(1))
    q=(efficiency=.5,potentials=(4.,4.));status=Int32[0]
    choose(decision,clock)=NativeLifetimeCredit.choose!(CPU(),64)(s,status,q,z,(1,),0,decision,.15f0,UInt32(71),UInt32(clock);ndrange=1)
    choose(true,1);drive=copy(z.drive);trace=copy(z.trace)
    @test z.ids==s.ids
    @test z.previous==[1.]
    s.inputs[1,1]=.9f0;choose(false,1)
    @test z.drive==drive
    @test z.trace==trace
    @test sum((s.weights.+s.marks).*reshape(s.inputs,1,2,1))≈drive[1] atol=1e-7
    choose(true,2);@test z.trace!=trace
    # Positive return, paid mark update, exact native work debit.
    z.ready[1]=1;z.previous[1]=.8;fill!(z.trace,1f0);e=s.energy[1]
    NativeLifetimeCredit.reward!(CPU(),64)(s,status,q,z,.1f0,.3f0,.001,true;ndrange=1)
    @test all(z.learned.>0)
    @test e-s.energy[1]==z.paid[1]
    @test Float64(s.ledger[1,1])+s.ledger_low[1,1]≈z.paid[1] atol=1e-18
    @test all(iszero,z.trace)
    @test z.previous[1]==e # next return includes learning cost
    # Unaffordable updates do not alter marks or invent credit.
    learned=copy(z.learned);s.energy[1]=0.;z.previous[1]=-1.;fill!(z.trace,1f0)
    NativeLifetimeCredit.reward!(CPU(),64)(s,status,q,z,.1f0,.3f0,.001,true;ndrange=1)
    @test z.learned==learned
    # A reused slot never inherits the previous resident's private learning.
    s.ids[1]=18;choose(false,2)
    @test all(iszero,z.learned)
    @test z.baseline==[0.]
    @test z.ready==[0]
    # First local reward calibrates a baseline, without assigning random credit.
    s.energy[1]=1.;fill!(z.trace,1f0)
    NativeLifetimeCredit.reward!(CPU(),64)(s,status,q,z,.1f0,.3f0,.001,true;ndrange=1)
    @test all(iszero,z.learned)
    @test z.baseline==[1.]
    @test z.ready==[1]
    @test z.paid==[0.]
    @test z.previous==[1.]
end
