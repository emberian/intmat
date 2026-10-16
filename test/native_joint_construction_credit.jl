using Interior, Test, KernelAbstractions
include(joinpath(@__DIR__,"../studies/native_joint_construction_credit.jl"))
@testset "joint construction feature and reserve boundaries" begin
    J=NativeJointConstructionCredit
    s=(active=Int32[1],ids=Int32[17],energy=[1.],metabolites=reshape([.2,.4],1,2),mass=[2.],
       inputs=zeros(Float32,1,22),hidden=reshape(Float32[.3],1,1),adaptation=reshape(Float32[-.2],1,1),actions=zeros(Float32,1,11))
    s.inputs[1,18:20].=Float32[-.4,.6,.2]
    q=(efficiency=.5,potentials=(4.,4.),energy_capacity=2f0,store_capacity=1f0)
    @test [J.feature(s,q,1,1,c,true) for c in 1:9]≈Float32[1,-.4,.6,.2,.25,.1,.2,.3,-.2]
    @test [J.feature(s,q,1,1,c,false) for c in 8:9]==zeros(Float32,2)
    z=(ids=Int32[0],ready=Int32[0],learned=zeros(Float32,1,9,3),trace=zeros(Float32,1,9,3),drive=zeros(Float32,1,3),
       baseline=zeros(1),previous=zeros(1),paid=zeros(1),last_reward=zeros(1),older_reward=zeros(1),critic=zeros(Float32,1,8),features=zeros(Float32,1,8),prediction=zeros(1))
    status=Int32[0]
    J.choose_joint!(CPU(),64)(s,status,q,z,1,true,0f0,UInt32(71),UInt32(1),true,true,true;ndrange=1)
    J.reserve!(CPU(),64)(s,status,z;ndrange=1)
    @test s.actions[1,11]≈-.6f0
    @test (s.actions[1,11]+1)*.5≈.2 atol=3e-8
    @test all(iszero,z.trace)
    J.choose_joint!(CPU(),64)(s,status,q,z,1,true,.3f0,UInt32(71),UInt32(2),false,true,false;ndrange=1)
    @test all(iszero,z.trace[:,8:9,:])
    @test any(!iszero,z.trace[:,1:7,:])
    drive=copy(z.drive);trace=copy(z.trace);s.hidden.=.9f0;s.energy[1]=.5
    J.choose_joint!(CPU(),64)(s,status,q,z,1,false,.3f0,UInt32(71),UInt32(2),false,true,true;ndrange=1)
    @test z.drive==drive
    @test z.trace==trace
    z.learned.=.1f0;s.ids[1]=18
    J.choose_joint!(CPU(),64)(s,status,q,z,1,false,0f0,UInt32(71),UInt32(2),false,true,true;ndrange=1)
    @test all(iszero,z.learned)
    @test z.ready==[0]
end
