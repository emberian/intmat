using Interior, Test, KernelAbstractions
include(joinpath(@__DIR__,"../studies/native_construction_credit.jl"))
@testset "construction action and forcing boundaries" begin
    N=NativeConstructionCredit
    s=(active=Int32[1],ids=Int32[17],energy=[1.],metabolites=zeros(1,2),inputs=zeros(Float32,1,19))
    s.inputs[1,16:18].=Float32[-.4,.6,.2]
    z=(ids=Int32[0],ready=Int32[0],learned=zeros(Float32,1,4,2),trace=zeros(Float32,1,4,2),drive=zeros(Float32,1,2),requested=zeros(1,2),
       baseline=zeros(1),previous=zeros(1),paid=zeros(1),last_reward=zeros(1),older_reward=zeros(1),critic=zeros(Float32,1,8),features=zeros(Float32,1,8),prediction=zeros(1))
    q=(efficiency=.5,potentials=(4.,4.));status=Int32[0]
    choose(decision,clock)=N.choose!(CPU(),64)(s,status,q,z,0,decision,.3f0,UInt32(71),UInt32(clock),true,true;ndrange=1)
    choose(true,1)
    @test z.ids==s.ids
    @test z.previous==[1.]
    @test z.trace[1,2,:]≈-.4f0.*z.trace[1,1,:]
    @test z.trace[1,3,:]≈.6f0.*z.trace[1,1,:]
    drive=copy(z.drive);trace=copy(z.trace)
    s.inputs[1,16:18].=1f0;choose(false,1)
    @test z.drive==drive
    @test z.trace==trace
    choose(true,2);@test z.trace!=trace
    z.learned.=.1f0;s.ids[1]=18;choose(false,2)
    @test all(iszero,z.learned)
    @test z.ready==[0]
    a=(expression=zeros(Float32,1,3),requests=zeros(Float32,1,2));mass=(active=s.active,mass=Float32[1])
    aq=(secretion_rate=.12f0,allocation_exponent=2f0)
    N.actuate!(CPU(),64)(a,mass,status,aq,z,.02f0,true;ndrange=1)
    @test 0<a.expression[1,1]<1
    @test a.expression[1,2]+a.expression[1,3]≈1
    @test all(a.requests.>0)
    N.actuate!(CPU(),64)(a,mass,status,aq,z,.02f0,false;ndrange=1)
    @test all(iszero,a.requests)
    original=Float32[1 3 0;3 1 0];masks=copy(original)
    N.feed_phase!(CPU(),64)(masks,original,true;ndrange=size(masks))
    @test masks==original[:,[2,1,3]]
    @test sum(masks;dims=1)==sum(original;dims=1)
    N.feed_phase!(CPU(),64)(masks,original,false;ndrange=size(masks))
    @test masks==original
end

struct CountingConstructionLaunch
    calls::Vector{Symbol}
end
function NativeConstructionCredit.L.policy_launch!(launch::CountingConstructionLaunch,h,kernel,args...;ndrange)
    push!(launch.calls,nameof(kernel))
    NativeConstructionCredit.L.policy_launch!(NativeConstructionCredit.L.NativePolicyLaunch(),h,kernel,args...;ndrange)
end
include(joinpath(@__DIR__,"..","studies","constructive_batch_fixture.jl"))
@testset "Native policy launch seam preserves coupled hooks" begin
    w,a=constructive_batch_fixture(71;strategy=:materialized,energy_type=Float64,inventory_type=Float64)
    v,b=constructive_batch_fixture(71;strategy=:materialized,energy_type=Float64,inventory_type=Float64)
    launch=CountingConstructionLaunch(Symbol[])
    native=NativeConstructionCredit.sidecar(a;dt=.1)
    observed=NativeConstructionCredit.sidecar(b;dt=.1,launcher=launch)
    for _ in 1:20
        Interior.device_cellular_step!(w,a.cells,.1;anabolism=native)
        Interior.device_cellular_step!(v,b.cells,.1;anabolism=observed)
    end
    @test constructive_snapshot_equal(Interior.device_habitat_snapshot(a),Interior.device_habitat_snapshot(b))
    @test native.state==observed.state
    @test native.options==observed.options
    @test native.steps==observed.steps==20
    @test count(==(:choose!),launch.calls)==20
    @test count(==(:actuate!),launch.calls)==20
    @test count(==(:feed_phase!),launch.calls)==1
    @test count(==(:reward!),launch.calls)==1
end
