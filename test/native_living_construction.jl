using Interior, Test, KernelAbstractions
include(joinpath(@__DIR__,"../studies/native_living_construction.jl"))
@testset "birth credit transfers and paid private inheritance" begin
    M=NativeLivingConstruction
    s=(active=Int32[1,1],ids=Int32[1,2],new_parent=Int32[0,1],energy=[1.,1.],metabolites=[.1 .15;.1 .15],ledger=zeros(Float32,2,9),ledger_low=zeros(Float32,2,9))
    z=(ids=Int32[1,17],ready=Int32[1,1],baseline=zeros(2),last_reward=zeros(2),older_reward=zeros(2),paid=zeros(2),birth_pending=zeros(Int32,2),invalidated=zeros(Int32,2),
        queue_ids=Int32[1 999;17 0],queue_value=[5. 6.;7. 8.],critic=ones(Float32,2,8),features=ones(Float32,2,8),current_features=ones(Float32,2,8),
        learned=fill(.1f0,2,2,3),trace=ones(Float32,2,2,3))
    ev=(child_id=zeros(Int32,8),parent_id=zeros(Int32,8),time=zeros(8),gift=zeros(8),copied_l1=zeros(8),copy_work=zeros(8),funded=zeros(Int32,8),pending=zeros(Int32,8))
    q=(efficiency=.65,potentials=(4.,4.),kappa=.001);status=Int32[0]
    gift=M.L.local_value(s,q,2)
    before=(energy=[2.],metabolites=reshape([.2,.3],1,2))
    @test M.L.local_value(before,q,1)==M.L.local_value(s,q,1)+gift
    a,b=deepcopy(z),deepcopy(z);sa,sb=deepcopy(s),deepcopy(s);ea,eb=deepcopy(ev),deepcopy(ev)
    M.birth!(CPU(),64)(sa,status,q,a,ea,1.,0f0,true;ndrange=2)
    M.birth!(CPU(),64)(sb,status,q,b,eb,1.,0f0,false;ndrange=2)
    @test sa==sb==s # adjusting credit has no physical inventory or ledger write
    @test a.queue_value[1,1]==5-gift
    @test a.queue_value[1,2]==6 # foreign identity cannot receive credit
    @test b.queue_value[1,:]==[5.,6.]
    @test ea.pending[2]==1
    @test all(iszero,a.queue_ids[2,:])
    @test a.invalidated[2]==1
    @test all(iszero,a.critic[2,:])
    @test a.ids[2]==2
    original=copy(a.queue_value)
    M.birth!(CPU(),64)(sa,status,q,a,ea,1.,0f0,true;ndrange=2)
    @test a.queue_value==original # one event, not one adjustment per hook invocation
    a=deepcopy(z);sa=deepcopy(s);ea=deepcopy(ev)
    M.birth!(CPU(),64)(sa,status,q,a,ea,1.,.5f0,true;ndrange=2)
    @test a.learned[2,:,:]==.5f0.*z.learned[1,:,:]
    @test ea.gift[2]==gift # specified before optional child-copy debit
    @test ea.copy_work[2]>=q.kappa*sum(abs,Float64.(a.learned[2,:,:]))
    @test s.energy[2]-sa.energy[2]==ea.copy_work[2]
    @test Float64(sa.ledger[2,1])+sa.ledger_low[2,1]≈ea.copy_work[2] atol=1e-18
    a=deepcopy(z);sa=deepcopy(s);ea=deepcopy(ev);sa.energy[2]=0.
    M.birth!(CPU(),64)(sa,status,q,a,ea,1.,.5f0,true;ndrange=2)
    @test all(iszero,a.learned[2,:,:])
    @test ea.funded[2]==0
    @test ea.copy_work[2]==0
    @test sa.energy[2]==0
end
