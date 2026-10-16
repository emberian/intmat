using Test
include(joinpath(@__DIR__,"../studies/native_development_living.jl"))
const W=NativeDevelopmentLiving
const V=W.V
@testset "developmental valuation and conservative birth" begin
    q=(efficiency=.65,potentials=(4.,4.));s=(energy=[1.],metabolites=reshape([.2,.3],1,2))
    @test V.L.local_value(s,q,1)==V.L.local_value(s,merge(q,(development_coefficient=0.,)),1)
    @test V.L.local_value(merge(s,(mass=[NaN],)),merge(q,(development_coefficient=0.,)),1)==V.L.local_value(s,q,1)
    @test V.L.local_value(merge(s,(mass=[2.],)),merge(q,(development_coefficient=3.,)),1)==V.L.local_value(s,q,1)+6
    physical=(active=Int32[1,1],ids=Int32[7,8],new_parent=Int32[0,1],mass=[.9,.9])
    z=(ids=Int32[7,8],queue_ids=Int32[7 99;0 0],queue_value=[10. 11.;0. 0.])
    events=(child_id=Int32[0,0,0,0,0,0,0,8],)
    cap=(mass=zeros(8),recorded=zeros(Int32,8));status=Int32[0]
    W.capital_birth!(W.KernelAbstractions.CPU(),64)(physical,status,z,events,cap,3.,true;ndrange=2)
    @test z.queue_value[1,1]≈7.3
    @test z.queue_value[1,2]==11.
    @test cap.mass[8]==.9 && cap.recorded[8]==1
    saved=copy(z.queue_value)
    W.capital_birth!(W.KernelAbstractions.CPU(),64)(physical,status,z,events,cap,3.,true;ndrange=2)
    @test z.queue_value==saved
end
@testset "zero valuation exact law and restart" begin
    profile=joinpath(@__DIR__,"../experiments/native/living-temporal-construction-16.toml")
    cfg=V.TOML.parsefile(profile);cfg["cellular"]["biomass"]=1.79
    path=tempname()*".toml";open(io->V.TOML.print(io,cfg),path,"w")
    f=V.N.fixture(path);g=deepcopy(f);d=deepcopy(f)
    kwargs=(rate=.1,decision_time=1.,observation_time=.2,horizon=4.,act_time=12.,reserve=.02)
    h=V.sidecar(f.p;kwargs...);j=W.sidecar(g.p;coefficient=0.,kwargs...);k=W.sidecar(d.p;kwargs...)
    V.L.advance!(f.w,f.p,h,9.;dt=.02);V.L.advance!(g.w,g.p,j,9.;dt=.02)
    @test isequal(V.Interior.device_habitat_snapshot(f.p),V.Interior.device_habitat_snapshot(g.p))
    @test h.state==j.state
    @test h.temporal==j.temporal
    V.L.advance!(d.w,d.p,k,.42;dt=.02)
    fork=deepcopy((w=d.w,p=d.p,h=k))
    V.L.advance!(d.w,d.p,k,8.58;dt=.02);V.L.advance!(fork.w,fork.p,fork.h,8.58;dt=.02)
    @test isequal(V.Interior.device_habitat_snapshot(d.p),V.Interior.device_habitat_snapshot(fork.p))
    @test k.state==fork.h.state && k.capital==fork.h.capital
    @test sum(k.capital.recorded)>0
    @test sum(k.state.actor_updates)>0 && sum(k.state.paid)>0
    @test k.coefficient==Float64(W.Value.replacement_coefficient(d.p.cells.parameters))
    # Original maturity consumes every valid queue even if the paid update is
    # unaffordable, or calibration suppresses the actor; invalid IDs are dropped.
    z=k.state;s=d.p.cells.state;i=findfirst(==(1),s.active);slot=1
    q=(efficiency=d.p.cells.parameters.efficiency,potentials=d.p.cells.parameters.potentials,
       critic_gain=k.options.critic_gain,critic_current=true,development_coefficient=k.coefficient)
    for ready in (0,1)
        s.energy[i]=0.;z.queue_ids[i,slot]=s.ids[i];z.queue_ready[i,slot]=ready;z.queue_value[i,slot]=0.
        old=copy(z.learned);paid=z.paid[i]
        W.V.D.mature!(W.KernelAbstractions.CPU(),64)(s,d.p.cells.environment.status,q,z,slot,.1f0;ndrange=length(s.active))
        view=(active=z.transaction,ids=s.ids,energy=s.energy,metabolites=s.metabolites,mass=s.mass,ledger=s.ledger,ledger_low=s.ledger_low)
        V.L.reward!(W.KernelAbstractions.CPU(),64)(view,d.p.cells.environment.status,q,z,.1f0,.3f0,.001,true,true,true;ndrange=length(s.active))
        @test z.queue_ids[i,slot]==0
        @test z.learned[i,:,:]==old[i,:,:] && z.paid[i]==paid && s.energy[i]==0
    end
    z.queue_ids[i,slot]=s.ids[i]+100000;previous=z.previous[i]
    W.V.D.mature!(W.KernelAbstractions.CPU(),64)(s,d.p.cells.environment.status,q,z,slot,.1f0;ndrange=length(s.active))
    @test z.queue_ids[i,slot]==0 && z.transaction[i]==0 && z.previous[i]==previous
end
