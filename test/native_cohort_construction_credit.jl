using Interior, Test, KernelAbstractions
include(joinpath(@__DIR__,"../studies/native_cohort_construction_credit.jl"))
@testset "decision-aligned cohort timing and paid transaction" begin
    D=NativeCohortConstructionCredit;C=D.C;L=D.L
    f=C.fixture(joinpath(@__DIR__,"../experiments/native/constructed-garden-audited-64.toml");reserve_fraction=.2)
    w=World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed);p=Interior.device_habitat_restore(w,f.checkpoint)
    h=D.sidecar(p;horizon=.4,act_time=1.2);origin=C.A.audit(p)
    L.advance!(w,p,h,.2)
    @test all(h.state.enqueued[1:4].==1)
    @test all(iszero,h.state.matured)
    prediction=copy(h.state.queue_prediction[:,1]);features=copy(h.state.queue_features[:,:,1])
    L.advance!(w,p,h,.2)
    @test all(h.state.matured[1:4].==1)
    @test all(h.state.skipped[1:4].==1)
    @test all(iszero,h.state.learned)
    @test h.state.prediction[1:4]==prediction[1:4]
    @test h.state.queue_prediction[:,1]==prediction
    @test h.state.queue_features[:,:,1]==features
    L.advance!(w,p,h,.8)
    @test all(h.state.enqueued[1:4].==6)
    @test all(h.state.matured[1:4].==5)
    @test all(h.state.skipped[1:4].==2)
    @test any(!iszero,h.state.learned)
    L.advance!(w,p,h,.4)
    @test all(h.state.enqueued[1:4].==6)
    @test all(h.state.matured[1:4].==6)
    @test all(iszero,h.state.queue_ids)
    @test all(iszero,h.state.invalidated)
    observed=C.observation(w,p,h,origin)
    @test observed.all_residents_valid
    @test abs(observed.potential)<1e-8
    @test abs(observed.material)<1e-8
    @test all(h.state.paid[1:4].>0)
    @test all(h.state.funded_increment.<=h.state.proposed_increment)
    # Explicit default gain must be bitwise identical to historical implicit .3.
    chemistry=(efficiency=p.cells.parameters.efficiency,potentials=p.cells.parameters.potentials)
    s=p.cells.state;v=(active=s.active,ids=s.ids,energy=s.energy,metabolites=s.metabolites,ledger=s.ledger,ledger_low=s.ledger_low)
    v1=deepcopy(v);v2=deepcopy(v);z1=deepcopy(h.state);z2=deepcopy(h.state)
    fill!(z1.features,.1f0);fill!(z2.features,.1f0)
    L.reward!(CPU(),64)(v1,Int32[0],chemistry,z1,.003f0,.3f0,.001,true,true,true;ndrange=256)
    L.reward!(CPU(),64)(v2,Int32[0],merge(chemistry,(critic_gain=.3f0,)),z2,.003f0,.3f0,.001,true,true,true;ndrange=256)
    @test v1==v2
    @test z1==z2
    # Actor keeps its old baseline; predictor fits current coefficients on φ.
    va=deepcopy(v);za=deepcopy(h.state);fill!(za.learned,0f0);fill!(za.critic,0f0)
    fill!(za.trace,1f0);fill!(za.ready,1);fill!(za.features,0f0);za.features[:,1].=1f0
    fill!(za.prediction,.5)
    for i in 1:4;za.previous[i]=L.local_value(va,chemistry,i)-.1;end
    L.reward!(CPU(),64)(va,Int32[0],merge(chemistry,(critic_gain=.03f0,critic_current=true)),za,.003f0,.3f0,.001,true,true,true;ndrange=256)
    @test all(za.learned[1:4,:,:].<0)
    @test all(za.critic[1:4,1].>0)
    @test all(isapprox.(za.prediction_error[1:4],-.4;atol=1e-12))
    h.state.queue_ids[1,1]=999
    before=h.state.invalidated[1]
    D.mature!(CPU(),64)((active=s.active,ids=s.ids,energy=s.energy,metabolites=s.metabolites),Int32[0],merge(chemistry,(critic_gain=.03f0,)),h.state,1,.003f0;ndrange=256)
    @test h.state.invalidated[1]==before+1
    @test h.state.transaction[1]==0
    @test h.state.queue_ids[1,1]==0
end
