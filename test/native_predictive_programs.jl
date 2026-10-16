using Test
include("../studies/native_predictive_programs.jl")
const P=NativePredictivePrograms
@testset "predictive programs" begin
    @test P.program_action(0)==(0f0,.5f0)
    @test P.program_action(1)==(.3f0,.9f0)
    @test P.concentration(P.receptor(.25f0))≈.25f0
    f=P.N.fixture(joinpath(@__DIR__,"../experiments/native/living-temporal-construction-16.toml"));w,p=f.w,f.p
    h=P.sidecar(p;forced=1);q=p.cells.parameters;hq=p.habitat.parameters;o=h.program_options
    @test isfinite(P.program_score(1,.1f0,.1f0,1f0,1f0,.2f0,.2f0,.5f0,.5f0,q,merge(hq,(decay_rate=0f0,)),o))
    @test P.program_score(1,.1f0,.1f0,1f0,1f0,0f0,0f0,.5f0,.5f0,q,hq,o)<0
    @test_throws ErrorException P.sidecar(p;capture=NaN)
    @test_throws ErrorException P.sidecar(p;dwell=.3)
    origin=P.N.C.A.audit(p)
    P.L.advance!(w,p,h,.42;dt=.02)
    obs=P.observation(w,p,h)
    @test all(c.phase==2 && c.investment==0 for c in obs.cells)
    @test sum(h.model.funded)>0
    @test sum(h.model.work)>0
    @test sum(h.model.updates)>0
    @test all(c.error_target_step==20 for c in obs.cells)
    @test all(c.prediction_for_step==30 for c in obs.cells)
    a=deepcopy((w=w,p=p,h=h));b=deepcopy((w=w,p=p,h=h))
    P.L.advance!(a.w,a.p,a.h,.4;dt=.02)
    P.L.advance!(b.w,b.p,b.h,.2;dt=.02);P.L.advance!(b.w,b.p,b.h,.2;dt=.02)
    @test P.observation(a.w,a.p,a.h)==P.observation(b.w,b.p,b.h)
    @test a.h.model.weights==b.h.model.weights
    # Unaffordable predictor write is atomic; identity replacement erases old evidence.
    c=deepcopy(a);oldweights=copy(c.h.model.weights);oldwork=copy(c.h.model.work)
    fill!(c.p.cells.state.energy,0.);c.h.model.obs[:,1].+=.01f0
    step=Int(c.h.model.target[1])
    P.predict_choose!(c.h.backend,64)(c.p.cells.state,c.p.cells.environment.status,c.p.cells.parameters,c.p.habitat.parameters,c.h.state,c.h.model,c.h.program_options,step;ndrange=Int(c.p.cells.geometry.capacity))
    @test c.h.model.weights==oldweights
    @test c.h.model.work==oldwork
    c.p.cells.state.ids[1]=9999
    P.predict_choose!(c.h.backend,64)(c.p.cells.state,c.p.cells.environment.status,c.p.cells.parameters,c.p.habitat.parameters,c.h.state,c.h.model,c.h.program_options,step+1;ndrange=Int(c.p.cells.geometry.capacity))
    @test all(iszero,c.h.model.weights[1,:,:])
    @test c.h.model.target[1]==(step÷c.h.program_options.stride+1)*c.h.program_options.stride
    frame=P.N.frame(a.w,a.p,a.h,origin)
    @test abs(frame.ledger.material_residual)<1e-8
    @test abs(frame.ledger.potential_residual)<1e-8
end
