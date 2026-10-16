using Test
include("../studies/homeostatic_temporal_memory.jl")
const M=HomeostaticTemporalMemory
@testset "input homeostasis preserves conditional memory poles" begin
    q=M.options()
    a=M.State();b=deepcopy(a);a.fast=.3;a.slow=-.4;b.fast=-.2;b.slow=.1
    df=a.fast-b.fast;ds=a.slow-b.slow
    for u in (.2,.7,-.8,0.,.6)
        M.observe!(a,u,q);M.observe!(b,u,q)
        df*=1-q.af;ds*=1-q.as
        @test a.fast-b.fast≈df atol=1e-14
        @test a.slow-b.slow≈ds atol=1e-14
        @test (a.center,a.variance,a.battery)==(b.center,b.variance,b.battery)
    end
    s=M.State(battery=0.);old=(s.center,s.variance)
    r=M.observe!(s,.8,q)
    @test (s.center,s.variance)==old
    @test r.paid==s.work==0.
    @test s.fast!=0. # losing update funding does not erase retained state
    for _ in 1:10000;M.observe!(a,.6,q);end
    @test isfinite(M.gain(a,q)) && M.gain(a,q)<=100
    @test abs(a.fast)<=1 && abs(a.slow)<=1
end
@testset "bounded changed-scale prediction and causal history intervention" begin
    a=M.replay();f=M.replay(adapt=false,yoked=a.payments);d=M.replay(battery=0.)
    @test a.work≈f.work atol=1e-15
    @test d.prediction==f.prediction
    @test a.metrics[2].late_mse<f.metrics[2].late_mse/5
    @test a.metrics[1].late_mse<.01 && f.metrics[1].late_mse<.01
    # Same instantaneous raw value, opposite preceding history/future phase.
    i=1501;j=1521 # t300 and t304
    @test a.inputs[i]≈a.inputs[j] atol=1e-14
    @test a.prediction[i]<-.4 && a.prediction[j]>.4
    @test f.prediction[i]<0 && f.prediction[j]<0
    # Frozen readout/physical input: replace only fast/slow state by the
    # opposite-history values. Feature norms are not matched; this tests
    # dependence, while scale robustness is the separately measured capability.
    swapped=a.weights[1]+a.weights[2]*a.features[j,2]+a.weights[3]*a.features[j,3]
    @test swapped≈a.prediction[j]
    @test a.weights[1]≈0 atol=.001 # erasing both memories removes phase forecast
end
