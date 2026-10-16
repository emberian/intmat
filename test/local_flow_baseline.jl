using Test, Statistics
include("../studies/local_flow_baseline.jl")
const F=LocalFlowBaseline
@testset "flow baseline respects history and partition" begin
    rows=[(id=1,step=0,usable=1.),(id=1,step=10,usable=.4),(id=2,step=10,usable=.6),(id=1,step=20,usable=.3)]
    events=[(parent_id=1,time=.2,child_usable_before_copy=.6)]
    v=F.adjusted(rows,events);f=F.features(rows,v)
    @test v≈[1.,1.,.6,.9]
    @test f[1:3]==zeros(3)
    @test f[4]<0 # growth/use loss remains; only partition is restored
    future=vcat(events,[(parent_id=1,time=1.,child_usable_before_copy=.2)])
    @test F.adjusted(rows,future)==v
    extended=vcat(rows,[(id=1,step=30,usable=.9)])
    @test F.features(extended,F.adjusted(extended,events))[1:4]==f
    @test_throws ErrorException F.features([rows[1],rows[1]],[1.,1.])
end
@testset "pre-action baseline preserves score expectation" begin
    sigma=.3;epsilon=[-sigma,sigma];g=1.2e-4;offset=-.06
    for baseline in (0.,-.05,-.06,.03)
        @test mean((offset.+g.*epsilon.-baseline).*epsilon./sigma^2)≈g atol=1e-15
    end
    # A baseline using this same innovation is invalid and can erase signal.
    @test mean((offset.+g.*epsilon.-g.*epsilon).*epsilon./sigma^2)≈0 atol=1e-15
end
