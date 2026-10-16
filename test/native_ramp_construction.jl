using Test
include("../studies/native_ramp_construction.jl")
const R=NativeRampConstruction
@testset "constant gain, observable current, and history sign" begin
    for balance in (-1f0,-.02f0,0f0,.02f0,1f0),trend in (-.1f0,0f0,.1f0)
        h=R.actor_features(balance,trend,false,false)
        r=R.actor_features(balance,trend,true,false)
        s=R.actor_features(balance,trend,false,true)
        @test isapprox(sum(abs2,h),2f0;atol=2f-7)
        @test sum(abs2,h)==sum(abs2,r)==sum(abs2,s)
        @test h[1:2]==r[1:2]==s[1:2]
        @test h[3]==-s[3]
        @test R.actor_features(balance,trend,true,true)==r
    end
    @test R.actor_features(0f0,0f0,false,false)==(1f0,0f0,1f0)
    @test R.actor_features(0f0,-.1f0,false,false)==(1f0,0f0,-1f0)
end
