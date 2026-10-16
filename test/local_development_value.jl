using Test,Interior
include(joinpath(@__DIR__,"../studies/local_development_value.jl"))
const DV=LocalDevelopmentValue
@testset "stoichiometric developmental accounting" begin
    q=(efficiency=.65,potentials=(4.,4.),growth_work=.4)
    @test DV.substrate_coefficient(q)==2.6
    @test DV.replacement_coefficient(q)==3.
    for T in (Float32,Float64),a0 in (.01,.2),b0 in (.03,.3),e0 in (.001,.4)
        a,b,m,e=T(a0),T(b0),T(1.),T(e0)
        before=DV.replacement_value(q,m,e,a,b);cap=DV.buildable_capacity(q,m,e,a,b)
        aa,bb,mm,ee,work,heat,g=Interior.device_anabolic_build!(nothing,1,a,b,m,e,T(.05),T(.4),T(4),T(4))
        @test isapprox(DV.replacement_value(q,mm,ee,aa,bb),before;atol=5e-7)
        @test isapprox(DV.buildable_capacity(q,mm,ee,aa,bb),cap;atol=5e-7)
        @test isapprox(mm+aa+bb,m+a+b;atol=5e-7)
        @test isapprox(ee+4aa+4bb+work+heat,e+4a+4b;atol=5e-7)
        @test DV.replacement_value(q,m/2,e/2,a/2,b/2)*2==before
        @test DV.buildable_capacity(q,m/2,e/2,a/2,b/2)*2==cap
    end
    @test DV.replacement_value(q,0.,0.,0.,0.)==0
end
