using Test
include(joinpath(@__DIR__,"..","studies","electrical_cdf_diagnostic.jl"))
@testset "archived propensity/CDF diagnostic" begin
 law=electrical_weight_law();occupied=falses(33,33,33);occupied[17,17,17]=true
 p=(potential=zeros(Float32,33,33,33),occupied=occupied,seed=UInt64(7331),events=84,eta=4.8,deposit_per_cell=.25f0)
 sites=[CartesianIndex(18,17,17),CartesianIndex(17,18,17),CartesianIndex(17,17,18)]
 values=Float32[.2,.5,.7];food=ones(Float32,3);dose=zeros(Float32,3);c=Dict("memory_gain"=>1.5)
 a=electrical_cdf_audit(law,p,sites,values,dose,food,c)
 expected=(1 .-Float64.(values)).^4.8
 @test a.weights≈expected
 @test sum(a.probabilities)≈1
 @test a.lower_margin>=0&&a.upper_margin>0
 @test a.cdf[end]≈1
 b=electrical_cdf_audit(law,p,sites,values.+1f-6,dose,food,c)
 comparison=electrical_cdf_compare(a,b)
 @test comparison.total_variation>0
 @test comparison.cdf_maximum_difference>0
 @test comparison.a_choice_certified_against_b
 @test comparison.same_choice
 @test a.uniform===Interior.electrical_uniform(p.seed,85)
 @test occursin("contacts*max",law.generated_source)
end
