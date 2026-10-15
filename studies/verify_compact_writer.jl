# Reproducible allocation-layout equivalence check, independent of rendering.
include(joinpath(@__DIR__,"refinement_dense.jl"))
using Test
n=isempty(ARGS) ? 96 : parse(Int,ARGS[1]); steps=length(ARGS)>1 ? parse(Int,ARGS[2]) : 120
@testset "compact writing preserves physical state" begin
    full=World(Grid(n,n,n)); compact=compact_writing_world(full.grid)
    pa=ThermalWriter(full); pb=ThermalWriter(compact)
    for _ in 1:steps; step!(full,pa,.01);step!(compact,pb,.01);end
    @test full.time==compact.time
    for name in (:nutrient,:deposit); @test getproperty(full.fields,name)==getproperty(compact.fields,name);end
    for name in (:temperature,:dose); @test getproperty(pa,name)==getproperty(pb,name);end
    @test isempty(compact.scratch)
    @test propertynames(compact.fields)==(:nutrient,:deposit)
end
