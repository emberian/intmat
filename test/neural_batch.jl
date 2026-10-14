using Test, Random, KernelAbstractions
# Also runnable before module integration.
isdefined(Interior,:neural_controller!) || Base.include(Interior,joinpath(@__DIR__,"../src/neural_batch.jl"))
@testset "batched recurrent controller reference" begin
    rng=MersenneTwister(512)
    for n in (0,1,37), h in (4,12)
        weights=randn(rng,Float32,h+9,h+32,n)
        marks=.1f0 .* randn(rng,Float32,size(weights))
        input=randn(rng,Float32,h+32,n)
        reference=hcat([tanh.((weights[:,:,i]+marks[:,:,i])*input[:,i]) for i in 1:n]...)
        out=zeros(Float32,h+9,n)
        Interior.neural_controller!(out,weights,marks,input)
        n>0 && @test out ≈ reference atol=2f-6 rtol=2f-5
        @test all(isfinite,out)
        Interior.neural_controller_backend!(out,weights,marks,input)
        n>0 && @test out ≈ reference atol=2f-6 rtol=2f-5
        @test all(isfinite,out)
    end
    @test_throws DimensionMismatch Interior.neural_controller!(zeros(Float32,2,1),zeros(Float32,2,3,1),zeros(Float32,2,3,1),zeros(Float32,4,1))
end
@testset "ecology controller drop-in" begin
    rng=MersenneTwister(513)
    weights=randn(rng,Float32,21,44,23); marks=.1f0 .*randn(rng,Float32,size(weights)); inputs=[randn(rng,Float32,44) for _ in 1:23]
    reference=[tanh.((weights[:,:,i]+marks[:,:,i])*inputs[i]) for i in 1:23]
    outputs=Interior.neural_controller_outputs(weights,marks,inputs)
    @test all(isapprox(x,y;atol=2f-6,rtol=2f-5) for (x,y) in zip(outputs,reference))
    @test all(x->all(isfinite,x),outputs)
end
