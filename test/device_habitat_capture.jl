using Interior,Test
include(joinpath(@__DIR__,"..","scripts","device_habitat_capture.jl"))
@testset "Device habitat block mean field units" begin
    g=(dims=(Int32(8),Int32(8),Int32(8)),volume=.125f0)
    food=fill(.5f0*g.volume,512,3)
    # A varying field verifies averaging rather than nearest-point sampling.
    food[:,1].=Float32.(1:512).*g.volume
    e=(fields=food,permeability=fill(.25f0,512))
    h=(enzyme_bits=fill(reinterpret(UInt32,g.volume),512,2),bound=fill(3f0*g.volume,512,2))
    r=(fields=fill(2f0*g.volume,512,3),)
    dst=zeros(Float32,64,9)
    habitat_coarse!(CPU(),128)(dst,e,h,r,g,Int32(4);ndrange=64)
    @test all(dst[:,2:9].==permutedims(Float32[.5,2,2,1,1,3,3,.25]))
    full=reshape(Float32.(1:512),8,8,8)
    expected=[sum(full[2x-1:2x,2y-1:2y,2z-1:2z])/8 for x in 1:4,y in 1:4,z in 1:4]
    @test reshape(dst[:,1],4,4,4)==expected
end
