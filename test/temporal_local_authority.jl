using Test
include(joinpath(@__DIR__,"../studies/temporal_local_authority.jl"))
const A=TemporalLocalAuthority
@testset "single resident pulse and exact restoration" begin
    s=(active=Int32[1,1,0],ids=Int32[7,8,9]);status=Int32[0]
    z=(drive=reshape(Float32.(1:9),3,3),);saved=zeros(Float32,3);original=copy(z.drive)
    A.pulse!(A.KernelAbstractions.CPU(),64)(s,status,z,saved,Int32(8),.3f0,false;ndrange=3)
    @test z.drive[2,2]==original[2,2]+.3f0
    @test z.drive[[1,3],:]==original[[1,3],:]
    A.pulse!(A.KernelAbstractions.CPU(),64)(s,status,z,saved,Int32(8),.3f0,true;ndrange=3)
    @test isequal(z.drive,original)
    A.pulse!(A.KernelAbstractions.CPU(),64)(s,status,z,saved,Int32(0),-.3f0,false;ndrange=3)
    @test z.drive[1:2,2]==original[1:2,2].-.3f0
    @test z.drive[3,:]==original[3,:]
    A.pulse!(A.KernelAbstractions.CPU(),64)(s,status,z,saved,Int32(0),-.3f0,true;ndrange=3)
    @test isequal(z.drive,original)
end
