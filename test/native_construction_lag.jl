using Interior, Test, KernelAbstractions
include(joinpath(@__DIR__,"../studies/native_construction_lag.jl"))
@testset "causal actuator pulse restores exact sampled drive" begin
    P=NativeConstructionLag
    s=(active=Int32[1,1],ids=Int32[17,18]);status=Int32[0]
    z=(drive=Float32[.1 .2 -.6931472; .3 .4 .5],);saved=zero(z.drive);before=copy(z.drive)
    for axis in Int32.(1:3),amp in Float32[-.15,0,.15]
        P.perturb_drive!(CPU(),64)(s,status,z,saved,axis,amp,Int32(17),false;ndrange=2)
        @test z.drive[1,axis]===before[1,axis]+amp
        @test z.drive[2,:]==before[2,:]
        P.perturb_drive!(CPU(),64)(s,status,z,saved,axis,amp,Int32(17),true;ndrange=2)
        @test z.drive==before
    end
    status[1]=8
    P.perturb_drive!(CPU(),64)(s,status,z,saved,Int32(1),.15f0,Int32(17),false;ndrange=2)
    @test z.drive==before
end
