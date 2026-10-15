using Test
include(joinpath(@__DIR__,"..","studies","electrical_patch_kernel.jl"))
@testset "single-workgroup electrical patch" begin
 for n in (17,33),center in ((2,2,2),(n-1,n-1,n-1),ntuple(_->cld(n,2),3)),sweeps in (0,1,16)
  p=process(n,specimen(n));mgsolve!(hierarchy(p.potential,p.occupied);tol=1.6e-5)
  I=CartesianIndex(center);Interior.electrical_mark!(p,I)
  expected=copy(p.potential);candidate=copy(expected)
  electrical_patch!(expected,p.occupied,I;radius=8,sweeps)
  stats=electrical_device_patch8!(candidate,UInt8.(p.occupied),I;sweeps)
  KernelAbstractions.synchronize(KernelAbstractions.get_backend(candidate))
  @test reinterpret(UInt32,vec(candidate))==reinterpret(UInt32,vec(expected))
  @test abs(stats.residual-Float64(electrical_patch_residual(expected,p.occupied,I,8)))<=eps(Float32)*max(1.,stats.residual)
 end
 @test_throws ArgumentError electrical_device_patch8!(zeros(Float32,17,17,17),zeros(UInt8,17,17,17),CartesianIndex(1,2,2))
 @test_throws ArgumentError electrical_device_patch8!(zeros(Float32,17,17,17),zeros(UInt8,17,17,17),CartesianIndex(2,2,2);sweeps=-1)
end
