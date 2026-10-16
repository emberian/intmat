include(joinpath(@__DIR__,"device_cellular_setup.jl"))
using AMDGPU, LinearAlgebra
const FFTW=Interior.FFTW
const fftlib=AMDGPU.librocfft
@assert AMDGPU.librocfft!=""
rows=[]
for shape in ((9,7,5),(16,16,16),(96,96,96))
 host=reshape(sin.(Float64.(1:prod(shape))),shape)
 x=ROCArray(host)
 setup=@elapsed begin
  f=FFTW.plan_rfft(x);y=f*x;r=FFTW.plan_brfft(y,shape[1]);z=r*y;z .*= inv(Float64(prod(shape)));AMDGPU.synchronize()
 end
 forward_error=maximum(abs,Array(y).-FFTW.rfft(host))
 roundtrip_error=maximum(abs,Array(z).-host)
 forward_error<1e-8&&roundtrip_error<1e-12||error("F64 native FFT parity failed")
 t=[@elapsed begin mul!(y,f,x);mul!(z,r,y);z .*= inv(Float64(prod(shape)));AMDGPU.synchronize() end for _ in 1:3]
 push!(rows,(;shape,setup_seconds=setup,forward_error,roundtrip_error,roundtrip_seconds=t,element_type=string(eltype(x))))
 println(JSON3.write(last(rows)));flush(stdout)
end
open(io->JSON3.pretty(io,(;version=string(AMDGPU.rocFFT.version()),library=fftlib,rows)),ARGS[1],"w")
