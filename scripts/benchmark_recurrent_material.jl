#!/usr/bin/env julia
using Interior
using Interior: JSON3
isdefined(Interior,:RecurrentMaterial) || Base.include(Interior,joinpath(@__DIR__,"..","src","recurrent_material.jl"))
using Interior: RecurrentMaterial
kind=ARGS[1];side=parse(Int,ARGS[2]);steps=parse(Int,ARGS[3]);out=ARGS[4]
kind in ("amd","metal") && side>0 && steps>0 || throw(ArgumentError("valid backend and positive side/steps required"))
ispath(out) && throw(ArgumentError("benchmark output exists; use a fresh destination"))
if kind=="amd"
    @eval using AMDGPU
    backend=AMDGPU.ROCBackend()
else
    @eval using Metal
    backend=Metal.MetalBackend()
end
function runbatch(w,p,n)
    for _ in 1:n;advance!(w,p,.05);end
    Interior.KernelAbstractions.synchronize(w.backend)
end
# Same seeded Julia law on CPU/GPU, including symmetric transport and memory.
a=World(Grid(8,8,8;extent=(8,8,8));seed=31);b=World(a.grid;backend,seed=31)
x=RecurrentMaterial(a);y=RecurrentMaterial(b)
runbatch(a,x,10);runbatch(b,y,10)
error=maximum(abs.(x.state.-Array(y.state)))
error<2e-5 || Base.error("CPU/GPU mismatch $error")
w=World(Grid(side,side,side;extent=Float32.((side,side,side)));backend,seed=31)
p=RecurrentMaterial(w;coupling=.28,feedback=2,writing=.5,erosion=.035,contrast=1.1)
runbatch(w,p,20)
bytes=@allocated runbatch(w,p,20)
elapsed=@elapsed for start in 1:20:steps;runbatch(w,p,min(20,steps-start+1));end
result=(backend=kind,side=side,steps=steps,seconds=elapsed,voxel_steps_per_second=side^3*steps/elapsed,host_bytes_per_step=bytes/20,cpu_gpu_max_error=error,finite=all(isfinite,Array(p.state)),material_mass=material_mass(w))
mkpath(dirname(out));open(out,"w") do io;JSON3.pretty(io,result);end
println(result)
