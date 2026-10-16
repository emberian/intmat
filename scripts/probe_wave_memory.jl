include(joinpath(@__DIR__,"../studies/wave_memory.jl"))
kind=get(ENV,"WAVE_BACKEND","cpu")
if kind=="metal";using Metal;backend=Metal.MetalBackend()
elseif kind=="amd";using AMDGPU;backend=AMDGPU.ROCBackend()
else;backend=CPU();end
function probe(backend;n=64,steps=100)
 a=memory_fixture(CPU(),16);b=memory_fixture(CPU(),16);c=memory_fixture(backend,16)
 for k=1:100
  wave_memory_step!(a...,.006;scalar=true);wave_memory_step!(b...,.006);wave_memory_step!(c...,.006)
 end
 fa=wave_memory_fields(a...);fb=wave_memory_fields(b...);fc=wave_memory_fields(c...)
 errors=Dict(string(k)=>maximum(abs.(getfield(fa,k).-getfield(fc,k))) for k in keys(fa));cpuerrors=Dict(string(k)=>maximum(abs.(getfield(fa,k).-getfield(fb,k))) for k in keys(fa))
 maximum(values(errors))<1e-4||error("device parity $errors");maximum(values(cpuerrors))==0||error("CPU operation parity $cpuerrors")
 w,p,m=memory_fixture(backend,n);for _=1:10;wave_memory_step!(w,p,m,.006);end
 times=Float64[]
 for _=1:3
  KernelAbstractions.synchronize(backend);elapsed=@elapsed for k=1:steps;wave_memory_step!(w,p,m,.006);end;KernelAbstractions.synchronize(backend);push!(times,elapsed/steps)
 end
 scalarseconds=nothing
 if backend isa CPU
  w,p,m=memory_fixture(backend,n);for _=1:5;wave_memory_step!(w,p,m,.006;scalar=true);end
  scalarseconds=(@elapsed for _=1:steps;wave_memory_step!(w,p,m,.006;scalar=true);end)/steps
 end
 result=(backend=kind,n=n,steps_per_repeat=steps,repeats=3,complete_step_seconds=times,scalar_complete_step_seconds=scalarseconds,cpu_exact_errors=cpuerrors,cpu_device_errors=errors,persistent_state_arrays=10,field_bytes=10*4*n^3,qualification="Full acoustic validation/reduction, production wave kernel, pressure mirror copy, material/dose kernel and clock update; warmed simulation excludes initialization, export and rendering; no asynchronous unsynchronized timing.")
 open(io->JSON3.pretty(io,result),ARGS[1],"w");println(JSON3.write(result))
end
probe(backend;n=parse(Int,get(ENV,"WAVE_N","64")),steps=parse(Int,get(ENV,"WAVE_STEPS","40")))
