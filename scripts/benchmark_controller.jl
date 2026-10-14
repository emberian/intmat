# Bounded controller component benchmark; CPU ecology remains outside this timing.
using Interior, Random
backend_name=isempty(ARGS) ? "metal" : ARGS[1]
backend=if backend_name=="metal"
    @eval using Metal
    Metal.functional() || error("Metal unavailable")
    Metal.MetalBackend()
elseif backend_name=="cuda"
    @eval using CUDA
    CUDA.functional() || error("CUDA unavailable")
    CUDA.CUDABackend()
elseif backend_name=="amd"
    @eval using AMDGPU
    AMDGPU.functional() || error("ROCm unavailable")
    AMDGPU.ROCBackend()
else
    error("select metal, cuda or amd")
end
KA=Interior.KernelAbstractions
rng=MersenneTwister(17)
results=[]
for n in (2400,10000)
    rows,cols=21,44
    weights=0.05f0.*randn(rng,Float32,rows,cols,n)
    marks=0.005f0.*randn(rng,Float32,rows,cols,n)
    inputs=randn(rng,Float32,cols,n); outputs=zeros(Float32,rows,n)
    gpu_weights=KA.zeros(backend,Float32,size(weights)...)
    gpu_marks=KA.zeros(backend,Float32,size(marks)...)
    gpu_inputs=KA.zeros(backend,Float32,size(inputs)...)
    gpu_outputs=KA.zeros(backend,Float32,size(outputs)...)
    copyto!(gpu_weights,weights); copyto!(gpu_marks,marks); copyto!(gpu_inputs,inputs)
    Interior.neural_controller!(outputs,weights,marks,inputs)
    Interior.neural_controller_backend!(gpu_outputs,gpu_weights,gpu_marks,gpu_inputs)
    error_max=maximum(abs.(Array(gpu_outputs).-outputs))
    error_max<1e-5 || error("Controller parity: $error_max")
    steps=20
    cpu_seconds=@elapsed for _ in 1:steps
        Interior.neural_controller!(outputs,weights,marks,inputs)
    end
    resident_seconds=@elapsed for _ in 1:steps
        Interior.neural_controller_backend!(gpu_outputs,gpu_weights,gpu_marks,gpu_inputs)
    end
    host_output=similar(outputs)
    transfer_seconds=@elapsed for _ in 1:steps
        copyto!(gpu_weights,weights); copyto!(gpu_marks,marks); copyto!(gpu_inputs,inputs)
        Interior.neural_controller_backend!(gpu_outputs,gpu_weights,gpu_marks,gpu_inputs)
        copyto!(host_output,gpu_outputs)
    end
    push!(results,(population=n,rows=rows,cols=cols,steps=steps,max_absolute_error=error_max,
        cpu_ms=cpu_seconds/steps*1000,device_resident_ms=resident_seconds/steps*1000,
        upload_all_and_download_ms=transfer_seconds/steps*1000))
end
println(Interior.JSON3.write((backend=backend_name,julia=string(VERSION),measurements=results,
    scope="controller only; excludes sensors, ecology, phase FFT and live state")))
