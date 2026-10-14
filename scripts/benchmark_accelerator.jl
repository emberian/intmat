#= Bounded dense-field benchmark; does not benchmark or change the live arena.
Usage: julia --project=ENV scripts/benchmark_accelerator.jl cpu|metal|cuda|amd [N=32] [STEPS=10]
The selected optional GPU package must already be installed in ENV.
=#
using Interior

backend_name = isempty(ARGS) ? "cpu" : lowercase(ARGS[1])
n = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 32
steps = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 10
4 <= n <= 128 || error("N must be between 4 and 128 for this bounded probe")
1 <= steps <= 100 || error("STEPS must be between 1 and 100")
backend = if backend_name == "cpu"
    CPU()
elseif backend_name == "metal"
    @eval using Metal
    Metal.functional() || error("Metal runtime unavailable")
    Metal.MetalBackend()
elseif backend_name == "cuda"
    @eval using CUDA
    CUDA.functional() || error("CUDA runtime unavailable")
    CUDA.CUDABackend()
elseif backend_name == "amd"
    @eval using AMDGPU
    AMDGPU.functional() || error("AMDGPU/ROCm runtime unavailable")
    AMDGPU.ROCBackend()
else
    error("backend must be cpu, metal, cuda or amd")
end
processes = (WaveForcing(), Excitable(), Diffusion(:nutrient, 0.01; material_feedback=2), Deposition())
# Match the default dense workload, but respect the most restrictive initial bound.
function prepare(b, size)
    w = World(Grid(size,size,size); backend=b)
    seed_activity!(w)
    w
end
reference = prepare(CPU(), 8)
candidate = prepare(backend, 8)
dt = min(0.005f0, 0.5f0 * minimum(stability_limit(reference,p) for p in processes))
for _ in 1:5
    step!(reference,processes,dt)
    step!(candidate,processes,dt)
end
errors = Dict(String(k) => maximum(abs.(Array(getproperty(reference.fields,k)) .- Array(getproperty(candidate.fields,k)))) for k in keys(reference.fields))
all(isfinite, values(errors)) && maximum(values(errors)) < 1e-4 || error("CPU/device parity failed: $errors")
w = prepare(backend,n)
dt = min(0.005f0, 0.5f0 * minimum(stability_limit(w,p) for p in processes))
step!(w,processes,dt) # Compile and warm the exact shape before timing.
seconds = @elapsed for _ in 1:steps
    step!(w,processes,dt)
end
println(Interior.JSON3.write((backend=backend_name, julia=string(VERSION), grid=[n,n,n],
    timed_steps=steps, dt=dt, seconds=seconds, cells_per_second=n^3*steps/seconds,
    parity_max_absolute_errors=errors, workload="dense fields; excludes continuous arena ecology and implicit FFT phase")))
