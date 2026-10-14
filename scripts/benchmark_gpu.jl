using Interior, Metal
const JSON3 = Interior.JSON3
Metal.functional() || error("Metal is not functional on this host")
processes=(WaveForcing(),Excitable(),Diffusion(:nutrient,0.01;material_feedback=2),Deposition())
function trial(backend,n,steps)
    w=World(Grid(n,n,n);backend)
    seed_activity!(w)
    step!(w,processes,0.005)
    seconds=@elapsed for _ in 1:steps
        step!(w,processes,0.005)
    end
    (world=w,seconds=seconds)
end
cpu=trial(CPU(),12,20)
gpu=trial(MetalBackend(),12,20)
errors=Dict(String(k)=>maximum(abs.(Array(getproperty(cpu.world.fields,k)).-Array(getproperty(gpu.world.fields,k)))) for k in keys(cpu.world.fields))
maximum(values(errors)) < 1e-4 || error("CPU / GPU discrepancy: $errors")
cpu_large=trial(CPU(),64,50)
gpu_large=trial(MetalBackend(),64,50)
result=(backend="Metal",comparison_grid=[12,12,12],comparison_steps=21,max_absolute_errors=errors,benchmark_grid=[64,64,64],timed_steps=50,cpu_seconds=cpu_large.seconds,gpu_seconds=gpu_large.seconds,speedup=cpu_large.seconds/gpu_large.seconds)
println(JSON3.write(result))
open(joinpath(@__DIR__,"..","gpu","verification.json"),"w") do io
    JSON3.write(io,result)
end
