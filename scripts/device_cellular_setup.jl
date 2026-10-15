# Shared standalone device experiment setup; no package/runtime mutation.
using Interior, Random
using Interior: JSON3
isdefined(Interior,:DeviceCellular) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_cellular.jl"))
# Refuse to archive a different repository revision than the loaded process.
for (role,file) in ((:model,"device_cellular.jl"),(:support,"device_cellular_support.jl"),(:initialization,"device_cellular_initialization.jl"),(:founder_initialization,"cellular_initialization.jl"),(:resource,"device_resource_exchange.jl"),(:kernels,"device_cellular_kernels.jl"),(:checkpoint,"device_cellular_checkpoint.jl"),(:protocol,"device_cellular_protocol.jl"))
    path=joinpath(@__DIR__,"..","src",file)
    actual=isfile(path) ? bytes2hex(Interior.SHA.sha256(read(path))) : nothing
    hasproperty(Interior.DEVICE_CELLULAR_SOURCE_PIN,role) && getproperty(Interior.DEVICE_CELLULAR_SOURCE_PIN,role)==actual ||
        error("loaded DeviceCellular source differs from this experiment checkout: $role")
end

kind=get(ENV,"DEVICE_BACKEND","amd")
get(ENV,"DEVICE_AMD_OPTLEVEL","2") in ("1","2") || error("DEVICE_AMD_OPTLEVEL must be1or2")
if kind=="amd"
    using AMDGPU
    backend=ROCBackend()
    if get(ENV,"DEVICE_AMD_OPTLEVEL","2")=="1"
        # Pinned ROCm runtime LLVM18 EarlyCSE/MemorySSA crash workaround.
        # Process-local compiler config only; no package/runtime files are modified.
        config=AMDGPU.compiler_config(AMDGPU.device();maxthreads=128)
        lock(AMDGPU.Compiler.compiler_config_lock) do
            for key in keys(AMDGPU.Compiler._compiler_configs)
                if AMDGPU.Compiler._compiler_configs[key]===config
                    AMDGPU.Compiler._compiler_configs[key]=AMDGPU.GPUCompiler.CompilerConfig(config;opt_level=1)
                end
            end
        end
        @assert AMDGPU.compiler_config(AMDGPU.device();maxthreads=128).opt_level==1
    end
elseif kind=="metal"
    using Metal
    backend=MetalBackend()
elseif kind=="cpu"
    backend=CPU()
else
    error("DEVICE_BACKEND must be amd|metal|cpu")
end

function device_cellular_runtime()
    project=Base.active_project();manifest=joinpath(dirname(project),"Manifest.toml")
    gpu=kind=="amd" ? string(AMDGPU.device()) : kind=="metal" ? string(Metal.device()) : "CPU"
    backend_version=kind=="amd" ? string(Base.pkgversion(AMDGPU)) : kind=="metal" ? string(Base.pkgversion(Metal)) : "KernelAbstractions"
    actual_optlevel=kind=="amd" ? Int(AMDGPU.compiler_config(AMDGPU.device();maxthreads=128).opt_level) : nothing
    (julia=string(VERSION),kernel=string(Sys.KERNEL),architecture=string(Sys.ARCH),
        launch_threads=Threads.nthreads(),backend=kind,device=gpu,backend_version=backend_version,
        actual_device_optlevel=actual_optlevel,
        project_sha256=bytes2hex(Interior.SHA.sha256(read(project))),
        manifest_sha256=isfile(manifest) ? bytes2hex(Interior.SHA.sha256(read(manifest))) : nothing)
end
