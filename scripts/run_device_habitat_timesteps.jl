#!/usr/bin/env julia
# Launch with AMDGPU preloaded before setup to avoid late-import world-age issues.
include("device_cellular_setup.jl")
include(joinpath(@__DIR__,"..","studies","device_habitat_timesteps.jl"))
function timestep_main()
    out=only(ARGS);ispath(out)&&error("fresh output required")
    mkpath(dirname(abspath(out)));mkdir(out)
    started=time_ns()
    result=habitat_timestep_study(backend;size=parse(Int,get(ENV,"TIMESTEP_SIZE","128")),
        founders=parse(Int,get(ENV,"TIMESTEP_FOUNDERS","512")),
        capacity=parse(Int,get(ENV,"TIMESTEP_CAPACITY","8192")),
        reference_dt=parse(Float64,get(ENV,"TIMESTEP_REFERENCE_DT",".0005")),
        reference_steps=parse(Int,get(ENV,"TIMESTEP_REFERENCE_STEPS","32")),
        factors=Tuple(parse.(Int,split(get(ENV,"TIMESTEP_FACTORS","1,2,4,8,16"),','))),
        support_strategy=Symbol(get(ENV,"TIMESTEP_SUPPORT","axis")),
        repetitions=parse(Int,get(ENV,"TIMESTEP_REPETITIONS","1")),
        reference_repetitions=parse(Int,get(ENV,"TIMESTEP_REFERENCE_REPETITIONS",get(ENV,"TIMESTEP_REPETITIONS","1"))),
        warmup=get(ENV,"TIMESTEP_WARMUP","false")=="true")
    whole_seconds=(time_ns()-started)/1e9
    report=(source=Interior.DEVICE_CELLULAR_SOURCE_PIN,runtime=device_cellular_runtime(),
        study_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"..","studies","device_habitat_timesteps.jl")))),
        driver_sha256=bytes2hex(Interior.SHA.sha256(read(@__FILE__))),
        whole_seconds,timing="per-candidate timer includes first-use compile and batch8 status synchronization; whole includes construction and full final downloads; not warmed throughput",result)
    open(joinpath(out,"timesteps.json"),"w") do io;JSON3.pretty(io,report);end
    println("completed fixed-duration native128³ timestep gate: ",whole_seconds," seconds")
end
timestep_main()
