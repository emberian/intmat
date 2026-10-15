# Tiny CPU/device gate before a separately authorized native-volume study.
# The launcher must own a bounded GPU queue slot. No CPU fallback bake occurs.
using Interior,Random,TOML,SHA
using Interior:JSON3,KernelAbstractions
include("living_tissue_turnover_run.jl")

function turnover_gate_fixture()
    rng=MersenneTwister(947);state=rand(rng,Float32,16,16,16,7)
    state[:,:,:,1].=.01f0.+.3f0.*state[:,:,:,1]
    state[:,:,:,2].=.05f0.+.8f0.*state[:,:,:,2]
    state[:,:,:,3].=.01f0.+.3f0.*state[:,:,:,3]
    state[:,:,:,4].*=.02f0
    state[:,:,:,5].=1.2f0.*(state[:,:,:,5].-.5f0)
    state[:,:,:,6:7].=.4f0.*(state[:,:,:,6:7].-.5f0)
    reserve=2f0.*rand(rng,Float32,16,16,16)
    for i in 1:11:length(reserve)
        state[i]=0;reserve[i]=0 # Include unfunded/partially funded donors.
    end
    state,reserve
end
function turnover_gpu_gate(recipe,output;backend_kind="amd")
    backend_kind in ("amd","metal") || error("native volume gate requires a GPU backend")
    backend=BASE.tissue_backend(backend_kind)
    Base.invokelatest(turnover_gpu_gate_backend,backend,recipe,output,backend_kind)
end
function turnover_gpu_gate_backend(backend,recipe,output,backend_kind)
    ispath(output)&&error("fresh output required")
    c=TOML.parsefile(recipe);base=turnover_base_parameters(;NamedTuple(Symbol(k)=>v for (k,v) in c["base"])...)
    phase=turnover_parameters(;NamedTuple(Symbol(k)=>v for (k,v) in c["phase"])...)
    host,reserve=turnover_gate_fixture();grid=Grid(16,16,16;extent=(6,8,10));dt=.02f0
    cpu=tissue_world(grid;seed=947);reference=TurnoverTissue(cpu;initial=host,reserve,base_parameters=base,parameters=phase)
    gpu=tissue_world(grid;backend,seed=947);device=TurnoverTissue(gpu;initial=host,reserve,base_parameters=base,parameters=phase)
    initial=turnover_statistics(cpu,reference)
    mkdir(output);cp(@__FILE__,joinpath(output,"gate-source.jl"));cp(recipe,joinpath(output,"recipe.toml"))
    for _ in 1:50;step!(cpu,reference,dt);step!(gpu,device,dt);end
    KernelAbstractions.synchronize(backend)
    a=Array(reference.tissue.state);b=Array(device.tissue.state)
    ra=Array(reference.reserve);rb=Array(device.reserve)
    state_error=maximum(abs.(a.-b));reserve_error=maximum(abs.(ra.-rb))
    cpu_stats=turnover_statistics(cpu,reference;state=a,reserve=ra)
    gpu_stats=turnover_statistics(gpu,device;state=b,reserve=rb)
    passed=state_error<1e-4 && reserve_error<1e-4 && gpu_stats.finite && gpu_stats.physical_min>=0 &&
        abs(gpu_stats.physical_material-initial.physical_material)<1e-5initial.physical_material
    record_state(cpu,reference,joinpath(output,"gate-cpu");step=50)
    record_state(gpu,device,joinpath(output,"gate-device");step=50)
    checkpoint=joinpath(output,"gate-checkpoint.jls");turnover_checkpoint(gpu,device,checkpoint)
    restored_world,restored=restore_turnover!(checkpoint;backend)
    exact_restore=Array(restored.tissue.state)==b && Array(restored.reserve)==rb && restored_world.time==gpu.time
    for _ in 1:11;step!(gpu,device,dt);step!(restored_world,restored,dt);end
    KernelAbstractions.synchronize(backend)
    exact_continuation=Array(restored.tissue.state)==Array(device.tissue.state) && Array(restored.reserve)==Array(device.reserve)
    gate=(passed=passed&&exact_restore&&exact_continuation,backend=backend_kind,dimensions=grid.size,extent=grid.extent,
        steps=50,dt,state_maxabs_error=state_error,reserve_maxabs_error=reserve_error,exact_restore,exact_continuation,
        initial=initial,reference=cpu_stats,device=gpu_stats,base_parameters=base,turnover_parameters=phase,
        study_sha256=LivingTissueTurnover.TURNOVER_SHA,base_sha256=BASE.MODEL_SHA,shared_sha256=BASE.SHARED_SHA,
        gate_sha256=bytes2hex(sha256(read(@__FILE__))),
        interpretation="Actual two-stage turnover law, seven tissue channels plus finite reserve. Dense randomized3D fixture includes zero-precursor/zero-reserve donors. No optical or biological interpretation follows from numerical parity.")
    open(io->JSON3.pretty(io,gate),joinpath(output,"gate.json"),"w")
    println(JSON3.write(gate));flush(stdout)
    gate.passed || error("turnover device parity/continuation failed; no native bake started")
    cpu=nothing;reference=nothing;gpu=nothing;device=nothing;restored_world=nothing;restored=nothing
    host=nothing;reserve=nothing;a=nothing;b=nothing;ra=nothing;rb=nothing;GC.gc(true)
    run_turnover_backend(recipe,joinpath(output,"volume"),backend,backend_kind)
    bytes=sum((filesize(joinpath(root,name)) for (root,_,files) in walkdir(output) for name in files);init=0)
    bytes<=3_000_000_000 || error("completed artifact exceeds3GB; preserved for inspection")
    open(io->JSON3.pretty(io,(artifact_bytes=bytes,budget_bytes=3_000_000_000,gate=gate.passed)),joinpath(output,"budget.json"),"w")
    println("VOLUME READY ",output," bytes=",bytes)
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS) in (2,3)||error("recipe fresh-output [amd|metal]")
    turnover_gpu_gate(ARGS[1],ARGS[2];backend_kind=length(ARGS)==3 ? ARGS[3] : "amd")
end
