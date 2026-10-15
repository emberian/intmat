# Run only in an explicitly granted bounded GPU slot. The CPU backend is exposed
# only by the inner function for tiny harness-plumbing checks, not the CLI.
using Interior,Random,TOML,SHA
include("living_tissue_packing_run.jl")
include("living_tissue_device_gate.jl")

function packing_gate_fixture(capacity)
    rng=MersenneTwister(1947);state=rand(rng,Float32,16,16,16,7)
    rho=(.05f0.+.85f0.*rand(rng,Float32,16,16,16)).*capacity
    partition=rand(rng,Float32,16,16,16)
    state[:,:,:,1].=.01f0.+.3f0.*state[:,:,:,1]
    state[:,:,:,2].=rho.*partition;state[:,:,:,3].=rho.*(1f0.-partition)
    state[:,:,:,4].*=.02f0
    state[:,:,:,5].=1.2f0.*(state[:,:,:,5].-.5f0)
    state[:,:,:,6:7].=.4f0.*(state[:,:,:,6:7].-.5f0)
    reserve=2f0.*rand(rng,Float32,16,16,16)
    for i in 1:11:length(reserve);state[i]=0;reserve[i]=0;end
    state[8,8,8,1]=8;state[8,8,8,2:3].=0
    for I in ((7,8,8),(9,8,8),(8,7,8),(8,9,8),(8,8,7),(8,8,9))
        i,j,k=I;state[i,j,k,1]=1;state[i,j,k,2]=.7f0*capacity
        state[i,j,k,3]=capacity-state[i,j,k,2];state[i,j,k,5]=-.9f0
    end
    state,reserve
end
function packing_gpu_gate(recipe,output;backend_kind="amd",bake=true)
    backend_kind in ("amd","metal") || error("device gate CLI requires a GPU backend")
    backend=BASE.tissue_backend(backend_kind)
    Base.invokelatest(packing_gpu_gate_backend,backend,recipe,output,backend_kind;bake)
end
function packing_gpu_gate_backend(backend,recipe,output,backend_kind;bake=true)
    ispath(output) && error("fresh output required")
    c=TOML.parsefile(recipe)
    base=packing_base_parameters(;NamedTuple(Symbol(k)=>v for (k,v) in c["base"])...)
    phase=packing_parameters(;NamedTuple(Symbol(k)=>v for (k,v) in c["phase"])...)
    host,reserve=packing_gate_fixture(phase.capacity);grid=Grid(16,16,16;extent=(6,8,10))
    cpu=LivingTissuePacking.tissue_world(grid;seed=1947)
    reference=PackingTissue(cpu;initial=host,reserve,base_parameters=base,parameters=phase)
    gpu=LivingTissuePacking.tissue_world(grid;backend,seed=1947)
    device=PackingTissue(gpu;initial=host,reserve,base_parameters=base,parameters=phase)
    gate=tissue_device_gate(cpu,reference,gpu,device,output;dt=.005f0,
        measure=packing_statistics,save_checkpoint=packing_checkpoint,restore=restore_packing!,backend_kind,
        provenance=(packing_sha256=LivingTissuePacking.PACKING_SHA,
            turnover_sha256=LivingTissuePacking.LT.TURNOVER_SHA,base_sha256=BASE.MODEL_SHA,
            shared_sha256=BASE.SHARED_SHA,base_parameters=base,packing_parameters=phase,
            wrapper_sha256=bytes2hex(sha256(read(@__FILE__))),
            fixture="Dense randomized3D with mixed B/M occupancy, zero precursor/reserve donors, six full mixed-species neighbors competing for one empty center and anisotropic physical spacing."))
    cp(@__FILE__,joinpath(output,"packing-gate-source.jl"));cp(recipe,joinpath(output,"recipe.toml"))
    gate.passed || error("packing gate failed; no requested bake started")
    cpu=nothing;reference=nothing;gpu=nothing;device=nothing;host=nothing;reserve=nothing;GC.gc(true)
    bake && run_packing_backend(recipe,joinpath(output,"volume"),backend,backend_kind)
    bytes=sum((filesize(joinpath(root,name)) for (root,_,files) in walkdir(output) for name in files);init=0)
    budget=get(c,"capture_budget_bytes",3_000_000_000)
    bytes<=budget || error("packing artifacts exceed declared budget; preserved")
    open(io->JSON3.pretty(io,(artifact_bytes=bytes,budget_bytes=budget,gate=gate.passed,bake)),joinpath(output,"budget.json"),"w")
    println("PACKING READY ",output," bytes=",bytes)
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS) in (2,3) || error("recipe fresh-output [amd|metal]")
    packing_gpu_gate(ARGS[1],ARGS[2];backend_kind=length(ARGS)==3 ? ARGS[3] : "amd")
end
