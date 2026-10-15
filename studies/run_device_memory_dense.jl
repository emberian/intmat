# Native storage gate: existing numerical processes, no algorithm replacement.
include(joinpath(@__DIR__,"../scripts/device_cellular_setup.jl"))
using Statistics
include("device_memory_dense.jl")
include(joinpath(@__DIR__,"../scripts/device_habitat_pilot_fixture.jl"))
const DMD=DeviceMemoryDense
const KA=Interior.KernelAbstractions
function vram_used()
    kind=="amd" || return nothing
    paths=filter(isfile,["/sys/class/drm/card$i/device/mem_info_vram_used" for i in 0:8])
    isempty(paths) ? nothing : sum(parse(Int,strip(read(p,String))) for p in paths)
end
function run_dense_gate()
    n=parse(Int,get(ENV,"DENSE_MEMORY_SIZE","16"));n<=128||error("bounded gate maximum128")
    samples=Any[]
    for variant in ("closed","finite_sources")
        host,hp=device_habitat_pilot_fixture(CPU();size=n,founders=8,capacity=64,variant)
        for _ in 1:3;step!(host,hp,.002);end
        cp=Interior.device_habitat_checkpoint(host,hp)
        emptyworld()=World(host.grid,backend,NamedTuple(),NamedTuple(),0.,host.seed)
        wa=emptyworld();a=Interior.device_habitat_restore(wa,cp)
        wb=emptyworld();b=DMD.restore(wb,cp)
        KA.synchronize(backend)
        resident_before=vram_used();ta=Float64[];tb=Float64[]
        for i in 1:10
            elapsed=@elapsed begin step!(wa,a,.002);KA.synchronize(backend);end
            i>3&&push!(ta,elapsed)
            elapsed=@elapsed begin step!(wb,b,.002);KA.synchronize(backend);end
            i>3&&push!(tb,elapsed)
            isequal(DMD.M.persistent_snapshot(a),DMD.M.persistent_snapshot(b))||error("persistent state mismatch $variant/$i")
        end
        saved=Interior.device_habitat_checkpoint(wb,b)
        wc=emptyworld();c=Interior.device_habitat_restore(wc,saved)
        for i in 1:3
            step!(wb,b,.002);step!(wc,c,.002)
            isequal(DMD.M.persistent_snapshot(b),DMD.M.persistent_snapshot(c))||error("restored future mismatch $variant/$i")
        end
        storage=x->x isa Array || (kind=="amd"&&x isa AMDGPU.ROCArray) || (kind=="metal"&&x isa Metal.MtlArray)
        push!(samples,(;variant,size=n,capacity=64,persistent_exact=true,restored_future_exact=true,
            baseline_bytes=DMD.retained_bytes(a,storage),specialized_bytes=DMD.retained_bytes(b,storage),
            baseline_median_seconds=median(ta),specialized_median_seconds=median(tb),
            vram_pair_before=resident_before,vram_three_after=vram_used(),
            memory_scope="reachable backing arrays counted once; VRAM includes runtime/pools and simultaneous reference/restored carriers"))
        a=nothing;b=nothing;c=nothing;GC.gc(true);KA.synchronize(backend)
    end
    helpers=("run_device_memory_dense.jl","device_memory_dense.jl","device_memory_specialization.jl","../scripts/device_cellular_setup.jl","../scripts/device_habitat_pilot_fixture.jl")
    pins=Dict(p=>bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,p)))) for p in helpers)
    report=(;schema="dense-memory-native-gate",runtime=device_cellular_runtime(),sources=DMD.SOURCE_PIN,helpers=pins,samples,
        direct_device_allocation=true,host_dense_restore_peak_remaining=true,paging_gpu=false)
    open(io->JSON3.pretty(io,report),ARGS[1],"w");report
end
run_dense_gate()
