using Interior,Random,Statistics
using Interior: JSON3
isdefined(@__MODULE__,:DeviceMemorySpecialization) || include("device_memory_specialization.jl")
include("../scripts/device_habitat_pilot_fixture.jl")
const M=DeviceMemorySpecialization
function run_memory_specialization()
    Threads.nthreads()==1 || error("CPU page prototype requires one Julia thread")
    rows=Any[]
    for variant in ("closed","finite_sources"), pages in (false,true)
        w,p=device_habitat_pilot_fixture(CPU();size=16,founders=8,capacity=64,variant)
        baseline=M.retained_array_bytes(p)
        world,copy_p=deepcopy((w,p));optimized=M.specialize(world,copy_p;pages)
        before=M.retained_array_bytes(optimized)
        exact=true;reference_times=Float64[];specialized_times=Float64[]
        for i in 1:16
            a=@elapsed step!(w,p,.002);b=@elapsed step!(world,optimized,.002)
            i>3 && (push!(reference_times,a);push!(specialized_times,b))
            exact &= isequal(M.persistent_snapshot(p),M.persistent_snapshot(optimized))
        end
        exact||error("persistent coupled state diverged")
        # Backing arrays are counted once; Dict/page metadata is additionally noted.
        arrays=optimized.habitat.arrays
        ledger=arrays.grid_ledger.data
        occupancy=ledger isa M.PagedMatrix ? length(ledger.pages)*ledger.page_rows/size(ledger,1) : 1.
        push!(rows,(;variant,pages,shape=w.grid.size,capacity=64,steps=16,baseline_bytes=baseline,
            specialized_initial_bytes=before,specialized_final_bytes=M.retained_array_bytes(optimized),specialized_host_summary_bytes=Base.summarysize(optimized),
            persistent_state_exact=exact,ledger_page_coverage=occupancy,reference_median_seconds=median(reference_times),specialized_median_seconds=median(specialized_times),
            base_voxel_bytes=84,habitat_voxel_bytes=132,feed_voxel_bytes=52,
            transient_contents_excluded=M.TRANSIENT_ENV))
    end
    report=(;schema="device-memory-specialization-cpu",prototype_sha=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"device_memory_specialization.jl")))),
        scope="existing coupled numerical kernels and schedule; CPU1 storage-only",julia=string(VERSION),threads=Threads.nthreads(),cellular_source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,gpu_validated=false,rows)
    println(JSON3.write(report));report
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    result=run_memory_specialization();!isempty(ARGS)&&open(io->JSON3.pretty(io,result),ARGS[1],"w")
end
