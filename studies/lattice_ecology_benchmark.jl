# Reproducible CPU1 complete-step allocation/throughput probe; no GPU launch.
using Interior,SHA
using Interior:JSON3
include("lattice_ecology.jl")
using .LatticeEcologyStudy
function lattice_benchmark_main(output)
    Threads.nthreads()==1||error("run this bounded reference with JULIA_NUM_THREADS=1")
    ispath(output)&&error("choose a new result path");mkpath(dirname(abspath(output)))
    rows=Any[]
    for n in (8,16,32,64)
        w=lattice_world(Grid(n,n,1;extent=(n,n,1)));p=LatticeEcology(w;occupancy=.5)
        for _ in 1:64;step!(w,p,.1);end
        # Exact same warmed future-determining state for allocation and time.
        v,q=deepcopy((w,p));initial=lattice_statistics(w,p)
        bytes=@allocated for _ in 1:100;step!(w,p,.1);end
        elapsed=@elapsed for _ in 1:100;step!(v,q,.1);end
        @assert v.time==w.time&&q.cells==p.cells&&q.genotype==p.genotype&&q.events==p.events
        final=lattice_statistics(w,p)
        push!(rows,(shape=(n,n,1),steps=100,dt=.1,warmup_steps=64,
            complete_step_allocated_bytes=bytes/100,complete_step_seconds=elapsed/100,
            voxel_steps_per_second=n*n*100/elapsed,population_initial=initial.population,population_final=final.population,
            resident_array_bytes=final.resident_bytes,shared_genotype_bytes=final.shared_genotype_bytes,
            hypothetical_per_voxel_genome_bytes=final.hypothetical_per_voxel_genome_bytes))
    end
    result=(kind="CPU1 complete physics step, separate identical-state allocation and timing windows",
        julia=string(VERSION),threads=Threads.nthreads(),model_source_sha256=LatticeEcologyStudy.MODEL_SHA,
        benchmark_source_sha256=bytes2hex(sha256(read(@__FILE__))),shared_sources=Interior.arena_sources(),rows=rows)
    open(output,"w") do io;JSON3.write(io,result);end
    println(JSON3.write(result))
end
if abspath(PROGRAM_FILE)==(@__FILE__)||(isdefined(Main,:INTERIOR_RUNNING_STUDY)&&Main.INTERIOR_RUNNING_STUDY=="lattice_ecology_benchmark")
    length(ARGS)==1||error("fresh-results.json");lattice_benchmark_main(first(ARGS))
end
