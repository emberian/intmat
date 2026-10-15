# Bounded Julia CPU/device comparison and complete-step timing. The caller must
# reserve the backend through the host coordinator before launching this study.
using Interior, Statistics, Random
using Interior: JSON3, KernelAbstractions
include("living_tissue.jl")
using .LivingTissueStudy

function tissue_dense_fixture(side;seed=913)
    rng=MersenneTwister(seed);state=rand(rng,Float32,side,side,side,7)
    state[:,:,:,1].=.05f0.+state[:,:,:,1]
    state[:,:,:,2].=.02f0.+.6f0.*state[:,:,:,2]
    state[:,:,:,3].=.01f0.+.3f0.*state[:,:,:,3]
    state[:,:,:,4].*=.02f0
    state[:,:,:,5:7].=.6f0.*(state[:,:,:,5:7].-.5f0)
    state
end

function tissue_gpu_probe(backend,output;side=128,steps=256)
    0<side<=256 && 0<steps<=10000 || error("bounded probe requires side 1:256 and steps 1:10000")
    ispath(output)&&error("fresh output required")
    q=tissue_parameters(taxis=.028,taxis_sensitivity=16.,barrier=8.,neural_gain=2.4)
    a=tissue_world(Grid(16,16,16;extent=(16,16,16));seed=913)
    p=LivingTissue(a;parameters=q,initial=tissue_dense_fixture(16))
    b=tissue_world(a.grid;backend,seed=913)
    r=LivingTissue(b;initial=Array(p.state),parameters=q)
    for _ in 1:50;step!(a,p,.04);step!(b,r,.04);end
    KernelAbstractions.synchronize(backend)
    device=Array(r.state);reference=Array(p.state)
    error=maximum(abs.(device.-reference))
    error<1e-4 || Base.error("CPU/device mismatch $error")
    fixture=(steps=50,dt=.04,maxabs=error,reference=tissue_statistics(a,p),device=tissue_statistics(b,r))
    cases=Any[]
    for name in ("sparse-lace","dense-tissue")
        w=tissue_world(Grid(side,side,side;extent=(side,side,side));backend,seed=913)
        initial=name=="dense-tissue" ? tissue_dense_fixture(side) : nothing
        tissue=LivingTissue(w;parameters=q,initial,pattern="lace");initial=nothing
        for _ in 1:10;step!(w,tissue,.04);end
        KernelAbstractions.synchronize(backend)
        measurements=Float64[]
        for trial in 1:3
            elapsed=@elapsed begin
                for _ in 1:steps;step!(w,tissue,.04);end
                KernelAbstractions.synchronize(backend)
            end
            push!(measurements,elapsed)
        end
        push!(cases,(name,timed_seconds=measurements,median_seconds=median(measurements),
            voxel_steps_per_second=side^3*steps/median(measurements),final=tissue_statistics(w,tissue)))
        w=nothing;tissue=nothing;GC.gc(true)
    end
    report=(backend=string(typeof(backend)),side,steps,fixture,warm_steps=10,
        cases,resident_array_bytes=56side^3,
        source_sha256=LivingTissueStudy.MODEL_SHA,shared_source_sha256=LivingTissueStudy.SHARED_SHA,
        boundary="Complete step! calls with all seven physical/neural channels, validation, buffer swaps and final synchronization. Excludes initialization, compilation, capture and rendering. Distinct continuum law; no cellular-arena speedup inferred.")
    open(io->JSON3.pretty(io,report),output,"w");println(JSON3.write(report))
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS) in (2,3,4)||error("cpu|amd|metal fresh-output.json [side] [steps]")
    backend=LivingTissueStudy.tissue_backend(ARGS[1])
    Base.invokelatest(tissue_gpu_probe,backend,ARGS[2];side=length(ARGS)>2 ? parse(Int,ARGS[3]) : 128,
        steps=length(ARGS)>3 ? parse(Int,ARGS[4]) : 256)
end
