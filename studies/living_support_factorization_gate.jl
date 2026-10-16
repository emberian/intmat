#!/usr/bin/env julia
include("living_gpu_profile.jl")
include("living_support_factorization.jl")
const F=LivingSupportFactorization
function gate(out)
    rows=Any[];steps=20
    for founders in (24,46)
        println((;stage="warm",founders));flush(stdout)
        cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","native","living-patch-recruitment-16.toml"))
        cfg["grid"]["shape"]=fill(64,3);cfg["cellular"]["founders"]=founders
        path=tempname()*".toml";open(io->TOML.print(io,cfg),path,"w")
        f=N.fixture(path;backend);h=N.sidecar(f.p;dt=.02,horizon=2.,act_time=100.,inheritance_fraction=.5,launcher=PolicyLaunch(nothing))
        original=(;w=f.w,p=f.p,h);advance!(original,220)
        sum(Array(h.state.actor_updates))>0||error("no paid actor")
        physical=Interior.device_habitat_checkpoint(original.w,original.p)
        warm=clone(original,physical);advance!(warm,steps;owner=nothing)
        warm=clone(original,physical);advance!(warm,steps;owner=F.Carrier(warm.p.cells))
        warm=nothing;GC.gc();samples=Any[];checks=Int[]
        for repeat in 1:5
            a=clone(original,physical);b=clone(original,physical);owner=F.Carrier(b.p.cells);timings=Dict{Symbol,Any}()
            for mode in (isodd(repeat) ? (:baseline,:candidate) : (:candidate,:baseline))
                x=mode==:baseline ? a : b;carrier=mode==:baseline ? nothing : owner
                GC.gc();KernelAbstractions.synchronize(backend)
                t=@timed advance!(x,steps;owner=carrier)
                t.compile_time==0||error("timed compilation $(t.compile_time)")
                timings[mode]=(;seconds=t.time,enqueue_seconds=t.value.enqueue_seconds,wait_seconds=t.value.wait_seconds,allocations_bytes=t.bytes,gc_seconds=t.gctime)
            end
            push!(checks,parity(a,b));push!(samples,timings)
        end
        base=median([s[:baseline].seconds for s in samples])/steps;candidate=median([s[:candidate].seconds for s in samples])/steps
        row=(;founders,population=count(==(1),Array(b.p.cells.state.active)),births=Array(b.p.cells.environment.metadata)[3],actor_updates=sum(Array(b.h.state.actor_updates)),steps,
            baseline_step_seconds=base,candidate_step_seconds=candidate,speedup=base/candidate,samples,checks,
            extra_scratch_bytes=sizeof(owner.sites)+sizeof(owner.weights)+sizeof(owner.totals))
        push!(rows,row)
        open(io->JSON3.write(io,(;schema="living-support-factorization-v1",runtime=device_cellular_runtime(),
            resources="CPU10–11 SMT siblings / one P-core; Julia2;8GiB;no swap;nice15; exclusive GPU",source=Interior.DEVICE_CELLULAR_SOURCE_PIN,
            semantics="parallel candidate geometry; stable ordered compaction and Float32 sum; parallel normalization; original compressed support layout and tails",
            study_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"living_support_factorization.jl")))),rows)),out,"w")
        println((;stage="done",founders,baseline_ms=1000base,candidate_ms=1000candidate,speedup=row.speedup));flush(stdout)
        original=nothing;physical=nothing;a=nothing;b=nothing;owner=nothing;GC.gc()
    end
end
gate(ARGS[1])
