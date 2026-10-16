#!/usr/bin/env julia
# Isolated CPU/AMD complete-step gate; use the verified hbox launcher environment.
include(joinpath(@__DIR__,"..","scripts","device_cellular_setup.jl"))
include("living_step_performance.jl")
function collect_arrays(x)
    d=Dict{String,Any}()
    groups=(cells=x.p.cells.state,environment=x.p.cells.environment,
        habitat=x.p.habitat.arrays,feed=x.p.habitat.feedstock.arrays,
        private=x.h.state,events=x.h.events)
    for group in keys(groups),key in keys(groups[group])
        d["$group.$key"]=Array(groups[group][key])
    end
    work=Interior.device_resource_workspace(x.p.cells)
    d["workspace.grants"]=Array(work.grants)
    for (i,level) in enumerate(work.levels),key in (:high,:low)
        d["workspace.level$i.$key"]=Array(level[key])
    end
    d["workspace.highwater"]=Array(Interior.device_lifecycle_workspace(x.p.cells).highwater)
    d
end
function full_compare(a,b;exact)
    count=0;errors=Dict{String,Float64}()
    a.w.time==b.w.time && a.h.steps==b.h.steps && a.p.cells.attempted_steps==b.p.cells.attempted_steps || error("clock mismatch")
    aa=collect_arrays(a);bb=collect_arrays(b)
    for key in keys(aa)
        av,bv=aa[key],bb[key];count+=length(av)
        if !isequal(av,bv)
            exact && error("exact mismatch $key")
            eltype(av)<:AbstractFloat || error("integer/identity mismatch $key")
            errorvalue=maximum(abs.(Float64.(av).-Float64.(bv));init=0.)
            tolerance=eltype(av)==Float64 ? 1e-10 : 2e-6
            isfinite(errorvalue)&&errorvalue<=tolerance||error("numeric mismatch $key: $errorvalue > $tolerance")
            errors[key]=errorvalue
        end
    end
    (;checked_array_elements=count,errors)
end
function clone_original(original,physical)
    kind=="cpu" && return deepcopy(original)
    w=World(original.w.grid,backend,NamedTuple(),NamedTuple(),original.w.time,original.w.seed)
    p=Interior.device_habitat_restore(w,physical)
    h=N.sidecar(p;dt=.02,horizon=2.,act_time=100.,inheritance_fraction=.5)
    h.options=original.h.options;h.steps=original.h.steps
    for key in keys(h.state);copyto!(h.state[key],original.h.state[key]);end
    for key in keys(h.events);copyto!(h.events[key],original.h.events[key]);end
    copyto!(h.original_masks,original.h.original_masks)
    (;w,p,h)
end
function native_run(out)
    rows=Any[];n=64;dt=.02
    cases=get(ENV,"LIVING_GATE_CASES","all")=="control" ? ((24,:full),) : ((24,:native),(46,:native),(24,:full))
    for (founders,coverage) in cases
        println((;stage="fixture",kind,founders,coverage));flush(stdout)
        cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","native","living-patch-recruitment-16.toml"))
        cfg["grid"]["shape"]=fill(n,3);cfg["cellular"]["founders"]=founders
        recipe=tempname()*".toml";open(io->TOML.print(io,cfg),recipe,"w")
        f=N.fixture(recipe;backend);h=N.sidecar(f.p;dt,horizon=2.,act_time=100.,inheritance_fraction=.5)
        original=(w=f.w,p=f.p,h=h);advance(original,false,220)
        sum(Array(original.h.state.actor_updates))>0 || error("fixture did not exercise paid actor updates")
        if coverage==:full
            bits=Array(original.p.habitat.arrays.enzyme_bits)
            values=reinterpret(Float32,bits)
            values .= max.(values,1f-7) # never reduce existing sorptive capacity
            copyto!(original.p.habitat.arrays.enzyme_bits,bits)
        end
        physical=kind=="cpu" ? nothing : Interior.device_habitat_checkpoint(original.w,original.p)
        a=clone_original(original,physical);b=clone_original(original,physical)
        advance(a,true,5);advance(b,false,5)
        checks=Any[full_compare(a,b;exact=kind=="cpu")]
        samples=Any[]
        for repeat in 1:5
            a=clone_original(original,physical);b=clone_original(original,physical)
            order=isodd(repeat) ? (true,false) : (false,true)
            timing=Dict{Bool,Any}()
            for baseline in order
                x=baseline ? a : b
                GC.gc();KernelAbstractions.synchronize(backend)
                measured=@timed advance(x,baseline,30)
                measured.compile_time==0 || error("timed block included compilation")
                timing[baseline]=(;seconds=measured.time,gc_seconds=measured.gctime,bytes=measured.bytes)
            end
            push!(checks,full_compare(a,b;exact=kind=="cpu"))
            push!(samples,(;baseline=timing[true],optimized=timing[false],order))
        end
        before=@timed collect_arrays(b)
        origin=N.C.A.audit(b.p)
        N.frame(b.w,b.p,b.h,origin)
        observer=@timed N.frame(b.w,b.p,b.h,origin)
        material=Array(original.p.habitat.arrays.enzyme_bits);bound=Array(original.p.habitat.arrays.bound)
        occupied=count(i->material[i,1]!=0||material[i,2]!=0||bound[i,1]!=0||bound[i,2]!=0,axes(bound,1))
        base=[s.baseline.seconds for s in samples];opt=[s.optimized.seconds for s in samples]
        row=(;resolution=n,founders,coverage,actor_updates=sum(Array(b.h.state.actor_updates)),learning_work=sum(Array(b.h.state.paid)),births=Array(b.p.cells.environment.metadata)[3],population=count(==(1),Array(b.p.cells.state.active)),occupied_fraction=occupied/n^3,
            baseline_step_seconds=median(base)/30,optimized_step_seconds=median(opt)/30,speedup=median(base)/median(opt),
            samples,checks,transfer_seconds=before.time,observer_seconds=observer.time,observer_origin_included=false,full_coverage_seed="max(existing_catalyst,1e-7)")
        push!(rows,row)
        report=(;schema="living-step-native-gate-v1",runtime=device_cellular_runtime(),baseline=BASELINE,cellular_source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
            kernel_sha256=bytes2hex(N.sha256(read(joinpath(@__DIR__,"..","src","device_habitat_kernels.jl")))),
            resource_pin="CPUs10–11;Julia2threads;OpenBLAS1;8GiB;no swap;nice15",
            cpu_affinity=only(filter(l->startswith(l,"Cpus_allowed_list:"),split(read("/proc/self/status",String),'\n'))),
            nice=ccall(:getpriority,Cint,(Cint,Cint),0,0),cgroup=read("/proc/self/cgroup",String),rows)
        open(io->JSON3.write(io,report),out,"w")
        println((;founders,coverage,speedup=row.speedup,baseline_ms=1000row.baseline_step_seconds,optimized_ms=1000row.optimized_step_seconds));flush(stdout)
        original=nothing;a=nothing;b=nothing;GC.gc()
    end
end
native_run(ARGS[1])
