#!/usr/bin/env julia
# Intrusive HIP-event attribution of a complete, paid-learning living trajectory.
include(joinpath(@__DIR__,"..","scripts","device_cellular_setup.jl"))
using KernelAbstractions, TOML, Statistics
include("native_living_construction.jl")
const N=NativeLivingConstruction
Base.include(Interior,joinpath(@__DIR__,"device_launch_profile.jl"))
mutable struct Events
    pool::Any
    labels::Vector{Symbol}
    cursor::Int
end
function observe!(f,o,kernel)
    o===nothing && return f()
    o.cursor+=1;i=o.cursor;i<=length(o.pool)||error("event pool exhausted")
    push!(o.labels,nameof(kernel));before,after=o.pool[i]
    AMDGPU.HIP.record(before);f();AMDGPU.HIP.record(after)
    nothing
end
struct PolicyLaunch
    observer::Any
end
N.L.policy_launch!(l::PolicyLaunch,h,k,args...;ndrange)=observe!(l.observer,k) do
    k(h.backend,64)(args...;ndrange)
end
struct Carrier{P}
    process::P
    observer::Any
end
Base.getproperty(p::Carrier,k::Symbol)=k in (:process,:observer) ? getfield(p,k) : getproperty(getfield(p,:process),k)
Interior.device_support_descriptor(p::Carrier)=Interior.device_support_descriptor(p.process)
function Interior.device_cellular_launch!(p::Carrier,k,args...;ndrange)
    observe!(p.observer,k) do
        Interior.device_cellular_launch!(p.process,k,args...;ndrange)
    end
end
function advance!(x,steps;observer=nothing,owner=nothing)
    enqueue_started=time_ns()
    for _ in 1:steps
        p=x.p.cells;Interior.device_cellular_validate(x.w,p);p.attempted_steps+=1
        Interior.device_cellular_enqueue!(x.w,owner===nothing ? Carrier(p,observer) : owner,Float32(.02),:projected,x.h,UInt32(p.attempted_steps))
    end
    enqueue_seconds=(time_ns()-enqueue_started)/1e9
    started=time_ns();Interior.device_cellular_status!(x.p.cells);wait_seconds=(time_ns()-started)/1e9
    x.w.time+=steps*.02
    (;enqueue_seconds,wait_seconds)
end
function clone(original,physical,observer=nothing)
    w=World(original.w.grid,backend,NamedTuple(),NamedTuple(),original.w.time,original.w.seed)
    p=Interior.device_habitat_restore(w,physical)
    h=N.sidecar(p;dt=.02,horizon=2.,act_time=100.,inheritance_fraction=.5,launcher=PolicyLaunch(observer))
    h.options=original.h.options;h.steps=original.h.steps
    for group in (:state,:events),k in keys(getproperty(h,group))
        copyto!(getproperty(h,group)[k],getproperty(original.h,group)[k])
    end
    copyto!(h.original_masks,original.h.original_masks)
    (;w,p,h)
end
function arrays(x)
    d=Dict{String,Any}()
    for (name,group) in pairs((cells=x.p.cells.state,environment=x.p.cells.environment,habitat=x.p.habitat.arrays,
            feed=x.p.habitat.feedstock.arrays,private=x.h.state,events=x.h.events))
        for k in keys(group);d["$name.$k"]=Array(group[k]);end
    end
    work=Interior.device_resource_workspace(x.p.cells);d["grants"]=Array(work.grants)
    for (i,level) in enumerate(work.levels),k in (:high,:low);d["level$i.$k"]=Array(level[k]);end
    d["highwater"]=Array(Interior.device_lifecycle_workspace(x.p.cells).highwater);d
end
function parity(a,b)
    aa=arrays(a);bb=arrays(b);n=0
    for k in keys(aa);isequal(aa[k],bb[k])||error("profile changed $k");n+=length(aa[k]);end
    a.w.time==b.w.time && a.h.steps==b.h.steps || error("clock mismatch");n
end
function run(out)
    rows=Any[];steps=8
    for founders in (24,46)
        println((;stage="warm",founders));flush(stdout)
        cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","native","living-patch-recruitment-16.toml"))
        cfg["grid"]["shape"]=fill(64,3);cfg["cellular"]["founders"]=founders
        path=tempname()*".toml";open(io->TOML.print(io,cfg),path,"w")
        f=N.fixture(path;backend);h=N.sidecar(f.p;dt=.02,horizon=2.,act_time=100.,inheritance_fraction=.5,launcher=PolicyLaunch(nothing))
        original=(;w=f.w,p=f.p,h);advance!(original,220)
        sum(Array(h.state.actor_updates))>0||error("no paid actor")
        physical=Interior.device_habitat_checkpoint(original.w,original.p)
        pool=Interior.device_profile_event_pool(AMDGPU,4096)
        warm=clone(original,physical,Events(pool,Symbol[],0));advance!(warm,steps;observer=warm.h.launcher.observer)
        warm=clone(original,physical);advance!(warm,steps;observer=nothing)
        warm=nothing;GC.gc()
        samples=Any[];lastcheck=0
        for repeat in 1:3
            plain=clone(original,physical);o=Events(pool,Symbol[],0);profile=clone(original,physical,o)
            timings=Dict{Symbol,Any}()
            for mode in (isodd(repeat) ? (:plain,:profile) : (:profile,:plain))
                x=mode==:plain ? plain : profile;obs=mode==:plain ? nothing : o
                GC.gc();KernelAbstractions.synchronize(backend)
                t=@timed advance!(x,steps;observer=obs)
                t.compile_time==0||error("timed compilation $(t.compile_time)")
                timings[mode]=(;wall_seconds=t.time,t.value...,allocations_bytes=t.bytes,gc_seconds=t.gctime)
            end
            lastcheck=parity(plain,profile)
            ks=Dict{Symbol,Any}()
            for i in 1:o.cursor
                key=o.labels[i];n,t=get(ks,key,(0,0.));ks[key]=(n+1,t+Float64(AMDGPU.HIP.elapsed(pool[i]...)))
            end
            labels=[(;kernel=k,category=Interior.device_profile_category(k),calls=v[1],seconds=v[2]) for (k,v) in ks]
            sort!(labels;by=x->x.seconds,rev=true)
            stream=Float64(AMDGPU.HIP.elapsed(pool[1][1],pool[o.cursor][2]))
            push!(samples,(;timings,launches=o.cursor,stream_seconds=stream,summed_event_seconds=sum(x.seconds for x in labels),labels))
        end
        row=(;founders,population=count(==(1),Array(original.p.cells.state.active)),births=Array(original.p.cells.environment.metadata)[3],
            paid_actor_updates=sum(Array(original.h.state.actor_updates)),steps,exact_elements=lastcheck,samples)
        push!(rows,row)
        open(io->JSON3.write(io,(;schema="living-gpu-profile-v1",runtime=device_cellular_runtime(),
            resource="CPU10–11 SMT siblings, one physical P-core; Julia2; nice15;8GiB/no swap;exclusive GPU",
            caveat="HIP event intervals are intrusive and include host submission gaps; compare paired native wall and enqueue/wait times. No observer or transfer inside step timings.",rows)),out,"w")
        println((;stage="done",founders,wall=[s.timings[:plain].wall_seconds/steps for s in samples],launches=samples[1].launches÷steps));flush(stdout)
        original=nothing;physical=nothing;pool=nothing;GC.gc()
    end
end
if abspath(PROGRAM_FILE)==(@__FILE__);run(ARGS[1]);end
