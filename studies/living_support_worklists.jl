module LivingSupportWorklists
using Interior, KernelAbstractions, TOML, Statistics
using Interior: JSON3
include("native_living_construction.jl")
const N=NativeLivingConstruction

mutable struct CPUKernelTimes
    seconds::Dict{Symbol,Float64}
    calls::Dict{Symbol,Int}
end
CPUKernelTimes()=CPUKernelTimes(Dict{Symbol,Float64}(),Dict{Symbol,Int}())
struct ProfileCarrier{P,T}
    process::P
    timer::T
end
@inline Base.getproperty(p::ProfileCarrier,k::Symbol)=k===:process || k===:timer ? getfield(p,k) : getproperty(getfield(p,:process),k)
Interior.device_support_descriptor(p::ProfileCarrier)=Interior.device_support_descriptor(p.process)
function Interior.device_cellular_launch!(p::ProfileCarrier,kernel,args...;ndrange)
    started=time_ns()
    Interior.device_cellular_launch!(p.process,kernel,args...;ndrange)
    elapsed=(time_ns()-started)/1e9;key=nameof(kernel)
    p.timer.seconds[key]=get(p.timer.seconds,key,0.)+elapsed
    p.timer.calls[key]=get(p.timer.calls,key,0)+1
end
function advance!(x,n;timer=nothing,dt=.02)
    for _ in 1:n
        Interior.device_cellular_validate(x.w,x.p.cells)
        delta=Float32(dt)
        0<delta<=Interior.device_cellular_timestep_limit(x.p.cells)||error("invalid timestep")
        x.p.cells.attempted_steps+=1
        carrier=timer===nothing ? x.p.cells : ProfileCarrier(x.p.cells,timer)
        Interior.device_cellular_enqueue!(x.w,carrier,delta,:projected,x.h,UInt32(x.p.cells.attempted_steps))
    end
    Interior.device_cellular_status!(x.p.cells);x.w.time+=n*dt
end
function fixture(n,founders)
    cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","native","living-patch-recruitment-16.toml"))
    cfg["grid"]["shape"]=fill(n,3);cfg["cellular"]["founders"]=founders
    recipe=tempname()*".toml";open(io->TOML.print(io,cfg),recipe,"w")
    f=N.fixture(recipe);h=N.sidecar(f.p;dt=.02,horizon=2.,act_time=100.,inheritance_fraction=.5)
    (;w=f.w,p=f.p,h)
end
function arrays(x)
    d=Dict{String,Any}()
    for (name,group) in pairs((cells=x.p.cells.state,environment=x.p.cells.environment,habitat=x.p.habitat.arrays,
            feed=x.p.habitat.feedstock.arrays,private=x.h.state,events=x.h.events))
        for key in keys(group);d["$name.$key"]=copy(group[key]);end
    end
    work=Interior.device_resource_workspace(x.p.cells);d["workspace.grants"]=copy(work.grants)
    for (i,level) in enumerate(work.levels),key in (:high,:low);d["workspace.level$i.$key"]=copy(level[key]);end
    d["workspace.highwater"]=copy(Interior.device_lifecycle_workspace(x.p.cells).highwater)
    d
end
function exact_compare(a,b)
    a.w.time==b.w.time && a.h.steps==b.h.steps && a.p.cells.attempted_steps==b.p.cells.attempted_steps||error("clock mismatch")
    aa=arrays(a);bb=arrays(b);count=0
    for key in keys(aa)
        isequal(aa[key],bb[key])||error("array mismatch $key");count+=length(aa[key])
    end
    count
end
function geometry_counts(x)
    s=x.p.cells.state;g=x.p.cells.geometry
    active=findall(==(1),s.active);sites=Set{Int32}();edges=0;contacts=0
    for i in active
        for k in 1:s.support_count[i];push!(sites,s.sites[i,k]);edges+=1;end
        for j in active
            j>i||continue
            d2=0f0
            for axis in 1:3
                delta=abs(s.position[i,axis]-s.position[j,axis]);delta=min(delta,g.extent[axis]-delta);d2+=delta^2
            end
            contacts+=d2<(s.radii[i]+s.radii[j])^2
        end
    end
    (;population=length(active),edges,unique_sites=length(sites),union_fraction=length(sites)/prod(Int.(g.dims)),
        actual_contacts=contacts,support_storage_bytes=sizeof(s.sites)+sizeof(s.support),actor_updates=sum(x.h.state.actor_updates),births=x.p.cells.environment.metadata[3])
end
function profile(out;resolution=64,founders_list=(24,46),steps=8)
    rows=Any[]
    for founders in founders_list
        println((;stage="fixture",resolution,founders));flush(stdout)
        original=fixture(resolution,founders);advance!(original,220)
        sum(original.h.state.actor_updates)>0||error("paid actor not exercised")
        disposable=deepcopy(original);advance!(disposable,1;timer=CPUKernelTimes())
        a=deepcopy(original);b=deepcopy(original);timer=CPUKernelTimes()
        GC.gc();plain=@timed advance!(a,steps)
        GC.gc();instrumented=@timed advance!(b,steps;timer)
        checked=exact_compare(a,b)
        labels=sort(collect(keys(timer.seconds));by=k->timer.seconds[k],rev=true)
        timings=[(;kernel=k,seconds=timer.seconds[k],calls=timer.calls[k],ms_per_step=1000timer.seconds[k]/steps) for k in labels]
        row=(;resolution,founders,steps,geometry=geometry_counts(a),plain_step_ms=1000plain.time/steps,
            instrumented_step_ms=1000instrumented.time/steps,profile_overhead_ratio=instrumented.time/plain.time,
            exact_array_elements=checked,timings)
        push!(rows,row);println(row);flush(stdout)
        open(io->JSON3.write(io,(;schema="living-support-profile-v1",threads=Threads.nthreads(),rows)),out,"w")
    end
end
end
