module LivingResourceFusion
using Interior, KernelAbstractions, Statistics
using Interior: JSON3
include("living_support_worklists.jl")
const B=LivingSupportWorklists
const N=B.N

# Preserve the authoritative voxel helper and each first-level column's original
# strided addition order. Two independent supplied/removed accumulators share
# the voxel read; no exact-endpoint repair or changed reduction tree is applied.
@kernel function fused_first!(e,q,dt,high,low)
    block,k=@index(Global,NTuple)
    if e.status[1]==0
        count=size(high,1);sh=zero(eltype(high));sl=zero(eltype(low));rh=sh;rl=sl
        for j in 0:127
            i=block+j*count
            if i<=size(e.fields,1)
                supplied,removed=Interior.device_resource_voxel!(e,q,i,k,dt)
                sh,sl=Interior.device_resource_pair_add_checked(sh,sl,supplied,e.status)
                rh,rl=Interior.device_resource_pair_add_checked(rh,rl,removed,e.status)
            end
        end
        high[block,k]=sh;low[block,k]=sl
        high[block,k+3]=rh;low[block,k+3]=rl
    end
end

"""Remove only reconstructible grant backing; all physical roles are borrowed.
The compact workspace is owned by this explicit experimental carrier. Core
constructors still allocate V×6 first, so this is retained-data, not peak-init,
memory evidence. Its checkpoint must be wrapped with this study's source pin.
"""
function compact_resource!(p)
    old=Interior.device_resource_workspace(p);saved=sizeof(old.grants)
    grants=similar(old.grants,0,6)
    contract=(grants=grants,levels=old.levels,grant_shape=size(grants),level_shapes=map(x->(size(x.high),size(x.low)),old.levels))
    work=Interior.DeviceResourceWorkspace(grants,old.levels,contract)
    stages=merge(p.ownership.stages,(resource=work,));ownership=merge(p.ownership,(stages=stages,))
    typeof(ownership)===typeof(p.ownership)||error("compact ownership type changed")
    p.ownership=ownership
    Interior.device_resource_workspace_validate(work)
    saved
end
struct FusedCarrier{P}
    process::P
end
@inline Base.getproperty(p::FusedCarrier,k::Symbol)=k===:process ? getfield(p,k) : getproperty(getfield(p,:process),k)
Interior.device_support_descriptor(p::FusedCarrier)=Interior.device_support_descriptor(p.process)
function Interior.device_cellular_launch!(p::FusedCarrier,kernel,args...;ndrange)
    work=Interior.device_resource_workspace(p)
    if kernel===Interior.dc_resource_tree_grants!
        e,q,dt,grants=args
        grants===work.grants||error("unexpected grants workspace")
        level=first(work.levels)
        Interior.device_cellular_launch!(p.process,fused_first!,e,q,dt,level.high,level.low;ndrange=(size(level.high,1),3))
    elseif kernel===Interior.dc_resource_tree_reduce! && args[3]===work.grants && args[4]===nothing
        args[1]===first(work.levels).high && args[2]===first(work.levels).low||error("unexpected first reduction")
        nothing
    else
        Interior.device_cellular_launch!(p.process,kernel,args...;ndrange)
    end
end
function advance!(x,n;owner=nothing,dt=.02)
    for _ in 1:n
        Interior.device_cellular_validate(x.w,x.p.cells)
        delta=Float32(dt)
        isfinite(delta)&&0<delta<=Interior.device_cellular_timestep_limit(x.p.cells)||error("invalid timestep")
        x.p.cells.attempted_steps+=1
        carrier=owner===nothing ? x.p.cells : owner
        Interior.device_cellular_enqueue!(x.w,carrier,delta,:projected,x.h,UInt32(x.p.cells.attempted_steps))
    end
    Interior.device_cellular_status!(x.p.cells);x.w.time+=n*dt
end
function compare(a,b)
    a.w.time==b.w.time && a.h.steps==b.h.steps && a.p.cells.attempted_steps==b.p.cells.attempted_steps||error("clock mismatch")
    aa=B.arrays(a);bb=B.arrays(b);count=0
    for key in keys(aa)
        key=="workspace.grants" && continue # intentionally replaced by implicit transient grants
        isequal(aa[key],bb[key])||error("array mismatch $key");count+=length(aa[key])
    end
    count
end
function source_pin()
    (;study=bytes2hex(Interior.SHA.sha256(read(@__FILE__))),base=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat=Interior.DEVICE_HABITAT_SOURCE_PIN)
end
function benchmark(out;resolution=64,founders_list=(24,46),steps=12,repeats=3)
    rows=Any[]
    for founders in founders_list
        println((;stage="fixture",resolution,founders));flush(stdout)
        original=B.fixture(resolution,founders);advance!(original,220)
        sum(original.h.state.actor_updates)>0||error("paid actor not exercised")
        disposable=deepcopy(original);advance!(disposable,steps;owner=nothing)
        disposable=deepcopy(original);compact_resource!(disposable.p.cells)
        advance!(disposable,steps;owner=FusedCarrier(disposable.p.cells))
        samples=Any[];checked=0;last_b=nothing
        for repeat in 1:repeats
            a=deepcopy(original);b=deepcopy(original)
            saved=compact_resource!(b.p.cells);owner=FusedCarrier(b.p.cells)
            timings=Dict{Bool,Any}()
            for fused in (isodd(repeat) ? (false,true) : (true,false))
                x=fused ? b : a;GC.gc()
                t=@timed advance!(x,steps;owner=fused ? owner : nothing)
                t.compile_time==0||error("timing included $(t.compile_time) seconds compilation in fused=$fused")
                timings[fused]=(;seconds=t.time,gc_seconds=t.gctime,bytes=t.bytes)
            end
            checked+=compare(a,b);last_b=b
            push!(samples,(;baseline=timings[false],fused=timings[true],saved_grants_bytes=saved))
        end
        baseline=[s.baseline.seconds for s in samples];fused=[s.fused.seconds for s in samples]
        row=(;resolution,founders,steps,repeats,geometry=B.geometry_counts(last_b),samples,
            baseline_step_ms=1000median(baseline)/steps,fused_step_ms=1000median(fused)/steps,
            speedup=median(baseline)/median(fused),exact_array_elements=checked,
            saved_grants_bytes=first(samples).saved_grants_bytes,transient_semantics="V×6 grants replaced by implicit direct first-level contributions; all other arrays and clocks exact")
        push!(rows,row);println(row);flush(stdout)
        open(io->JSON3.write(io,(;schema="living-resource-fusion-v1",source=source_pin(),threads=Threads.nthreads(),rows)),out,"w")
    end
end
end
