# Opt-in source-pinned experiment: exact endpoint pairs and fused first reduction.
# The underlying production process still owns its unused V×6 grants allocation.
# This study proves avoided traffic; it does not yet claim retained-memory savings.
@inline function resource_endpoint_difference(a::Float32,b::Float32)
    hi=a-b;part=hi-a
    hi,(a-(hi-part))+(-b-part)
end
@inline function resource_endpoint_gt(a::Float32,b::Float32,request::Float32)
    hi,lo=resource_endpoint_difference(a,b)
    hi==request ? lo>0f0 : hi>request ? hi-request > -lo : lo>request-hi
end
function resource_fused_install!()
    path=joinpath(@__DIR__,"..","src","device_resource_exchange.jl")
    source=read(path,String)
    bytes2hex(SHA.sha256(source))==DEVICE_CELLULAR_SOURCE_PIN.resource || error("resource source differs from loaded pin")
    matched=match(r"@inline function device_resource_voxel!\(.*?\nend"s,source)
    matched===nothing && error("authoritative voxel helper missing")
    definition=matched.match
    replacements=("device_resource_voxel!"=>"device_resource_endpoint_voxel!",
        "after-before>give"=>"resource_endpoint_gt(after,before,give)",
        "supplied=after-before"=>"supplied,supplied_low=resource_endpoint_difference(after,before)",
        "after-next>take"=>"resource_endpoint_gt(after,next,take)",
        "removed=after-next"=>"removed,removed_low=resource_endpoint_difference(after,next)",
        "!isfinite(supplied)||!isfinite(removed)"=>"!isfinite(supplied)||!isfinite(supplied_low)||!isfinite(removed)||!isfinite(removed_low)",
        "supplied,removed\n"=>"supplied,supplied_low,removed,removed_low\n")
    for (before,after) in replacements
        length(findall(before,definition))==1 || error("authoritative voxel helper changed at $before")
        definition=replace(definition,before=>after)
    end
    Core.eval(@__MODULE__,Meta.parse(definition))
end
resource_fused_install!()
const DEVICE_RESOURCE_FUSED_SOURCE=(study=bytes2hex(SHA.sha256(read(@__FILE__))),base=DEVICE_CELLULAR_SOURCE_PIN)

@kernel function dc_resource_fused_first!(e,q,dt,high,low)
    block,k=@index(Global,NTuple)
    if e.status[1]==0
        count=size(high,1);sh=0f0;sl=0f0;rh=0f0;rl=0f0
        for j in 0:127
            i=block+j*count
            if i<=size(e.fields,1)
                supplied,supplied_low,removed,removed_low=device_resource_endpoint_voxel!(e,q,i,k,dt)
                sh,sl=device_resource_pair_add_checked(sh,sl,supplied,e.status)
                if supplied_low!=0f0
                    sh,sl=device_resource_pair_add_checked(sh,sl,supplied_low,e.status)
                end
                rh,rl=device_resource_pair_add_checked(rh,rl,removed,e.status)
                if removed_low!=0f0
                    rh,rl=device_resource_pair_add_checked(rh,rl,removed_low,e.status)
                end
            end
        end
        high[block,k]=sh;low[block,k]=sl
        high[block,k+3]=rh;low[block,k+3]=rl
    end
end
struct DeviceResourceFusedCarrier{P}
    process::P
end
@inline Base.getproperty(p::DeviceResourceFusedCarrier,k::Symbol)=k===:process ? getfield(p,k) : getproperty(getfield(p,:process),k)
function Base.setproperty!(p::DeviceResourceFusedCarrier,k::Symbol,value)
    k===:attempted_steps || throw(ArgumentError("only attempted_steps may be forwarded"))
    setproperty!(p.process,k,value)
end
device_cellular_validate(w,p::DeviceResourceFusedCarrier)=device_cellular_validate(w,p.process)
device_cellular_timestep_limit(p::DeviceResourceFusedCarrier)=device_cellular_timestep_limit(p.process)
device_cellular_status!(p::DeviceResourceFusedCarrier)=device_cellular_status!(p.process)
device_support_descriptor(p::DeviceResourceFusedCarrier)=device_support_descriptor(p.process)
function device_cellular_launch!(p::DeviceResourceFusedCarrier,kernel,args...;ndrange)
    work=device_resource_workspace(p)
    if kernel===dc_resource_tree_grants!
        e,q,dt,grants=args
        grants===work.grants || error("unexpected grant workspace")
        level=first(work.levels)
        device_cellular_launch!(p.process,dc_resource_fused_first!,e,q,dt,level.high,level.low;ndrange=(size(level.high,1),3))
    elseif kernel===dc_resource_tree_reduce! && args[3]===work.grants && args[4]===nothing
        args[1]===first(work.levels).high && args[2]===first(work.levels).low || error("unexpected first reduction")
        # Already computed above. The authoritative shared phase stays unchanged.
        nothing
    else
        device_cellular_launch!(p.process,kernel,args...;ndrange)
    end
end
function device_resource_fused_step!(w,p::DeviceResourceFusedCarrier,dt;
        execution=:projected,anabolism=nothing,check_status=true,lifecycle=:prefix)
    device_cellular_validate(w,p)
    execution in (:staged,:projected) && lifecycle in (:prefix,:dense) || throw(ArgumentError("unknown execution recipe"))
    delta=Float32(dt)
    isfinite(delta)&&0<delta<=device_cellular_timestep_limit(p) || throw(ArgumentError("invalid timestep"))
    p.attempted_steps+=1
    device_cellular_enqueue!(w,p,delta,execution,anabolism,UInt32(p.attempted_steps),lifecycle)
    check_status && device_cellular_status!(p)
    p
end
# Experimental schedule identity must accompany its physical checkpoint.
function device_resource_fused_checkpoint(w,composition,p::DeviceResourceFusedCarrier)
    p.process===composition.cells || throw(ArgumentError("checkpoint carrier does not own these cells"))
    (format="Interior.ExperimentalResourceFused",version=1,source=DEVICE_RESOURCE_FUSED_SOURCE,
        checkpoint=device_habitat_checkpoint(w,composition))
end
function device_resource_fused_restore(w,x)
    x.format=="Interior.ExperimentalResourceFused" && x.version==1 && x.source==DEVICE_RESOURCE_FUSED_SOURCE ||
        throw(ArgumentError("fused experiment checkpoint/source differs"))
    composition=device_habitat_restore(w,x.checkpoint)
    composition,DeviceResourceFusedCarrier(composition.cells)
end
