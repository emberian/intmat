# Opt-in source-derived lifecycle experiment. Include into Interior.
# No production definitions are replaced. The high-water buffer is study-owned.
function lifecycle_prefix_definition(path,name)
    function search(x)
        x isa Expr || return nothing
        if x.head===:function && x.args[1] isa Expr && x.args[1].head===:call && x.args[1].args[1]===name
            return deepcopy(x)
        end
        for a in x.args
            result=search(a);result===nothing || return result
        end
        nothing
    end
    result=search(Meta.parseall(read(path,String)))
    result===nothing && error("missing lifecycle source definition")
    result
end
function lifecycle_prefix_install!()
    path=joinpath(@__DIR__,"..","src","device_cellular_kernels.jl")
    bytes2hex(SHA.sha256(read(path)))==DEVICE_CELLULAR_SOURCE_PIN.kernels || error("lifecycle source differs from loaded pin")
    fn=lifecycle_prefix_definition(path,:dc_lifecycle!)
    fn.args[1].args[1]=:dc_lifecycle_prefix!;push!(fn.args[1].args,:highwater)
    scans=0;insertions=0;updates=0
    function rewrite!(x)
        x isa Expr || return
        if x.head===:for && x.args[1] isa Expr && x.args[1].head===:(=) &&
                x.args[1].args[1]===:i && x.args[1].args[2]==:(1:Int(g.capacity))
            x.args[1].args[2]=:(1:Int(highwater[1]));scans+=1
        end
        if x.head===:block
            at=findfirst(a->a isa Expr && a.head===:(=) && a.args[1]===:births,x.args)
            if at!==nothing
                insert!(x.args,at,:(for i in Int(highwater[1])+1:min(Int(g.capacity),Int(highwater[1])+count)
                    freecount+=1;e.free[freecount]=i
                end));insertions+=1
            end
            committed=findfirst(a->a isa Expr && a.head===:+= && a.args[1]==:(e.metadata[3]) && a.args[2]===:births,x.args)
            if committed!==nothing
                insert!(x.args,committed+1,:(if births>0
                    highwater[1]=max(highwater[1],e.free[births])
                end));updates+=1
            end
        end
        for a in x.args;rewrite!(a);end
    end
    rewrite!(fn.args[2]);scans==insertions==updates==1 || error("lifecycle source schedule changed")
    Core.eval(@__MODULE__,Expr(:macrocall,Symbol("@kernel"),LineNumberNode(0),fn))
end
lifecycle_prefix_install!()
const DEVICE_LIFECYCLE_PREFIX_SOURCE=(study=bytes2hex(SHA.sha256(read(@__FILE__))),base=DEVICE_CELLULAR_SOURCE_PIN)
struct DeviceLifecyclePrefixCarrier{P,A}
    process::P
    highwater::A
end
function DeviceLifecyclePrefixCarrier(p)
    # One explicit setup transfer, not a per-step observation. IDs survive deaths;
    # include new_parent to make reconstruction conservative for partial fixtures.
    ids=Array(p.state.ids);active=Array(p.state.active);parents=Array(p.state.new_parent)
    h=findlast(i->ids[i]!=0 || active[i]!=0 || parents[i]!=0,eachindex(ids))
    highwater=KernelAbstractions.zeros(p.backend,Int32,1)
    copyto!(highwater,Int32[something(h,0)])
    DeviceLifecyclePrefixCarrier(p,highwater)
end
@inline Base.getproperty(p::DeviceLifecyclePrefixCarrier,k::Symbol)=
    k in (:process,:highwater) ? getfield(p,k) : getproperty(getfield(p,:process),k)
function Base.setproperty!(p::DeviceLifecyclePrefixCarrier,k::Symbol,v)
    k===:attempted_steps || throw(ArgumentError("only attempt bookkeeping may be forwarded"))
    setproperty!(p.process,k,v)
end
device_support_descriptor(p::DeviceLifecyclePrefixCarrier)=device_support_descriptor(p.process)
device_cellular_validate(w,p::DeviceLifecyclePrefixCarrier)=device_cellular_validate(w,p.process)
device_cellular_timestep_limit(p::DeviceLifecyclePrefixCarrier)=device_cellular_timestep_limit(p.process)
device_cellular_status!(p::DeviceLifecyclePrefixCarrier)=device_cellular_status!(p.process)
function device_cellular_launch!(p::DeviceLifecyclePrefixCarrier,kernel,args...;ndrange)
    if kernel===dc_lifecycle!
        device_cellular_launch!(p.process,dc_lifecycle_prefix!,args...,p.highwater;ndrange)
    else
        device_cellular_launch!(p.process,kernel,args...;ndrange)
    end
end
function device_lifecycle_prefix_step!(w,p::DeviceLifecyclePrefixCarrier,dt;anabolism=nothing,execution=:projected,check_status=true)
    device_cellular_validate(w,p)
    execution in (:staged,:projected) || throw(ArgumentError("invalid prefix execution"))
    delta=Float32(dt)
    isfinite(delta)&&0<delta<=device_cellular_timestep_limit(p) || throw(ArgumentError("invalid prefix timestep"))
    p.attempted_steps+=1
    device_cellular_enqueue!(w,p,delta,execution,anabolism,UInt32(p.attempted_steps))
    check_status && device_cellular_status!(p)
    p
end
