# Historical atomic resource oracle; explicitly loaded by studies/tests only.
@kernel function dc_resource_exchange!(e,q,g,dt)
    i,k=@index(Global,NTuple)
    if e.status[1]==0
        before=e.fields[i,k];give=dt*q.source_rates[k]*e.source_masks[i,k]*e.quota[k,1]
        after=before+give
        after-before>give && (after=prevfloat(after))
        supplied=after-before
        take=min(after,dt*q.sink_rates[k]*e.sink_masks[i,k]*e.quota[k,2]);next=max(0f0,after-take)
        after-next>take && (next=nextfloat(next))
        removed=after-next;e.fields[i,k]=next
        oldsource,newsource=KernelAbstractions.@atomic e.resource_ledger[k,5] + supplied
        oldsink,newsink=KernelAbstractions.@atomic e.resource_ledger[k,6] + removed
        ((supplied>0 && newsource==oldsource) || (removed>0 && newsink==oldsink)) &&
            (KernelAbstractions.@atomic e.status[1]|=Int32(16))
        (!isfinite(newsource)||!isfinite(newsink)) && (KernelAbstractions.@atomic e.status[1]|=Int32(8))
        !isfinite(next) && (KernelAbstractions.@atomic e.status[1]|=Int32(8))
    end
end
@kernel function dc_resource_finish!(e)
    k=@index(Global,Linear)
    if e.status[1]==0
        dc_resource_add!(e,k,1,-e.resource_ledger[k,5])
        dc_resource_add!(e,k,2,-e.resource_ledger[k,6])
        dc_resource_add!(e,k,3,e.resource_ledger[k,5])
        dc_resource_add!(e,k,4,e.resource_ledger[k,6])
    end
end

# Oracle-only schedule adapter: ordinary production never includes this file.
struct DeviceAtomicResourceCarrier{P}
    process::P
end
@inline Base.getproperty(p::DeviceAtomicResourceCarrier,key::Symbol)=key===:process ? getfield(p,key) : getproperty(getfield(p,:process),key)
function Base.setproperty!(p::DeviceAtomicResourceCarrier,key::Symbol,value)
    key===:attempted_steps || throw(ArgumentError("only attempted_steps may be forwarded"))
    setproperty!(p.process,key,value)
end
device_cellular_validate(w,p::DeviceAtomicResourceCarrier)=device_cellular_validate(w,p.process)
device_cellular_timestep_limit(p::DeviceAtomicResourceCarrier)=device_cellular_timestep_limit(p.process)
device_cellular_status!(p::DeviceAtomicResourceCarrier)=device_cellular_status!(p.process)
device_support_descriptor(p::DeviceAtomicResourceCarrier)=device_support_descriptor(p.process)
function device_cellular_launch!(p::DeviceAtomicResourceCarrier,kernel,args...;ndrange)
    device_cellular_launch!(p.process,kernel,args...;ndrange)
end
function device_resource_phase!(p::DeviceAtomicResourceCarrier,e,q,g,dt)
    device_cellular_launch!(p,dc_resource_quota!,e,q,dt;ndrange=3)
    device_cellular_launch!(p,dc_resource_exchange!,e,q,g,dt;ndrange=size(e.fields))
    device_cellular_launch!(p,dc_resource_finish!,e;ndrange=3)
end
function device_atomic_resource_step!(w,p::DeviceAtomicResourceCarrier,dt;anabolism=nothing)
    device_cellular_validate(w,p)
    delta=Float32(dt)
    isfinite(delta)&&0<delta<=device_cellular_timestep_limit(p) || throw(ArgumentError("invalid timestep"))
    p.attempted_steps+=1
    device_cellular_enqueue!(w,p,delta,:projected,anabolism,UInt32(p.attempted_steps))
    device_cellular_status!(p)
end
