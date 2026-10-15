# Compatibility carrier for the integrated lifecycle optimization.
# Historical source-derived versions remain in their immutable source archives.
# Current experiments delegate to the production shared lifecycle implementation.
isdefined(@__MODULE__,:DeviceLifecycleWorkspace) || error("this helper requires integrated lifecycle workspace; use the historical pinned archive for older runtimes")
const DEVICE_LIFECYCLE_PREFIX_SOURCE=(study=bytes2hex(SHA.sha256(read(@__FILE__))),base=DEVICE_CELLULAR_SOURCE_PIN,
    implementation="production shared dense/prefix lifecycle")
struct DeviceLifecyclePrefixCarrier{P,A}
    process::P
    highwater::A
end
function DeviceLifecyclePrefixCarrier(p)
    DeviceLifecyclePrefixCarrier(p,device_lifecycle_workspace(p).highwater)
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
    device_cellular_launch!(p.process,kernel,args...;ndrange)
end
function device_lifecycle_prefix_step!(w,p::DeviceLifecyclePrefixCarrier,dt;kwargs...)
    device_cellular_step!(w,p.process,dt;lifecycle=:prefix,kwargs...)
    p
end
