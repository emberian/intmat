"""Reconstructible lifecycle watermark; no physical inventory or checkpoint payload.

Only the lifecycle kernel writes this owned array. Direct caller edits of IDs,
active flags or new_parent require rebuilding the process ownership contract.
"""
struct DeviceLifecycleWorkspace{A,B,C}
    highwater::A
    backend::B
    capacity::Int
    contract::C
end
function DeviceLifecycleWorkspace(backend,capacity::Integer;initial=0,highwater=nothing)
    0<capacity<typemax(Int32) && 0<=initial<=capacity || throw(ArgumentError("invalid lifecycle watermark capacity/value"))
    if highwater===nothing
        highwater=KernelAbstractions.zeros(backend,Int32,1)
        fill!(highwater,Int32(initial))
    end
    size(highwater)==(1,) && eltype(highwater)===Int32 && KernelAbstractions.get_backend(highwater)==backend ||
        throw(ArgumentError("lifecycle watermark backend/shape/type differs"))
    0<=only(Array(highwater))<=capacity || throw(ArgumentError("lifecycle watermark exceeds capacity"))
    DeviceLifecycleWorkspace(highwater,backend,Int(capacity),(highwater=highwater,shape=size(highwater)))
end
"""Use already-host checkpoint arrays; IDs survive deaths and retain the watermark."""
function device_lifecycle_host_highwater(state)
    length(state.ids)==length(state.active)==length(state.new_parent) || throw(ArgumentError("lifecycle host roles differ"))
    something(findlast(i->state.ids[i]!=0 || state.active[i]!=0 || state.new_parent[i]!=0,eachindex(state.ids)),0)
end
function device_lifecycle_workspace_validate(work::DeviceLifecycleWorkspace,backend=work.backend,capacity=work.capacity)
    work.highwater===work.contract.highwater && size(work.highwater)==work.contract.shape==(1,) &&
        eltype(work.highwater)===Int32 && work.backend==backend && work.capacity==capacity &&
        KernelAbstractions.get_backend(work.highwater)==backend || throw(ArgumentError("lifecycle workspace role/backend/shape changed"))
    true
end
function device_lifecycle_workspace_nonalias(work,state,env,resource)
    arrays=(values(state)...,values(env)...,resource.grants,
        (a for level in resource.levels for a in (level.high,level.low))...)
    any(a->Base.mightalias(work.highwater,a),arrays) && throw(ArgumentError("lifecycle watermark aliases another process role"))
    true
end
device_lifecycle_workspace(p)=p.ownership.stages.lifecycle
