# Shared resident source/sink phase. Fixed compensated reduction replaces
# global per-voxel atomic credits. This changes Float32 summation order.
@inline function device_resource_voxel!(e,q,i,k,dt)
    before=e.fields[i,k];give=dt*q.source_rates[k]*e.source_masks[i,k]*e.quota[k,1]
    after=before+give
    after-before>give && (after=prevfloat(after))
    supplied=after-before
    take=min(after,dt*q.sink_rates[k]*e.sink_masks[i,k]*e.quota[k,2]);next=max(0f0,after-take)
    after-next>take && (next=nextfloat(next))
    removed=after-next;e.fields[i,k]=next
    (!isfinite(supplied)||!isfinite(removed)||!isfinite(next)) &&
        (KernelAbstractions.@atomic e.status[1]|=Int32(8))
    supplied,removed
end
@kernel function dc_resource_tree_grants!(e,q,dt,grants)
    i,k=@index(Global,NTuple)
    if e.status[1]==0
        supplied,removed=device_resource_voxel!(e,q,i,k,dt)
        grants[i,k]=supplied;grants[i,k+3]=removed
    end
end
# Split wide endpoint transfers before feeding the compensated Float32 tree.
@inline function device_resource_pair_add_checked(hi::Float32,lo::Float32,amount::Float64,status)
    high=Float32(amount);low=Float32(amount-Float64(high))
    hi,lo=device_resource_pair_add_checked(hi,lo,high,status)
    device_resource_pair_add_checked(hi,lo,low,status)
end
@inline function device_resource_pair_add_checked(hi,lo,amount,status)
    next,next_low=dc_pair_add(hi,lo,amount)
    if !(isfinite(next)&&isfinite(next_low))
        KernelAbstractions.@atomic status[1]|=Int32(8)
    elseif amount!=0f0 && next==hi && next_low==lo
        KernelAbstractions.@atomic status[1]|=Int32(16)
    end
    next,next_low
end
# Strided chunks coalesce adjacent lanes at each read; order is deterministic,
# unlike the old global per-voxel atomic sum. Carry each partial as a Float32 pair.
@kernel function dc_resource_tree_reduce!(high,low,input_high,input_low,status,input_count)
    block,column=@index(Global,NTuple)
    if status[1]==0
        count=size(high,1);hi=zero(eltype(high));lo=zero(eltype(low))
        for j in 0:127
            index=block+j*count
            if index<=input_count
                hi,lo=device_resource_pair_add_checked(hi,lo,input_high[index,column],status)
                if input_low!==nothing
                    hi,lo=device_resource_pair_add_checked(hi,lo,input_low[index,column],status)
                end
            end
        end
        high[block,column]=hi;low[block,column]=lo
        (!isfinite(hi)||!isfinite(lo)) && (KernelAbstractions.@atomic status[1]|=Int32(8))
    end
end
@inline function device_resource_budget_add!(ledger,ledger_low,status,k,column,dh,dl)
    before=ledger[k,column];before_low=ledger_low[k,column]
    hi,lo=dc_pair_add(before,before_low,dh)
    hi,lo=dc_pair_add(hi,lo,dl)
    if !(isfinite(hi)&&isfinite(lo))
        KernelAbstractions.@atomic status[1]|=Int32(8)
    elseif ((dh!=0f0||dl!=0f0)&&hi==before&&lo==before_low) || (column<=2&&(hi<0f0||(hi==0f0&&lo<0f0)))
        KernelAbstractions.@atomic status[1]|=Int32(16)
    else
        ledger[k,column]=hi;ledger_low[k,column]=lo
    end
end
@kernel function dc_resource_tree_finish!(e,high,low)
    k=@index(Global,Linear)
    if e.status[1]==0
        supplied=high[1,k];supplied_low=low[1,k]
        removed=high[1,k+3];removed_low=low[1,k+3]
        device_resource_budget_add!(e.resource_ledger,e.resource_ledger_low,e.status,k,1,-supplied,-supplied_low)
        device_resource_budget_add!(e.resource_ledger,e.resource_ledger_low,e.status,k,2,-removed,-removed_low)
        device_resource_budget_add!(e.resource_ledger,e.resource_ledger_low,e.status,k,3,supplied,supplied_low)
        device_resource_budget_add!(e.resource_ledger,e.resource_ledger_low,e.status,k,4,removed,removed_low)
        e.resource_ledger[k,5]=supplied;e.resource_ledger_low[k,5]=supplied_low
        e.resource_ledger[k,6]=removed;e.resource_ledger_low[k,6]=removed_low
    end
end

"""Reconstructible transient scratch, never physical inventory or checkpoint state."""
struct DeviceResourceWorkspace{G,L,C}
    grants::G
    levels::L
    contract::C
end
function DeviceResourceWorkspace(backend,n::Integer;grants=nothing,inventory_type=Float32)
    n>0 || throw(ArgumentError("positive resource workspace size required"))
    grants===nothing && (grants=KernelAbstractions.zeros(backend,inventory_type,n,6))
    size(grants)==(n,6) && eltype(grants)===inventory_type && KernelAbstractions.get_backend(grants)==backend ||
        throw(ArgumentError("resource grant scratch backend/shape/type differs"))
    levels=NamedTuple[]
    while true
        n=cld(n,128)
        push!(levels,(high=KernelAbstractions.zeros(backend,inventory_type,n,6),low=KernelAbstractions.zeros(backend,inventory_type,n,6)))
        n==1 && break
    end
    typed=Tuple(levels)
    contract=(grants=grants,levels=typed,grant_shape=size(grants),level_shapes=map(x->(size(x.high),size(x.low)),typed))
    DeviceResourceWorkspace(grants,typed,contract)
end
function device_resource_workspace_validate(work::DeviceResourceWorkspace)
    c=work.contract
    work.grants===c.grants && work.levels===c.levels && size(work.grants)==c.grant_shape &&
        map(x->(size(x.high),size(x.low)),work.levels)==c.level_shapes ||
        throw(ArgumentError("device resource workspace arrays or shapes changed"))
    true
end
function device_resource_workspace_nonalias(work,state,env;allowed_transients=())
    all(k->k in (:scratch,:demand,:scale,:returns,:occupancy),allowed_transients) ||
        throw(ArgumentError("only documented post-resource transients may borrow grant storage"))
    arrays=(work.grants,(a for level in work.levels for a in (level.high,level.low))...)
    for i in eachindex(arrays)
        for j in 1:i-1
            Base.mightalias(arrays[i],arrays[j]) && throw(ArgumentError("resource scratch roles alias"))
        end
        for physical in values(state)
            Base.mightalias(arrays[i],physical) && throw(ArgumentError("resource scratch aliases persistent state"))
        end
        for (role,physical) in pairs(env)
            if Base.mightalias(arrays[i],physical)
                i==1 && role in allowed_transients || throw(ArgumentError("resource scratch aliases a forbidden process role: $role"))
            end
        end
    end
    true
end
device_resource_workspace(p)=p.ownership.stages.resource
function device_resource_phase!(p,e,q,g,dt)
    work=device_resource_workspace(p);V=size(e.fields,1)
    device_cellular_launch!(p,dc_resource_quota!,e,q,dt;ndrange=3)
    device_cellular_launch!(p,dc_resource_tree_grants!,e,q,dt,work.grants;ndrange=(V,3))
    input_high=work.grants;input_low=nothing
    for level in work.levels
        device_cellular_launch!(p,dc_resource_tree_reduce!,level.high,level.low,input_high,input_low,e.status,size(input_high,1);ndrange=size(level.high))
        input_high=level.high;input_low=level.low
    end
    device_cellular_launch!(p,dc_resource_tree_finish!,e,input_high,input_low;ndrange=3)
    nothing
end
