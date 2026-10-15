isdefined(@__MODULE__,:dc_resource_exchange!) || include(joinpath(@__DIR__,"device_resource_atomic_oracle.jl"))
# Opt-in deterministic compensated resource reduction; production stays unchanged.
# Physical voxel arithmetic is extracted from the pinned exchange kernel.
function resource_reduce_definition(path,name)
    function find(x)
        x isa Expr || return nothing
        if x.head===:function
            sig=x.args[1]
            sig isa Expr && sig.head===:call && sig.args[1]===name && return deepcopy(x)
        end
        for arg in x.args
            found=find(arg);found===nothing || return found
        end
        nothing
    end
    found=find(Meta.parseall(read(path,String)))
    found===nothing && error("resource reduction source definition missing")
    found
end
function resource_reduce_install!()
    path=joinpath(@__DIR__,"..","src","device_cellular_kernels.jl")
    bytes2hex(SHA.sha256(read(path)))==DEVICE_CELLULAR_SOURCE_PIN.kernels || error("resource reduction source pin mismatch")
    definition=resource_reduce_definition(joinpath(@__DIR__,"device_resource_atomic_oracle.jl"),:dc_resource_exchange!)
    definition.args[1].args[1]=:dc_resource_grants!
    push!(definition.args[1].args,:grants)
    guards=filter(x->x isa Expr && x.head===:if,definition.args[2].args)
    length(guards)==1 || error("exchange latch source changed")
    body=guards[1].args[2]
    firstatomic=findfirst(x->x isa Expr && x.head===:(=) && x.args[1] isa Expr && x.args[1].head===:tuple,body.args)
    firstatomic===nothing && error("exchange accumulation source changed")
    resize!(body.args,firstatomic-1)
    append!(body.args,(quote
        grants[i,k]=supplied;grants[i,k+3]=removed
        (!isfinite(supplied)||!isfinite(removed)||!isfinite(next)) &&
            (KernelAbstractions.@atomic e.status[1]|=Int32(8))
    end).args)
    Core.eval(@__MODULE__,Expr(:macrocall,Symbol("@kernel"),LineNumberNode(0),definition))
end
resource_reduce_install!()
const DEVICE_RESOURCE_REDUCTION_SOURCE=(study=bytes2hex(SHA.sha256(read(@__FILE__))),base=DEVICE_CELLULAR_SOURCE_PIN)

@inline function resource_reduce_add_checked(hi,lo,amount,status)
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
@kernel function dc_resource_reduce!(high,low,input_high,input_low,status,input_count)
    block,column=@index(Global,NTuple)
    if status[1]==0
        count=size(high,1);hi=0f0;lo=0f0
        for j in 0:127
            index=block+j*count
            if index<=input_count
                hi,lo=resource_reduce_add_checked(hi,lo,input_high[index,column],status)
                if input_low!==nothing
                    hi,lo=resource_reduce_add_checked(hi,lo,input_low[index,column],status)
                end
            end
        end
        high[block,column]=hi;low[block,column]=lo
        (!isfinite(hi)||!isfinite(lo)) && (KernelAbstractions.@atomic status[1]|=Int32(8))
    end
end
@inline function resource_reduced_add!(ledger,ledger_low,status,k,column,dh,dl)
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
@kernel function dc_resource_reduced_finish!(e,high,low)
    k=@index(Global,Linear)
    if e.status[1]==0
        supplied=high[1,k];supplied_low=low[1,k]
        removed=high[1,k+3];removed_low=low[1,k+3]
        resource_reduced_add!(e.resource_ledger,e.resource_ledger_low,e.status,k,1,-supplied,-supplied_low)
        resource_reduced_add!(e.resource_ledger,e.resource_ledger_low,e.status,k,2,-removed,-removed_low)
        resource_reduced_add!(e.resource_ledger,e.resource_ledger_low,e.status,k,3,supplied,supplied_low)
        resource_reduced_add!(e.resource_ledger,e.resource_ledger_low,e.status,k,4,removed,removed_low)
        e.resource_ledger[k,5]=supplied;e.resource_ledger_low[k,5]=supplied_low
        e.resource_ledger[k,6]=removed;e.resource_ledger_low[k,6]=removed_low
    end
end
struct DeviceResourceReductionCarrier{P,G,L}
    process::P
    grants::G
    levels::L
end
function DeviceResourceReductionCarrier(p)
    n=size(p.environment.fields,1)
    grants=KernelAbstractions.zeros(p.backend,Float32,n,6)
    levels=NamedTuple[]
    while true
        n=cld(n,128)
        push!(levels,(high=KernelAbstractions.zeros(p.backend,Float32,n,6),low=KernelAbstractions.zeros(p.backend,Float32,n,6)))
        n==1 && break
    end
    DeviceResourceReductionCarrier(p,grants,Tuple(levels))
end
@inline Base.getproperty(p::DeviceResourceReductionCarrier,key::Symbol)=key in (:process,:grants,:levels) ? getfield(p,key) : getproperty(getfield(p,:process),key)
function Base.setproperty!(p::DeviceResourceReductionCarrier,key::Symbol,value)
    key===:attempted_steps || throw(ArgumentError("only attempted_steps may be forwarded"))
    setproperty!(p.process,key,value)
end
device_cellular_validate(w,p::DeviceResourceReductionCarrier)=device_cellular_validate(w,p.process)
device_cellular_timestep_limit(p::DeviceResourceReductionCarrier)=device_cellular_timestep_limit(p.process)
device_cellular_status!(p::DeviceResourceReductionCarrier)=device_cellular_status!(p.process)
device_support_descriptor(p::DeviceResourceReductionCarrier)=device_support_descriptor(p.process)
function device_cellular_launch!(p::DeviceResourceReductionCarrier,kernel,args...;ndrange)
    if kernel===dc_resource_exchange!
        device_cellular_launch!(p.process,dc_resource_grants!,args...,p.grants;ndrange)
    elseif kernel===dc_resource_finish!
        e=only(args);input_high=p.grants;input_low=nothing
        for level in p.levels
            device_cellular_launch!(p.process,dc_resource_reduce!,level.high,level.low,input_high,input_low,e.status,size(input_high,1);ndrange=size(level.high))
            input_high=level.high;input_low=level.low
        end
        device_cellular_launch!(p.process,dc_resource_reduced_finish!,e,input_high,input_low;ndrange)
    else
        device_cellular_launch!(p.process,kernel,args...;ndrange)
    end
end
function device_resource_reduction_step!(w,p::DeviceResourceReductionCarrier,dt;execution=:projected,anabolism=nothing,check_status=true)
    device_cellular_validate(w,p)
    execution in (:staged,:projected) || throw(ArgumentError("unknown execution recipe"))
    delta=Float32(dt)
    isfinite(delta)&&0<delta<=device_cellular_timestep_limit(p) || throw(ArgumentError("invalid timestep"))
    p.attempted_steps+=1
    device_cellular_enqueue!(w,p,delta,execution,anabolism,UInt32(p.attempted_steps))
    check_status && device_cellular_status!(p)
    p
end
