isdefined(@__MODULE__,:dc_resource_exchange!) || include(joinpath(@__DIR__,"device_resource_atomic_oracle.jl"))
# Opt-in exact-zero RMW elimination. No production methods are replaced.
# Zero increments use atomic LOADS, never racing ordinary reads, so finite checks
# remain in the same kernel even alongside threads with positive grants.
function resource_zero_definition(path,name)
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
    found===nothing && error("resource experiment source definition missing")
    found
end
function resource_zero_install!()
    path=joinpath(@__DIR__,"..","src","device_cellular_kernels.jl")
    bytes2hex(SHA.sha256(read(path)))==DEVICE_CELLULAR_SOURCE_PIN.kernels ||
        error("resource experiment source differs from loaded kernel pin")
    definition=resource_zero_definition(joinpath(@__DIR__,"device_resource_atomic_oracle.jl"),:dc_resource_exchange!)
    definition.args[1].args[1]=:dc_resource_zero_exchange!
    replaced=0
    for (index,expr) in pairs(definition.args[2].args)
        # Kernel body wraps the arithmetic in its failure-latch conditional.
        expr isa Expr && expr.head===:if || continue
        body=expr.args[2]
        for (j,line) in pairs(body.args)
            line isa Expr && line.head===:(=) && line.args[1] isa Expr && line.args[1].head===:tuple || continue
            names=line.args[1].args
            if names==[:oldsource,:newsource] || names==[:oldsink,:newsink]
                source=names[1]===:oldsource
                amount=source ? :supplied : :removed
                column=source ? 5 : 6
                old,new=names
                atomic_load=:(KernelAbstractions.@atomic e.resource_ledger[k,$column])
                body.args[j]=quote
                    if $amount==0f0
                        $old=$atomic_load
                        $new=$old
                    else
                        $line
                    end
                end
                replaced+=1
            end
        end
    end
    replaced==2 || error("resource experiment expected exactly two atomic increments")
    Core.eval(@__MODULE__,Expr(:macrocall,Symbol("@kernel"),LineNumberNode(0),definition))
end
resource_zero_install!()
const DEVICE_RESOURCE_ZERO_SOURCE=(study=bytes2hex(SHA.sha256(read(@__FILE__))),base=DEVICE_CELLULAR_SOURCE_PIN)
struct DeviceResourceZeroCarrier{P}
    process::P
end
@inline Base.getproperty(p::DeviceResourceZeroCarrier,key::Symbol)=key===:process ? getfield(p,key) : getproperty(getfield(p,:process),key)
function Base.setproperty!(p::DeviceResourceZeroCarrier,key::Symbol,value)
    key===:attempted_steps || throw(ArgumentError("only attempted_steps may be forwarded"))
    setproperty!(p.process,key,value)
end
device_cellular_validate(w,p::DeviceResourceZeroCarrier)=device_cellular_validate(w,p.process)
device_cellular_timestep_limit(p::DeviceResourceZeroCarrier)=device_cellular_timestep_limit(p.process)
device_cellular_status!(p::DeviceResourceZeroCarrier)=device_cellular_status!(p.process)
device_support_descriptor(p::DeviceResourceZeroCarrier)=device_support_descriptor(p.process)
function device_cellular_launch!(p::DeviceResourceZeroCarrier,kernel,args...;ndrange)
    replacement=kernel===dc_resource_exchange! ? dc_resource_zero_exchange! : kernel
    device_cellular_launch!(p.process,replacement,args...;ndrange)
end
function device_resource_zero_step!(w,p::DeviceCellular,dt;execution=:projected,anabolism=nothing,check_status=true)
    device_cellular_validate(w,p)
    execution in (:staged,:projected) || throw(ArgumentError("unknown execution recipe"))
    delta=Float32(dt)
    isfinite(delta)&&0<delta<=device_cellular_timestep_limit(p) || throw(ArgumentError("invalid timestep"))
    p.attempted_steps+=1
    device_cellular_enqueue!(w,DeviceResourceZeroCarrier(p),delta,execution,anabolism,UInt32(p.attempted_steps))
    check_status && device_cellular_status!(p)
    p
end
function device_resource_zero_step!(w,p::DeviceConstructiveEcology,dt;kwargs...)
    check_dt(w,p,dt)
    device_resource_zero_step!(w,p.cells,dt;anabolism=p.habitat,kwargs...)
    p
end
