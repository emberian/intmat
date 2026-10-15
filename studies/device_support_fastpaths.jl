# Opt-in experiment; include in Interior after the device habitat implementation.
# Derive numerical bodies from the loaded source, changing only named call edges
# and exact-zero guards. No production methods or persistent state are replaced.
function support_fast_definition(path,name)
    function search(x)
        x isa Expr || return nothing
        if x.head===:function
            signature=x.args[1]
            while signature isa Expr && signature.head===:where;signature=signature.args[1];end
            signature isa Expr && signature.head===:call && signature.args[1]===name && return deepcopy(x)
        end
        for child in x.args
            found=search(child);found===nothing || return found
        end
        nothing
    end
    found=search(Meta.parseall(read(path,String)))
    found===nothing && error("support experiment cannot find source definition $name")
    found
end
function support_fast_rewrite(x,replacements)
    x isa Symbol && return get(replacements,x,x)
    x isa Expr || return x
    Expr(x.head,(support_fast_rewrite(a,replacements) for a in x.args)...)
end
function support_fast_install!()
    root=@__DIR__
    kernels=joinpath(root,"..","src","device_cellular_kernels.jl")
    habitat=joinpath(root,"..","src","device_habitat.jl")
    bytes2hex(SHA.sha256(read(kernels)))==DEVICE_CELLULAR_SOURCE_PIN.kernels &&
        bytes2hex(SHA.sha256(read(habitat)))==DEVICE_HABITAT_SOURCE_PIN.helper ||
        throw(ArgumentError("support experiment source files differ from loaded runtime pins"))
    waste=support_fast_definition(kernels,:dc_return_waste!)
    waste=support_fast_rewrite(waste,Dict(:dc_return_waste! => :dc_support_fast_return_waste!))
    guarded=filter(x->x isa Expr && x.head===:if,waste.args[2].args)
    length(guarded)==1 || error("support experiment waste guard source changed")
    guarded[1].args[1]=Expr(:&&,guarded[1].args[1],:(s.waste[i]!=0f0))
    Core.eval(@__MODULE__,Expr(:macrocall,Symbol("@kernel"),LineNumberNode(0),waste))
    exportfn=support_fast_definition(habitat,:dh_export!)
    exportfn=support_fast_rewrite(exportfn,Dict(:dh_export! => :dh_support_fast_export!))
    pushfirst!(exportfn.args[2].args,:(amount==0. && return 0.))
    Core.eval(@__MODULE__,exportfn)
    build=support_fast_definition(habitat,:device_anabolic_build!)
    build=support_fast_rewrite(build,Dict(:device_anabolic_build! => :device_support_fast_anabolic_build!,
        :dh_export! => :dh_support_fast_export!))
    Core.eval(@__MODULE__,build)
    metabolism=support_fast_definition(kernels,:dc_metabolism!)
    metabolism=support_fast_rewrite(metabolism,Dict(:dc_metabolism! => :dc_support_fast_metabolism!,
        :device_anabolic_build! => :device_support_fast_anabolic_build!))
    Core.eval(@__MODULE__,Expr(:macrocall,Symbol("@kernel"),LineNumberNode(0),metabolism))
    nothing
end
device_support_fast_anabolic_build!(args...)=device_anabolic_build!(args...)
support_fast_install!()
const DEVICE_SUPPORT_FASTPATH_SOURCE=(
    study=bytes2hex(SHA.sha256(read(@__FILE__))),base=DEVICE_CELLULAR_SOURCE_PIN,habitat=DEVICE_HABITAT_SOURCE_PIN)

struct DeviceSupportFastCarrier{P}
    process::P
end
@inline Base.getproperty(p::DeviceSupportFastCarrier,key::Symbol)=
    key===:process ? getfield(p,key) : getproperty(getfield(p,:process),key)
function Base.setproperty!(p::DeviceSupportFastCarrier,key::Symbol,value)
    key===:attempted_steps || throw(ArgumentError("only attempt bookkeeping is forwarded by support experiment"))
    setproperty!(p.process,key,value)
end
device_cellular_validate(w,p::DeviceSupportFastCarrier)=device_cellular_validate(w,p.process)
device_cellular_timestep_limit(p::DeviceSupportFastCarrier)=device_cellular_timestep_limit(p.process)
device_cellular_status!(p::DeviceSupportFastCarrier)=device_cellular_status!(p.process)
device_support_descriptor(p::DeviceSupportFastCarrier)=device_support_descriptor(p.process)
function device_cellular_launch!(p::DeviceSupportFastCarrier,kernel,args...;ndrange)
    replacement=kernel===dc_return_waste! ? dc_support_fast_return_waste! :
        kernel===dc_metabolism! ? dc_support_fast_metabolism! : kernel
    device_cellular_launch!(p.process,replacement,args...;ndrange)
end
function device_support_fast_step!(w,p::DeviceCellular,dt;execution=:projected,anabolism=nothing,check_status=true)
    device_cellular_validate(w,p)
    execution in (:staged,:projected) || throw(ArgumentError("unknown support experiment execution recipe"))
    delta=Float32(dt)
    isfinite(delta)&&0<delta<=device_cellular_timestep_limit(p) || throw(ArgumentError("invalid support experiment timestep"))
    p.attempted_steps+=1
    device_cellular_enqueue!(w,DeviceSupportFastCarrier(p),delta,execution,anabolism,UInt32(p.attempted_steps))
    check_status && device_cellular_status!(p)
    p
end
function device_support_fast_step!(w,p::DeviceConstructiveEcology,dt;kwargs...)
    check_dt(w,p,dt)
    device_support_fast_step!(w,p.cells,dt;anabolism=p.habitat,kwargs...)
    p
end
