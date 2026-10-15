# Experimental enqueue observer. It changes no device kernel or physical state.
# Include in Interior; the caller owns an exclusive GPU window and output.
struct DeviceProfileCarrier{P,O}
    process::P
    observer::O
end
@inline function Base.getproperty(p::DeviceProfileCarrier,key::Symbol)
    key===:process || key===:observer ? getfield(p,key) : getproperty(getfield(p,:process),key)
end
device_support_descriptor(p::DeviceProfileCarrier)=device_support_descriptor(getfield(p,:process))
mutable struct DeviceProfileTrace
    labels::Vector{Symbol}
end
mutable struct DeviceProfileEvents{E,R}
    events::E
    labels::Vector{Symbol}
    cursor::Int
    runtime::R
    completed::Bool
end
function device_profile_launch!(trace::DeviceProfileTrace,p,kernel,args...;ndrange)
    push!(trace.labels,nameof(kernel))
    nothing
end
function device_profile_launch!(observer::DeviceProfileEvents,p,kernel,args...;ndrange)
    cursor=observer.cursor+1
    cursor<=length(observer.events) || throw(ArgumentError("profile event pool exhausted"))
    observer.labels[cursor]===nameof(kernel) || throw(ArgumentError("profile launch schedule differs from plan"))
    observer.cursor=cursor
    observer.completed=false
    before,after=observer.events[cursor]
    observer.runtime.HIP.record(before)
    device_cellular_launch!(p,kernel,args...;ndrange)
    observer.runtime.HIP.record(after)
    nothing
end
function device_cellular_launch!(p::DeviceProfileCarrier,kernel,args...;ndrange)
    device_profile_launch!(getfield(p,:observer),getfield(p,:process),kernel,args...;ndrange)
end
function device_profile_plan(w,p,dt;execution=:projected,anabolism=nothing)
    device_cellular_validate(w,p)
    execution in (:staged,:projected) || throw(ArgumentError("invalid profile execution recipe"))
    delta=Float32(dt)
    isfinite(delta)&&0<delta<=device_cellular_timestep_limit(p) || throw(ArgumentError("invalid profile timestep"))
    trace=DeviceProfileTrace(Symbol[])
    carrier=DeviceProfileCarrier(p,trace)
    # Current hooks only enqueue kernels. This trace performs no launches,
    # allocates no device state, and advances no host or device clock.
    device_cellular_enqueue!(w,carrier,delta,execution,anabolism,UInt32(p.attempted_steps+1))
    trace.labels
end
function device_profile_event_pool(runtime,count::Integer)
    count>0 || throw(ArgumentError("positive profile event count required"))
    [(runtime.HIP.HIPEvent(runtime.stream();do_record=false,timing=true),
      runtime.HIP.HIPEvent(runtime.stream();do_record=false,timing=true)) for _ in 1:count]
end
function DeviceProfileEvents(plan,events,runtime;steps=1)
    steps>0 && length(events)==steps*length(plan) || throw(ArgumentError("event pool size differs from planned launches"))
    DeviceProfileEvents(events,repeat(plan,steps),0,runtime,false)
end
function device_profile_step!(w,p,observer,dt;execution=:projected,anabolism=nothing)
    device_cellular_validate(w,p)
    execution in (:staged,:projected) || throw(ArgumentError("invalid profile execution recipe"))
    delta=Float32(dt)
    isfinite(delta)&&0<delta<=device_cellular_timestep_limit(p) || throw(ArgumentError("invalid profile timestep"))
    p.attempted_steps+=1
    device_cellular_enqueue!(w,DeviceProfileCarrier(p,observer),delta,execution,anabolism,UInt32(p.attempted_steps))
    p
end
function device_profile_category(kernel::Symbol)
    kernel===:dc_support! && return :support_construction
    kernel===:dc_sense! && return :support_sensing_gradient_gather
    kernel===:dc_regulate! && return :support_uptake_scatter
    kernel===:dc_gather! && return :support_uptake_gather
    kernel in (:dc_return_waste!,:dc_porosity_scatter!) && return :support_waste_porosity_scatter
    kernel===:dc_diffuse! && return :dense_diffusion
    kernel in (:dc_resource_quota!,:dc_resource_exchange!,:dc_resource_finish!,:dc_resource_tree_grants!,:dc_resource_tree_reduce!,:dc_resource_tree_finish!,:dc_commit_resources!,
        :dc_allocate_sites!,:dc_clear_grid!,:dc_commit_waste!,:dc_clear_occupancy!,:dc_porosity_commit!) && return :dense_resource_passes
    kernel in (:dc_clear_heads!,:dc_geometry!) && return :bins_geometry
    kernel in (:dc_forces!,:dc_move!,:dc_commit_positions!) && return :mechanics
    kernel in (:dc_propose_bonds!,:dc_commit_bonds!) && return :contacts
    kernel in (:dc_lifecycle!,:dc_divide!) && return :lifecycle
    kernel===:dc_controller! && return :controller
    kernel in (:dc_exchange_scale!,:dc_exchange_commit!,:dc_commit_metabolites!) && return :junctions
    kernel===:dc_metabolism! && return :metabolism
    kernel===:dh_prepare! && return :habitat_support_readout_gather
    kernel===:dh_react_sorb! && return :dense_habitat_reaction_sorption
    kernel===:dh_permeability! && return :dense_habitat_permeability
    kernel===:dh_inherit! && return :habitat_lifecycle
    kernel in (:dc_validate_complete!,:dh_validate_complete!) && return :validation
    kernel===:dc_complete_step! && return :clock
    :other
end
function device_profile_wait!(observer::DeviceProfileEvents)
    observer.cursor==length(observer.events) || throw(ArgumentError("profile did not consume the complete event pool"))
    # Call inside the wall-time boundary, before allocating report summaries.
    if !observer.completed
        observer.runtime.HIP.synchronize(last(observer.events)[2])
        observer.completed=true
    end
    nothing
end
function device_profile_finish(observer::DeviceProfileEvents)
    device_profile_wait!(observer)
    samples=[(index=i,kernel=observer.labels[i],category=device_profile_category(observer.labels[i]),
        seconds=Float64(observer.runtime.HIP.elapsed(pair[1],pair[2]))) for (i,pair) in enumerate(observer.events)]
    kernels=Dict{Symbol,Tuple{Int,Float64}}();categories=Dict{Symbol,Tuple{Int,Float64}}()
    for sample in samples
        for (table,key) in ((kernels,sample.kernel),(categories,sample.category))
            count,seconds=get(table,key,(0,0.));table[key]=(count+1,seconds+sample.seconds)
        end
    end
    (;samples,kernels,categories,
      stream_seconds=Float64(observer.runtime.HIP.elapsed(first(observer.events)[1],last(observer.events)[2])),
      summed_event_seconds=sum(s->s.seconds,samples;init=0.),
      caveat="Intrusive event intervals include event-record/driver gaps. Compare with an uninstrumented matched batch; do not treat them as isolated instruction time.")
end
