# Experimental HIP graph owner. AMDGPU is passed by the configured probe runner;
# this file adds no dependency to Interior's portable CPU installation.
isdefined(@__MODULE__,:DeviceGraphStepIndex) || include("device_graph_clock.jl")
struct DeviceCellularGraph{W,P,H,G,E,R,O,HO}
    world::W
    process::P
    anabolism::H
    graph::G
    executable::E
    runtime::R
    delta::Float32
    execution::Symbol
    process_ownership::O
    habitat_ownership::HO
end
DeviceCellularGraph(w,p,h,graph,executable,runtime,delta,execution)=
    DeviceCellularGraph(w,p,h,graph,executable,runtime,delta,execution,p.ownership,h===nothing ? nothing : h.ownership)
function device_graph_preflight(w,p,anabolism)
    device_cellular_validate(w,p);device_cellular_status!(p)
    anabolism===nothing || validate(w,DeviceConstructiveEcology(p,anabolism))
    completed=Int(Array(p.environment.metadata)[6])
    completed==p.attempted_steps || throw(ArgumentError("graph requires a complete native state"))
    completed
end
function device_cellular_capture_graph(w,p,dt,runtime;execution=:projected,anabolism=nothing)
    nameof(typeof(w.backend))===:ROCBackend || throw(ArgumentError("HIP capture requires AMD ROCBackend"))
    execution in (:staged,:projected) || throw(ArgumentError("unknown graph execution recipe"))
    delta=Float32(dt)
    isfinite(delta)&&0<delta<=device_cellular_timestep_limit(p) || throw(ArgumentError("invalid graph timestep"))
    anabolism===nothing || delta<=stability_limit(w,DeviceConstructiveEcology(p,anabolism)) ||
        throw(ArgumentError("invalid coupled graph timestep"))
    completed=device_graph_preflight(w,p,anabolism)
    clock=DeviceGraphStepIndex(p.environment.metadata)
    # Caller warms exact specializations on a separate disposable state. Never
    # execute/retry a stateful body automatically when capture fails.
    graph=runtime.HIP.capture() do
        device_cellular_enqueue!(w,p,delta,execution,anabolism,clock)
    end
    executable=runtime.HIP.instantiate(graph)
    device_graph_preflight(w,p,anabolism)==completed || error("capture unexpectedly advanced state")
    DeviceCellularGraph(w,p,anabolism,graph,executable,runtime,delta,execution)
end
function device_cellular_graph_replay!(owner::DeviceCellularGraph,count::Integer=1;check_status=true)
    count>0 || throw(ArgumentError("positive replay count required"))
    p=owner.process
    p.ownership===owner.process_ownership || throw(ArgumentError("captured cellular buffers/ownership were replaced"))
    (owner.anabolism===nothing || owner.anabolism.ownership===owner.habitat_ownership) ||
        throw(ArgumentError("captured habitat buffers/ownership were replaced"))
    completed=device_graph_preflight(owner.world,p,owner.anabolism)
    count<=typemax(Int32)-completed || throw(ArgumentError("device complete-step counter would overflow"))
    for _ in 1:count
        # The ordinary API also records the attempt before its first launch.
        # A failed or latched replay cannot be checkpointed as a complete state.
        p.attempted_steps+=1
        owner.runtime.HIP.launch(owner.executable)
    end
    check_status && device_cellular_status!(p)
    p
end
