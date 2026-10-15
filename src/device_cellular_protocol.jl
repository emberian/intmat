# The generic process protocol is a synchronous observation/execution boundary.
# The native queued API retains its cheaper contract-only preflight.
function validate(w::World,p::DeviceCellular)
    device_cellular_validate(w,p)
    device_cellular_status!(p)
    nothing
end
stability_limit(w::World,p::DeviceCellular)=device_cellular_timestep_limit(p)
# Arena needs a separate observation/energy/portable-checkpoint adapter.
arena_process_supported(::DeviceCellular)=false

"""Strict inventory observation without downloading genomes or controller scratch.
Device resources are independent of World's nutrient/deposit inventories.
"""
function owned_material_mass(p::DeviceCellular)
    device_cellular_status!(p)
    s=p.state;e=p.environment
    active=findall(==(1),Array(s.active))
    mass=Array(s.mass);metabolites=Array(s.metabolites)
    sum(Float64,Array(e.fields))+sum(Float64,Array(e.inert))+sum(Float64,Array(e.inert_low))+
        sum(Float64,mass[active])+sum(Float64,metabolites[active,:])
end

"""Full floating-state finiteness audit by backend reductions, without array downloads.
This scans genomes and scratch at observation cadence; native compact monitoring
remains the cheaper choice when full diagnostics are not required.
"""
function process_state_finite(p::DeviceCellular)
    all(a->!(eltype(a)<:AbstractFloat)||all(isfinite,a),values(p.state)) &&
        all(a->!(eltype(a)<:AbstractFloat)||all(isfinite,a),values(p.environment))
end

process_ports(::DeviceCellular)=ProcessPorts(owned=(:state,:environment),
    drives="explicit finite sources and sinks in an independent three-species device pool",
    memory="device-resident development, genomes, bodies, contact graph and compensated ledgers",
    law="Float32 capacity-major cellular physiology; synchronous donor-limited junctions, mutual bonds and checked mechanical substeps")

"""Full detached observation with explicit array-role names; not a restart payload.
Use device_cellular_checkpoint/device_cellular_restore for source-pinned restart.
"""
function record_process_fields(p::DeviceCellular)
    x=device_cellular_snapshot(p)
    (; schema="interior-device-cellular-observation",schema_version=1,
       attempted_steps=x.attempted_steps,completed_steps=Int(x.environment.metadata[6]),
       cell_axis="first array dimension; persistent device slots; active==1 selects living cells",
       position_coordinates="physical periodic [0, extent); subtract extent/2 for centered rendering",
       resource_axes="linear voxel index in Julia column-major order, then species (two precursors and inert waste)",
       parameter_names=Tuple(keys(x.parameters)),parameters=Tuple(values(x.parameters)),
       geometry_names=Tuple(keys(x.geometry)),geometry=Tuple(values(x.geometry)),
       source_names=Tuple(keys(DEVICE_CELLULAR_SOURCE_PIN)),source=Tuple(values(DEVICE_CELLULAR_SOURCE_PIN)),
       (Symbol("state_",name)=>array for (name,array) in pairs(x.state))...,
       (Symbol("environment_",name)=>array for (name,array) in pairs(x.environment))...)
end
