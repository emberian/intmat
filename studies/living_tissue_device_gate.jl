# Numerical/protocol gate for processes exposing the shared native tissue and
# reserve ports. Construct worlds/fixtures outside this helper; it adds no law.
using Interior,SHA
using Interior:JSON3,KernelAbstractions
function tissue_gate_valid(s,initial)
    s.finite && s.physical_min>=0 && s.neural_maxabs<=1.00001 &&
        abs(s.physical_material-initial.physical_material)<1e-5initial.physical_material &&
        (!hasproperty(s,:capacity_excess) || s.capacity_excess<=16eps(Float32)*s.capacity)
end
function tissue_device_gate(reference_world,reference,device_world,device,output;
        dt,steps=50,continuation_steps=11,measure,save_checkpoint,restore,
        backend_kind,provenance=NamedTuple())
    ispath(output) && error("fresh gate output required")
    backend=device_world.backend
    initial=measure(reference_world,reference)
    mkdir(output);cp(@__FILE__,joinpath(output,"device-gate-source.jl"))
    for _ in 1:steps
        step!(reference_world,reference,dt);step!(device_world,device,dt)
    end
    KernelAbstractions.synchronize(backend)
    a=Array(reference.tissue.state);b=Array(device.tissue.state)
    ra=Array(reference.reserve);rb=Array(device.reserve)
    state_error=maximum(abs.(a.-b));reserve_error=maximum(abs.(ra.-rb))
    cpu_stats=measure(reference_world,reference;state=a,reserve=ra)
    device_stats=measure(device_world,device;state=b,reserve=rb)
    numerical=state_error<1e-4 && reserve_error<1e-4 &&
        tissue_gate_valid(cpu_stats,initial) && tissue_gate_valid(device_stats,initial)
    record_state(reference_world,reference,joinpath(output,"gate-cpu");step=steps)
    record_state(device_world,device,joinpath(output,"gate-device");step=steps)
    checkpoint=joinpath(output,"gate-checkpoint.jls");save_checkpoint(device_world,device,checkpoint)
    restored_world,restored=restore(checkpoint;backend)
    exact_restore=Array(restored.tissue.state)==b && Array(restored.reserve)==rb && restored_world.time==device_world.time
    for _ in 1:continuation_steps
        step!(device_world,device,dt);step!(restored_world,restored,dt)
    end
    KernelAbstractions.synchronize(backend)
    exact_continuation=Array(restored.tissue.state)==Array(device.tissue.state) &&
        Array(restored.reserve)==Array(device.reserve) && restored_world.time==device_world.time
    continuation=measure(device_world,device)
    result=(passed=numerical&&exact_restore&&exact_continuation&&tissue_gate_valid(continuation,initial),
        backend=backend_kind,backend_type=string(typeof(backend)),dimensions=reference_world.grid.size,
        extent=reference_world.grid.extent,dt,steps,continuation_steps,
        state_maxabs_error=state_error,reserve_maxabs_error=reserve_error,
        exact_restore,exact_continuation,initial,reference=cpu_stats,device=device_stats,continuation,
        gate_source_sha256=bytes2hex(sha256(read(@__FILE__))),provenance,
        interpretation="Numerical parity, material/capacity budget, recording and exact same-backend continuation only. CPU plumbing is not an accelerator result; this does not establish morphology or biological behavior.")
    open(io->JSON3.pretty(io,result),joinpath(output,"gate.json"),"w")
    println(JSON3.write(result));flush(stdout)
    result
end
