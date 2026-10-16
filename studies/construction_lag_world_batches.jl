isdefined(@__MODULE__,:SidecarBatchPolicy)||include("learning_world_batches.jl")

"""Register an exact held-command intervention around the shared joint policy.

Pulse timing is cohort-common; axis, amplitude, producer and memory are private
member values. Both saved commands and original source masks survive restart.
"""
function register_construction_lag_world_batches!(P;name=:construction_lag)
    joint=register_joint_construction_world_batches!(P.J)
    path=joinpath(@__DIR__,"native_construction_lag.jl")
    kernels=(joint.kernels...,(P.perturb_drive!,P,path))
    policy=SidecarBatchPolicy(name,P.PulseHabitat,P.L,kernels,joint.member_options,
        (:memory,:axis,:amplitude,:producer),(:origin_step,:duration_steps),(:original_masks,:saved_drive),
        (h,physical,state,backend,options,cache,properties,launcher)->begin
            base=joint.repack(h.base,physical,state,backend,options,cache,properties,launcher)
            P.PulseHabitat(base,properties.axis,properties.amplitude,properties.producer,cache.saved_drive,
                properties.origin_step,properties.duration_steps)
        end,
        (p,x)->begin
            base=joint.restore(p,x)
            P.PulseHabitat(base,x.properties.axis,x.properties.amplitude,x.properties.producer,zero(base.state.drive),
                x.properties.origin_step,x.properties.duration_steps)
        end,
        merge(joint.source,(pulse=bytes2hex(sha256(read(path))),)))
    register_sidecar_batch_policy!(policy)
end
