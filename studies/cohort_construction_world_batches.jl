# Execution-only registration for decision-aligned return queues.
isdefined(@__MODULE__,:SidecarBatchPolicy)||include("learning_world_batches.jl")
function register_cohort_construction_world_batches!(Q;name=:cohort_construction)
    joint=register_joint_construction_world_batches!(Q.J)
    path=joinpath(@__DIR__,"native_cohort_construction_credit.jl")
    basekernels=Tuple(x for x in joint.kernels if x[1]!==Q.J.choose_joint!)
    kernels=(basekernels...,(Q.choose_cohort!,Q,path),(Q.mature!,Q,path),(Q.transaction_diagnostics!,Q,path))
    policy=SidecarBatchPolicy(name,Q.CohortHabitat,Q.L,kernels,(joint.member_options...,:critic_gain),
        joint.member_properties,joint.cohort_properties,joint.cache_names,
        (h,physical,state,backend,options,cache,properties,launcher)->Q.CohortHabitat(joint.repack(h.base,physical,state,backend,options,cache,properties,launcher)),
        (p,x)->Q.sidecar(p;memory=x.properties.memory,horizon=Float64(x.options.dt)*x.options.horizon_stride,
            act_time=Float64(x.options.dt)*x.options.act_stride,construction_sidecar_keywords(x.options)...),
        merge(joint.source,(cohort=bytes2hex(sha256(read(path))),)))
    register_sidecar_batch_policy!(policy)
end
