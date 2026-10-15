# Experimental only: a device-resident RNG index for captured full-step queues.
# Include in Interior so existing lifecycle kernels specialize on this own type.
struct DeviceGraphStepIndex{A}
    metadata::A
end
@inline Base.:*(step::DeviceGraphStepIndex,factor::UInt32)=
    (UInt32(step.metadata[6])+UInt32(1))*factor
KernelAbstractions.Adapt.adapt_structure(to,step::DeviceGraphStepIndex)=
    DeviceGraphStepIndex(KernelAbstractions.Adapt.adapt(to,step.metadata))
