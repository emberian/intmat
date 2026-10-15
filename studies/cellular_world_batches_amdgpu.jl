# Optional native-device slice adapter. Load only in an authorized AMD window.
using AMDGPU
isdefined(Main,:CellularWorldBatch) || include("cellular_world_batches.jl")
# Preserve AMDGPU's ROCDeviceArray indexing and native floating atomic dispatch.
# Construct a borrowed descriptor; no allocation, upload, or radius change.
@inline function batch_world(a::WorldBatchArray{T,N,A},w) where {T,N,A<:AMDGPU.ROCDeviceArray}
    AMDGPU.ROCDeviceArray(a.shape,pointer(a.storage,(Int(w)-1)*length(a)+1))
end
