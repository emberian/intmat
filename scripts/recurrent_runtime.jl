# Import only the requested backend, before entering the simulation world age.
function prepare_recurrent_backend!(kind)
    if kind=="amd"
        Core.eval(@__MODULE__,:(using AMDGPU))
    elseif kind=="metal"
        Core.eval(@__MODULE__,:(using Metal))
    elseif kind!="cpu"
        throw(ArgumentError("backend must be cpu, metal or amd"))
    end
    nothing
end
function recurrent_backend(kind)
    if kind=="amd"
        AMDGPU.functional() || error("AMDGPU unavailable")
        AMDGPU.ROCBackend()
    elseif kind=="metal"
        Metal.functional() || error("Metal unavailable")
        Metal.MetalBackend()
    else
        CPU()
    end
end
function recurrent_capture_steps(steps,count)
    count>=0 || throw(ArgumentError("capture count must be nonnegative"))
    count==0 && return Set{Int}()
    Set(max.(1,round.(Int,steps.*collect(1:min(count,steps))./min(count,steps))))
end
