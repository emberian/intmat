#!/usr/bin/env julia
# Compatibility gate only. Run exclusively in an assigned AMD window.
include("device_cellular_setup.jl")
using Interior: KernelAbstractions
using Interior.KernelAbstractions: @kernel, @index
backend isa AMDGPU.ROCBackend || error("native AMD backend required")
isdefined(Interior,:DeviceGraphStepIndex) || Base.include(Interior,joinpath(@__DIR__,"..","studies","device_graph_clock.jl"))
length(ARGS)==1 || error("result JSON path required")
@kernel function graph_capture_rng!(out,status,step)
    i=@index(Global,Linear)
    if status[1]==0
        out[i]+=Interior.dc_random(UInt32(911),i,step,11)
    end
end
@kernel function graph_capture_complete!(metadata,status)
    i=@index(Global,Linear)
    if i==1 && status[1]==0
        metadata[6]+=Int32(1)
    end
end
out=KernelAbstractions.zeros(backend,Float32,256)
metadata=KernelAbstractions.zeros(backend,Int32,7)
status=KernelAbstractions.zeros(backend,Int32,1)
clock=Interior.DeviceGraphStepIndex(metadata)
function enqueue_graph_probe!(out,metadata,status,clock)
    graph_capture_rng!(backend,128)(out,status,clock;ndrange=256)
    graph_capture_complete!(backend,128)(metadata,status;ndrange=1)
    nothing
end
# Warm both exact specializations outside capture, then restore explicit state.
enqueue_graph_probe!(out,metadata,status,clock);AMDGPU.synchronize()
fill!(out,0);fill!(metadata,0);AMDGPU.synchronize()
started=time_ns()
graph=AMDGPU.HIP.capture() do
    enqueue_graph_probe!(out,metadata,status,clock)
end
exec=AMDGPU.HIP.instantiate(graph)
capture_seconds=(time_ns()-started)/1e9
all(iszero,Array(out)) && all(iszero,Array(metadata)) || error("capture executed device work unexpectedly")
for _ in 1:5
    AMDGPU.HIP.launch(exec)
end
AMDGPU.synchronize()
actual=Array(out);complete=Array(metadata)[6]
expected=zeros(Float32,256)
for step in UInt32(1):UInt32(5),i in 1:256
    expected[i]+=Interior.dc_random(UInt32(911),i,step,11)
end
actual==expected && complete==5 || error("captured RNG clock repeats or diverges")
fill!(status,8);AMDGPU.synchronize()
AMDGPU.HIP.launch(exec);AMDGPU.synchronize()
Array(out)==actual && Array(metadata)[6]==5 || error("failure latch did not suppress replay")
graph_path=joinpath(dirname(pathof(AMDGPU)),"hip","graph.jl")
result=(runtime=device_cellular_runtime(),capture_seconds,replays=5,nodes_per_replay=2,
    capture_did_not_execute=true,dynamic_rng_bitwise=true,failure_latch=true,
    note="Compatibility gate only; no complete-step acceleration claim.",
    graph_api_source_sha256=bytes2hex(Interior.SHA.sha256(read(graph_path))))
open(io->JSON3.pretty(io,result),ARGS[1],"w");println(JSON3.write(result))
