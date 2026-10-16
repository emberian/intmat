# Explicit F64 chemical storage, original shared phase/transport/reaction laws.
include(joinpath(@__DIR__,"..","studies","native_brine_producer.jl"))
using Serialization
c=TOML.parsefile(ARGS[1]);rows=Any[]
for n in (16,32)
 cpu=phase_runtime_fixture(c,n;chemical_type=Float64);gpu=producer_upload(cpu)
 initial=(phase=sum(Float64,cpu.p.phi),chemical=sum(cpu.a)+sum(cpu.b))
 for _ in 1:200
  coupled_runtime_step!(cpu,c,.01f0,true);coupled_runtime_step!(gpu,c,.01f0,true)
 end
 errors=(phi=maximum(abs,cpu.p.phi.-Array(gpu.p.phi)),a=maximum(abs,cpu.a.-Array(gpu.a)),b=maximum(abs,cpu.b.-Array(gpu.b)),mineral=maximum(abs,cpu.w.fields.deposit.-Array(gpu.w.fields.deposit)))
 maximum(errors)<2e-5||error("full coupled parity failed $errors")
 cs=producer_scalars(cpu,initial);gs=producer_scalars(gpu,initial)
 abs(cs.chemical_relative_error)<1e-12&&abs(gs.chemical_relative_error)<1e-12||error("F64 chemical mass gate failed")
 cp=Interior.implicit_phase_checkpoint(gpu.w,gpu.p);io=IOBuffer();serialize(io,cp);seekstart(io);w,p=Interior.implicit_phase_restore(deserialize(io);backend)
 eltype(w.fields.deposit)===Float64||error("restart narrowed chemical storage")
 maximum(abs,Array(w.fields.deposit).-Array(gpu.w.fields.deposit))==0||error("restart changed mineral")
 w.time==gpu.w.time||error("restart clock changed")
 push!(rows,(;n,steps=200,errors,cpu=cs,gpu=gs));println(JSON3.write(last(rows)));flush(stdout)
 cpu=nothing;gpu=nothing;w=nothing;p=nothing;cp=nothing;GC.gc()
end
open(io->JSON3.pretty(io,(;rows,chemical_storage="Float64",phase_storage="Float32",shared_law=true)),ARGS[2],"w")
