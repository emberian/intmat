#!/usr/bin/env julia
# Numerical regression only: read a completed checkpoint's arrays; no restore,
# cell step, source transfer, trajectory continuation or writes to that run.
include("device_cellular_setup.jl")
using Serialization, Statistics, KernelAbstractions
length(ARGS)==3 || error("CHECKPOINT EXPECTED_SHA OUTPUT_JSON required")
path,expected,out=ARGS
actual=bytes2hex(Interior.SHA.sha256(read(path)));actual==expected || error("historical checkpoint hash differs")
x=deserialize(path);snap=x.cells.snapshot;g=snap.geometry;q=snap.parameters
steps=4000;dt=.001f0
results=Dict{String,Any}("checkpoint_sha256"=>actual,"checkpoint_source"=>x.cells.source,"snapshot_step"=>snap.attempted_steps,
    "purpose"=>"Repeated diffusion-only numerical fixture. Does not resume the ecological trajectory.","steps"=>steps,"dt"=>dt,
    "current_source"=>Interior.DEVICE_CELLULAR_SOURCE_PIN,"runtime"=>device_cellular_runtime())
for T in (Float32,Float64)
    fields=KernelAbstractions.zeros(backend,T,size(snap.environment.fields)...);copyto!(fields,T.(snap.environment.fields))
    scratch=similar(fields);permeability=KernelAbstractions.zeros(backend,Float32,length(snap.environment.permeability));copyto!(permeability,snap.environment.permeability)
    status=KernelAbstractions.zeros(backend,Int32,1);e=(;fields,scratch,permeability,status)
    before=vec(sum(Float64.(Array(fields));dims=1));V=size(fields,1)
    Interior.dc_diffuse!(backend,128)(e,q,g,dt;ndrange=(V,3));KernelAbstractions.synchronize(backend) # compile without committing
    started=time_ns()
    for _ in 1:steps
        Interior.dc_diffuse!(backend,128)(e,q,g,dt;ndrange=(V,3))
        Interior.dc_commit_resources!(backend,128)(e;ndrange=(V,3))
    end
    KernelAbstractions.synchronize(backend);seconds=(time_ns()-started)/1e9
    final=Array(fields);after=vec(sum(Float64.(final);dims=1));residual=after.-before
    results[string(T)]=(before,after,species_residual=residual,material_residual=sum(residual),potential_residual=sum(residual.*collect(q.potentials)),seconds,
        min=minimum(final),finite=all(isfinite,final))
    all(isfinite,final)&&minimum(final)>=0 || error("invalid diffusion result")
    open(io->JSON3.pretty(io,results),out,"w")
end
maximum(abs,results["Float64"].species_residual)<1e-8 || error("wide diffusion residual exceeded gate")
println(JSON3.write(results))
