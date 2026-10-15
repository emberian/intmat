#!/usr/bin/env julia
# Full device pipeline benchmark; includes bounded status/observation transfers.
include("device_cellular_setup.jl")
include("device_cellular_population.jl")
mechanical_substeps=parse(Int,get(ENV,"DEVICE_MECHANICAL_SUBSTEPS","16"))
side=parse(Int,get(ENV,"DEVICE_CELLS_SIDE","8"));steps=parse(Int,get(ENV,"DEVICE_STEPS","100"));batch=parse(Int,get(ENV,"DEVICE_BATCH","8"))
min(side,steps,batch)>0 || error("positive side, steps, batch required")
out=length(ARGS)>0 ? ARGS[1] : "outputs/device-cellular-benchmark"
ispath(out) && error("fresh output required")
# Compile the identical parameter/array types on a tiny separate living population.
warmw,warmp=device_cellular_population(2,backend;mechanical_substeps);Interior.device_cellular_step!(warmw,warmp,.01)
Interior.device_cellular_statistics(warmp);warmw=nothing;warmp=nothing;GC.gc()
w,p=device_cellular_population(side,backend;mechanical_substeps)
initial=Interior.device_cellular_audit(Interior.device_cellular_snapshot(p));initial_stats=Interior.device_cellular_statistics(p)
function timed_run!(w,p,steps,batch)
    observations=Any[];transfer_bytes=0;elapsed=0.
    allocation_bytes=@allocated begin
        elapsed=@elapsed begin
            for step in 1:steps
                Interior.device_cellular_step!(w,p,.01;check_status=false)
                if step%batch==0 || step==steps
                    push!(observations,Interior.device_cellular_statistics(p));transfer_bytes+=last(observations).observation_bytes
                end
            end
        end
    end
    observations,transfer_bytes,elapsed,allocation_bytes
end
observations,transfer_bytes,elapsed,allocation_bytes=timed_run!(w,p,steps,batch)
final=Interior.device_cellular_audit(Interior.device_cellular_snapshot(p));stats=last(observations)
material_error=final.material-initial.material-final.supplied+initial.supplied+final.removed-initial.removed
potential_error=final.potential-initial.potential-final.supplied_potential+initial.supplied_potential+final.removed_potential-initial.removed_potential+final.work-initial.work+final.heat-initial.heat
resident=sum(sizeof(a) for a in values(p.state))+sum(sizeof(a) for a in values(p.environment))
result=(runtime=device_cellular_runtime(),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,backend=kind,amd_optlevel=get(ENV,"DEVICE_AMD_OPTLEVEL","2"),law="experimental fixed-substep donor-limited DeviceCellular",side=side,capacity=Int(p.geometry.capacity),grid=p.geometry.dims,steps=steps,dt=.01,batch=batch,mechanical_substeps=mechanical_substeps,host_allocation_bytes=allocation_bytes,
    seconds=elapsed,full_steps_per_second=steps/elapsed,observation_bytes=transfer_bytes,process_resident_bytes=resident,
    initial=initial,initial_stats=initial_stats,final=final,final_stats=stats,material_error=material_error,potential_error=potential_error,
    timing_boundary="all submitted full steps, batch failure checks and compact observations; excludes compilation, construction and full snapshot audits",
    model_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"..","src","device_cellular.jl")))),
    kernels_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"..","src","device_cellular_kernels.jl")))))
mkpath(dirname(abspath(out)));mkdir(out)
open(joinpath(out,"benchmark.json"),"w") do io;JSON3.pretty(io,result);end
println(result)
