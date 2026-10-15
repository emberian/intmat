#!/usr/bin/env julia
include("device_cellular_setup.jl")
include("device_cellular_population.jl")
side=parse(Int,get(ENV,"DEVICE_CELLS_SIDE","8"));steps=parse(Int,get(ENV,"DEVICE_STEPS","100"))
min(side,steps)>0 || error("positive side/steps required")
out=length(ARGS)>0 ? ARGS[1] : "outputs/device-cellular-substeps"
ispath(out) && error("fresh output required")
mkpath(dirname(abspath(out)));mkdir(out)
function sample_substeps(side,steps,n,backend)
    w,p=device_cellular_population(side,backend;mechanical_substeps=n)
    initial=Interior.device_cellular_audit(Interior.device_cellular_snapshot(p))
    elapsed=0.;allocation=0;failure=nothing
    allocation=@allocated elapsed=@elapsed try
        for step in 1:steps
            Interior.device_cellular_step!(w,p,.01;check_status=false)
            step%8==0 && Interior.device_cellular_status!(p)
        end
        Interior.device_cellular_status!(p)
    catch exception
        exception isa ArgumentError || rethrow()
        failure=sprint(showerror,exception)
    end
    snap=Interior.device_cellular_snapshot(p;check_status=false)
    code=Int(only(snap.environment.status));complete=Int(snap.environment.metadata[6])
    stats=code==0 ? Interior.device_cellular_statistics(p) : nothing
    audit=code==0 ? Interior.device_cellular_audit(snap) : nothing
    result=(substeps=n,seconds=elapsed,host_allocation_bytes=allocation,completed_steps=complete,
        full_steps_per_second=code==0 ? steps/elapsed : nothing,status=code,failure=failure,stats=stats,audit=audit,
        material_error=code==0 ? audit.material-initial.material-audit.supplied+initial.supplied+audit.removed-initial.removed : nothing)
    result,snap
end
# Compile before sampling; same types, independently constructed small population.
w,p=device_cellular_population(2,backend;mechanical_substeps=16)
Interior.device_cellular_step!(w,p,.01);p=nothing;w=nothing;GC.gc()
results=Any[];snapshots=Dict{Int,Any}()
for n in (1,2,4,8,16)
    result,snap=sample_substeps(side,steps,n,backend)
    push!(results,result)
    # Keep compact actual final physical state, not full capacity genomes/scratch.
    result.status==0 && (snapshots[n]=(active=snap.state.active,ids=snap.state.ids,position=snap.state.position,
        mass=snap.state.mass,metabolites=snap.state.metabolites,energy=snap.state.energy,fields=snap.environment.fields,
        bonds=snap.state.bonds,bond_ids=snap.state.bond_ids))
    snap=nothing;GC.gc()
end
reference=snapshots[16];comparisons=Any[]
for n in sort!(collect(keys(snapshots)))
    x=snapshots[n];same_ids=x.active==reference.active && x.ids==reference.ids
    active=findall(==(1),reference.active)
    edges(y)=Set((min(Int(y.ids[i]),Int(y.bond_ids[i,b])),max(Int(y.ids[i]),Int(y.bond_ids[i,b]))) for i in active for b in axes(y.bonds,2) if y.bonds[i,b]>0)
    a=edges(x);b=edges(reference)
    push!(comparisons,(substeps=n,same_active_identities=same_ids,
        maximum_position_error=same_ids ? maximum(abs.(Float64.(x.position[active,:]).-Float64.(reference.position[active,:]))) : nothing,
        rms_position_error=same_ids ? sqrt(sum(abs2,Float64.(x.position[active,:]).-Float64.(reference.position[active,:]))/(3length(active))) : nothing,
        maximum_resource_error=maximum(abs.(Float64.(x.fields).-Float64.(reference.fields))),
        bond_jaccard=isempty(union(a,b)) ? 1. : length(intersect(a,b))/length(union(a,b))))
end
result=(source=Interior.DEVICE_CELLULAR_SOURCE_PIN,runtime=device_cellular_runtime(),side=side,steps=steps,dt=.01,
    timing_boundary="full asynchronous pipeline and status check every8steps; excludes compile/construction and full final snapshots",results=results,comparisons_to_16=comparisons)
open(joinpath(out,"substeps.json"),"w") do io;JSON3.pretty(io,result);end
println(result)
