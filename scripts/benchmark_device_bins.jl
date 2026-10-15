#!/usr/bin/env julia
# Paired full-step comparison; wide bins are the same process with a coarser index.
include("device_cellular_setup.jl")
include("device_cellular_population.jl")
steps=parse(Int,get(ENV,"DEVICE_STEPS","200"));pairs=parse(Int,get(ENV,"DEVICE_PAIRS","5"))
sides=parse.(Int,split(get(ENV,"DEVICE_SIDES","8,16"),','))
min(steps,pairs)>0 && all(>(0),sides) || error("positive work required")
out=length(ARGS)>0 ? ARGS[1] : "outputs/device-cellular-bin-pairs"
ispath(out) && error("fresh output required")
mkpath(dirname(abspath(out)));mkdir(out)
function bin_sample(side,steps,width)
    w,p=device_cellular_population(side,backend;bin_width=width)
    initial=Interior.device_cellular_audit(Interior.device_cellular_snapshot(p));stats=nothing;elapsed=0.
    allocated=@allocated elapsed=@elapsed for step in 1:steps
        Interior.device_cellular_step!(w,p,.01;execution=:projected,check_status=false)
        if step%8==0 || step==steps;stats=Interior.device_cellular_statistics(p);end
    end
    snap=Interior.device_cellular_snapshot(p);final=Interior.device_cellular_audit(snap)
    result=(side=side,bin_width=width,steps=steps,seconds=elapsed,full_steps_per_second=steps/elapsed,
        host_allocation_bytes=allocated,stats=stats,audit=final,
        material_error=final.material-initial.material-final.supplied+initial.supplied+final.removed-initial.removed,
        potential_error=final.potential-initial.potential-final.supplied_potential+initial.supplied_potential+
            final.removed_potential-initial.removed_potential+final.work-initial.work+final.heat-initial.heat)
    result,snap
end
function bin_trajectory_difference(a,b)
    sa=a.state;sb=b.state;ia=findall(==(1),sa.active);ib=findall(==(1),sb.active)
    same=sa.active==sb.active && sa.ids==sb.ids
    function edges(s,active)
        live=Set(s.ids[active]);result=Set{Tuple{Int32,Int32}}()
        for i in active,k in axes(s.bonds,2)
            id=s.bond_ids[i,k]
            s.bonds[i,k]>0 && id in live && push!(result,minmax(s.ids[i],id))
        end
        result
    end
    ea=edges(sa,ia);eb=edges(sb,ib)
    (same_active_identities=same,maximum_position_error=same ? maximum(abs.(Float64.(sa.position[ia,:]).-Float64.(sb.position[ia,:]))) : nothing,
        rms_position_error=same ? sqrt(sum(abs2,Float64.(sa.position[ia,:]).-Float64.(sb.position[ia,:]))/(3length(ia))) : nothing,
        maximum_hidden_error=same ? maximum(abs.(Float64.(sa.hidden[ia,:]).-Float64.(sb.hidden[ia,:]))) : nothing,
        maximum_resource_error=maximum(abs.(Float64.(a.environment.fields).-Float64.(b.environment.fields))),
        bond_jaccard=isempty(union(ea,eb)) ? 1. : length(intersect(ea,eb))/length(union(ea,eb)))
end
# Same type signatures on a tiny population; no timed compilation.
for width in (4.15f0,1f0)
    w,p=device_cellular_population(2,backend;bin_width=width)
    Interior.device_cellular_step!(w,p,.01);Interior.device_cellular_statistics(p)
    w=nothing;p=nothing;GC.gc()
end
samples=Any[];comparisons=Any[]
for side in sides,pair in 1:pairs
    states=Dict{Float32,Any}();times=Dict{Float32,Float64}()
    for width in (isodd(pair) ? (4.15f0,1f0) : (1f0,4.15f0))
        sample,snapshot=bin_sample(side,steps,width)
        push!(samples,merge((pair=pair,),sample));states[width]=snapshot;times[width]=sample.seconds
        GC.gc()
    end
    push!(comparisons,(side=side,pair=pair,paired_speedup=times[4.15f0]/times[1f0],
        trajectory=bin_trajectory_difference(states[4.15f0],states[1f0])))
    empty!(states);GC.gc()
end
result=(source=Interior.DEVICE_CELLULAR_SOURCE_PIN,runtime=device_cellular_runtime(),mechanical_substeps=16,
    timing_boundary="full projected steps including batch8compact observations/failure checks; excludes construction, compilation, full audit downloads",samples=samples,comparisons=comparisons)
open(joinpath(out,"bin-pairs.json"),"w") do io;JSON3.pretty(io,result);end
println((source=result.source,comparisons=comparisons))
