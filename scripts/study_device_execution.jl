#!/usr/bin/env julia
include("device_cellular_setup.jl")
include("device_cellular_population.jl")
steps=parse(Int,get(ENV,"DEVICE_STEPS","100"));sides=parse.(Int,split(get(ENV,"DEVICE_SIDES","8,16"),','))
steps>0 && all(>(0),sides) || error("positive work required")
out=length(ARGS)>0 ? ARGS[1] : "outputs/device-cellular-execution"
ispath(out) && error("fresh output required")
mkpath(dirname(abspath(out)));mkdir(out)
function execution_sample(side,steps,recipe,bin_width)
    w,p=device_cellular_population(side,backend;bin_width)
    initial=Interior.device_cellular_audit(Interior.device_cellular_snapshot(p));stats=nothing;elapsed=0.
    allocated=@allocated elapsed=@elapsed for i in 1:steps
        Interior.device_cellular_step!(w,p,.01;execution=recipe,check_status=false)
        if i%8==0 || i==steps;stats=Interior.device_cellular_statistics(p);end
    end
    snap=Interior.device_cellular_snapshot(p);audit=Interior.device_cellular_audit(snap)
    (side=side,recipe=recipe,bin_width=bin_width,steps=steps,seconds=elapsed,full_steps_per_second=steps/elapsed,
        host_allocation_bytes=allocated,stats=stats,audit=audit,
        material_error=audit.material-initial.material-audit.supplied+initial.supplied+audit.removed-initial.removed),snap
end
function compare_execution(reference,actual)
    errors=Dict{String,Float64}();integers=String[]
    for group in (:state,:environment),key in keys(getproperty(reference,group))
        a=getproperty(getproperty(reference,group),key);b=getproperty(getproperty(actual,group),key)
        if eltype(a)<:AbstractFloat
            errors["$group.$key"]=isempty(a) ? 0. : maximum(abs.(Float64.(a).-Float64.(b)))
        elseif a!=b
            push!(integers,"$group.$key")
        end
    end
    (maximum_float_error=maximum(values(errors)),integer_mismatches=integers,errors=errors)
end
for recipe in (:staged,:projected)
    w,p=device_cellular_population(2,backend)
    Interior.device_cellular_step!(w,p,.01;execution=recipe)
    w=nothing;p=nothing;GC.gc()
end
results=Any[];comparisons=Any[]
for side in sides,bin_width in parse.(Float32,split(get(ENV,"DEVICE_BIN_WIDTHS","4.15,1"),','))
    reference=nothing
    for recipe in (:staged,:projected)
        result,snap=execution_sample(side,steps,recipe,bin_width);push!(results,result)
        if recipe==:staged;reference=snap
        else;push!(comparisons,(side=side,bin_width=bin_width,recipe=recipe,comparison=compare_execution(reference,snap)))
        end
        snap=nothing;GC.gc()
    end
    reference=nothing;GC.gc()
end
result=(source=Interior.DEVICE_CELLULAR_SOURCE_PIN,runtime=device_cellular_runtime(),
    timing_boundary="all full steps with batch8compactobservations/failure checks; excludes compilation/construction/full snapshots",results=results,comparisons=comparisons)
open(joinpath(out,"execution.json"),"w") do io;JSON3.pretty(io,result);end
println(result)
