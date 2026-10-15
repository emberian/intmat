# Fixed-duration discretization study; no production law changes.
using Interior, Random, Statistics
include(joinpath(@__DIR__,"..","scripts","device_habitat_pilot_fixture.jl"))

function habitat_timestep_difference(a,b,extent)
    position=abs.(Float64.(a.cells.state.position).-Float64.(b.cells.state.position))
    for d in 1:3
        position[:,d].=min.(position[:,d],abs.(extent[d].-position[:,d]))
    end
    err(x,y)=(maxabs=maximum(abs.(Float64.(x).-Float64.(y))),
        relative_l2=sqrt(sum(abs2,Float64.(x).-Float64.(y))/max(sum(abs2,Float64.(y)),eps())))
    (identities_equal=a.cells.state.ids==b.cells.state.ids&&a.cells.state.active==b.cells.state.active,
        position_max=maximum(position),hidden=err(a.cells.state.hidden,b.cells.state.hidden),
        metabolites=err(a.cells.state.metabolites,b.cells.state.metabolites),
        free=err(a.cells.environment.fields,b.cells.environment.fields),
        raw=err(a.feedstock.fields,b.feedstock.fields),
        enzyme=err(reinterpret(Float32,a.arrays.enzyme_bits),reinterpret(Float32,b.arrays.enzyme_bits)),
        bound=err(a.arrays.bound,b.arrays.bound))
end

"""Run fresh same-seed recipes at exactly matched rounded Float32 durations.
The first candidate is the reference. Timings include batch status checks, exclude
construction and detached final observations. A cold first run is not warm speed.
"""
function habitat_timestep_study(backend;size=16,founders=8,capacity=128,
        reference_dt=.0005,reference_steps=32,factors=(1,2,4,8,16),seed=911)
    all(f->f>0&&reference_steps%f==0,factors)||error("factors must divide reference_steps")
    first(factors)==1||error("first factor must be reference")
    rows=NamedTuple[];reference=nothing;reference_audit=nothing
    for factor in factors
        w,p=device_habitat_pilot_fixture(backend;size,founders,capacity,seed)
        dt=Float32(reference_dt)*factor;steps=reference_steps÷factor
        limit=stability_limit(w,p);dt<=limit||error("candidate exceeds declared bound")
        initial=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
        elapsed=@elapsed begin
            for k in 1:steps
                Interior.device_cellular_step!(w,p.cells,dt;anabolism=p.habitat,check_status=false)
                k%8==0&&Interior.device_cellular_status!(p.cells)
            end
            Interior.device_cellular_status!(p.cells)
        end
        snapshot=Interior.device_habitat_snapshot(p);audit=Interior.device_habitat_audit(snapshot)
        if reference===nothing
            reference=snapshot;reference_audit=audit
        end
        duration=steps*Float64(dt)
        push!(rows,(factor,dt=Float64(dt),steps,duration,declared_limit=limit,seconds=elapsed,
            steps_per_second=steps/elapsed,simulated_time_per_wall_second=duration/elapsed,
            material_error=audit.material-initial.material-audit.supplied+audit.removed,
            potential_error=audit.potential+audit.work+audit.heat-initial.potential-audit.supplied_potential+audit.removed_potential,
            work_difference=audit.work-reference_audit.work,heat_difference=audit.heat-reference_audit.heat,
            audit,difference=habitat_timestep_difference(snapshot,reference,p.cells.geometry.extent)))
        p=nothing;w=nothing;snapshot=nothing;GC.gc()
    end
    (recipe=(size,founders,capacity,seed,reference_dt,reference_steps),rows)
end
