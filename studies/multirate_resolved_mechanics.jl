module MultirateResolvedMechanics
using Interior, Random, Statistics
include("multirate_material_ecology.jl")
const Shared=MultirateMaterialEcology.BlockReservoirEcology
include("multirate_cached_controller.jl")
export resolved_fixture, resolved_step!, resolved_comparison

"""Experimental slow biology with original fast contact/motor resolution.
Controller, metabolism and lifecycle advance at H; passive mechanics and funded
motor debit retain delta/16, and persistent-bond capture is checked each delta.
Optional cached-chemical control retains fast neural dynamics. Biological
event times intentionally change; ratio1 is the exact native control.
"""
struct ResolvedCarrier{P}
    process::P
    fast_dt::Float32
    ratio::Int
    fast_controller::Bool
end
Base.getproperty(p::ResolvedCarrier,k::Symbol)=k in (:process,:fast_dt,:ratio,:fast_controller) ? getfield(p,k) : getproperty(getfield(p,:process),k)
Interior.device_support_descriptor(p::ResolvedCarrier)=Interior.device_support_descriptor(p.process)
function Interior.device_cellular_launch!(p::ResolvedCarrier,kernel,args...;ndrange)
    if p.fast_controller && kernel===Interior.dc_regulate!
        s,e,q,g,H=args
        Interior.device_cellular_launch!(p.process,mr_fast_neural!,s,e,q,g,p.fast_dt;ndrange)
        Interior.device_cellular_launch!(p.process,mr_request_only!,s,e,q,g,H;ndrange)
    else
        Interior.device_cellular_launch!(p.process,kernel,args...;ndrange)
    end
end
function Interior.device_cellular_mechanics!(p::ResolvedCarrier,H,mode::Val)
    native=p.process;g=native.geometry;C=Int(g.capacity)
    for j in 1:p.ratio
        if p.fast_controller && j>1
            Interior.device_cellular_launch!(native,mr_cached_sense!,native.state,native.environment,native.parameters,g;ndrange=C)
            controller=native.ownership.stages.controller
            Interior.device_cellular_launch!(native,Interior.dc_controller!,controller.state,controller.environment,controller.geometry;ndrange=(C,Int(g.hidden)+CELLULAR_ACTIONS))
            Interior.device_cellular_launch!(native,mr_fast_neural!,native.state,native.environment,native.parameters,g,p.fast_dt;ndrange=C)
        end
        Interior.device_cellular_mechanics!(native,p.fast_dt,mode)
        if j<p.ratio
            Interior.device_cellular_launch!(native,Interior.dc_propose_bonds!,native.state,native.environment,native.parameters,g;ndrange=C)
            Interior.device_cellular_launch!(native,Interior.dc_commit_bonds!,native.state,native.environment,native.parameters,g;ndrange=C)
        end
    end
    nothing
end
function resolved_fixture(;size=16,founders=8,capacity=128,seed=911,ratio=4,fast_dt=.001,fast_controller=false,backend=CPU(),energy_type=Float32)
    ratio>0||error("positive ratio required")
    w,p=Shared.device_habitat_pilot_fixture(backend;size,founders,capacity,seed,support_strategy=:materialized,energy_type)
    delta=Float32(fast_dt);H=delta*ratio
    isfinite(delta)&&delta>0&&H<=stability_limit(w,p)||error("invalid clocks")
    (world=w,process=p,carrier=ResolvedCarrier(p.cells,delta,ratio,fast_controller),H=H)
end
function resolved_step!(r;check_status=true)
    p=r.process.cells
    Interior.device_cellular_validate(r.world,p);check_status&&Interior.device_cellular_status!(p)
    p.attempted_steps+=1
    Interior.device_cellular_enqueue!(r.world,r.carrier,r.H,:projected,r.process.habitat,UInt32(p.attempted_steps))
    check_status&&Interior.device_cellular_status!(p)
    r
end
function resolved_comparison(;size=32,founders=32,capacity=128,fast_steps=32,fast_dt=.001,ratios=(1,2,4,8),seed=911,fast_controller=false,backend=CPU(),energy_type=Float32)
    fast_steps>0&&all(r->r>0&&fast_steps%r==0,ratios)||error("positive cadence-aligned duration required")
    ww,pp=Shared.device_habitat_pilot_fixture(backend;size,founders,capacity,seed,support_strategy=:materialized,energy_type)
    Interior.device_cellular_step!(ww,pp.cells,fast_dt;anabolism=pp.habitat)
    pp=nothing;ww=nothing
    for ratio in ratios
        r=resolved_fixture(;size,founders,capacity,seed,ratio,fast_dt,fast_controller,backend,energy_type);resolved_step!(r)
    end
    GC.gc()
    w,p=Shared.device_habitat_pilot_fixture(backend;size,founders,capacity,seed,support_strategy=:materialized,energy_type)
    initial=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
    baseline_seconds=@elapsed for _ in 1:fast_steps;Interior.device_cellular_step!(w,p.cells,fast_dt;anabolism=p.habitat);end
    reference=Interior.device_habitat_snapshot(p);rows=NamedTuple[]
    for ratio in ratios
        r=resolved_fixture(;size,founders,capacity,seed,ratio,fast_dt,fast_controller,backend,energy_type)
        seconds=@elapsed for _ in 1:fast_steps÷ratio;resolved_step!(r);end
        state=Interior.device_habitat_snapshot(r.process);audit=Interior.device_habitat_audit(state)
        exact=state.cells.state==reference.cells.state&&state.cells.environment==reference.cells.environment&&state.arrays==reference.arrays&&state.feedstock==reference.feedstock
        position=abs.(Float64.(state.cells.state.position).-Float64.(reference.cells.state.position));position=min.(position,4 .-position)
        push!(rows,(ratio,seconds,complete_step_speedup=baseline_seconds/seconds,duration=fast_steps*Float64(r.carrier.fast_dt),
            neural_updates=fast_controller ? fast_steps : fast_steps÷ratio,
            hidden_relative_l2=sqrt(sum(abs2,Float64.(state.cells.state.hidden).-Float64.(reference.cells.state.hidden))/max(sum(abs2,Float64.(reference.cells.state.hidden)),eps())),
            actual_mechanical_substep=Float64(r.carrier.fast_dt)/p.cells.parameters.mechanical_substeps,
            biological_events=fast_steps÷ratio,bond_capture_events=fast_steps,
            all_arrays_exact=exact,identities_equal=state.cells.state.ids==reference.cells.state.ids&&state.cells.state.active==reference.cells.state.active,
            audit,material_error=audit.material-initial.material,potential_error=audit.potential+audit.work+audit.heat-initial.potential,
            position_max=maximum(position),
            free_error=Shared.field_error(state.cells.environment.fields,reference.cells.environment.fields,prod(p.cells.geometry.spacing)),
            bound_error=Shared.field_error(state.arrays.bound,reference.arrays.bound,prod(p.cells.geometry.spacing))))
    end
    (recipe=(size,founders,capacity,fast_steps,fast_dt,seed,fast_controller),cached_controller_source_sha256=CACHED_CONTROLLER_SOURCE_SHA256,source=Interior.DEVICE_CELLULAR_SOURCE_PIN,baseline_seconds,rows)
end
end
