module MultirateMaterialEcology
using Interior, Random, Statistics
using Interior.KernelAbstractions: @kernel, @index
# Reuse only its checked shared fixture / fresh no-diffusion carrier helper.
# No block-grid approximation or block transfer is executed in this study.
include("block_reservoir_ecology.jl")
const fixture=BlockReservoirEcology.device_habitat_pilot_fixture
export multirate_fixture, multirate_step!, multirate_comparison

mutable struct SlowMaterialHook{H}
    habitat::H
    finish::Bool
    interval::Float32
end
@kernel function restore_held_material!(a,e)
    site=@index(Global,Linear)
    if e.status[1]==0;e.permeability[site]*=a.permeability[site];end
end
Interior.device_anabolic_argument(h::SlowMaterialHook)=Interior.device_anabolic_argument(h.habitat)
Interior.device_anabolic_begin!(::SlowMaterialHook,w,p,dt)=nothing
Interior.device_anabolic_sense!(h::SlowMaterialHook,w,p,dt)=Interior.device_anabolic_sense!(h.habitat,w,p,dt)
Interior.device_anabolic_prepare!(h::SlowMaterialHook,w,p,dt)=Interior.device_anabolic_prepare!(h.habitat,w,p,dt)
Interior.device_anabolic_inherit!(h::SlowMaterialHook,w,p,step)=Interior.device_anabolic_inherit!(h.habitat,w,p,step)
function Interior.device_anabolic_finish!(h::SlowMaterialHook,w,p,dt)
    native=h.habitat;V=prod(Int.(p.geometry.dims))
    if h.finish
        Interior.device_cellular_launch!(p,Interior.dh_react_sorb!,native.arrays,native.feedstock.arrays,
            p.environment,native.parameters,native.feedstock.parameters,p.geometry,h.interval;ndrange=V)
        Interior.device_cellular_launch!(p,Interior.dh_permeability!,native.arrays,p.environment,p.geometry,native.parameters;ndrange=V)
    else
        Interior.device_cellular_launch!(p,restore_held_material!,native.arrays,p.environment;ndrange=V)
    end
    Interior.device_cellular_launch!(p,Interior.dh_validate_complete!,native.arrays,native.feedstock.arrays,
        p.environment,p.geometry;ndrange=max(V,Int(p.geometry.capacity)))
    nothing
end

mutable struct MultirateEcology{W,P,H,Q,F}
    world::W
    process::P
    hook::H
    free_parameters::Q
    raw_parameters::F
    ratio::Int
    fast_dt::Float32
    complete_outer_steps::Int
    fresh_catalyst_release_bound::Float64
end
function multirate_fixture(;size=16,founders=8,capacity=128,seed=911,ratio=4,fast_dt=.001)
    ratio>0||error("positive cadence ratio required")
    w,p=fixture(CPU();size,founders,capacity,seed,support_strategy=:materialized)
    q=p.cells.parameters;fq=p.habitat.feedstock.parameters
    all(iszero,(q.source_rates...,q.sink_rates...,fq.source_rates...,fq.sink_rates...))||error("closed-only study")
    dt=Float32(fast_dt);H=dt*ratio
    isfinite(dt)&&dt>0&&isfinite(H)&&H>0||error("positive represented clocks required")
    dt<=Interior.device_cellular_timestep_limit(p.cells)||error("fast clock exceeds bound")
    H<=stability_limit(w,p)||error("outer transport clock exceeds diffusion/neural guard")
    split=BlockReservoirEcology.without_diffusion(p)
    hook=SlowMaterialHook(split.habitat,false,H)
    MultirateEcology(w,split,hook,q,fq,ratio,dt,0,0.)
end
function secretion_total(p)
    a=p.habitat.arrays
    sum(Float64,a.cell_ledger[:,7:8])+sum(Float64,a.cell_ledger_low[:,7:8])
end
function multirate_step!(m)
    p=m.process.cells;h=m.process.habitat;g=p.geometry;V=prod(Int.(g.dims));H=m.hook.interval
    Interior.device_cellular_status!(p)
    before=secretion_total(m.process)
    # Same ordering as native begin: raw transport, then free transport.
    Interior.device_resource_phase!(p,h.feedstock.argument,m.raw_parameters,g,H)
    Interior.device_cellular_launch!(p,Interior.dc_diffuse!,h.feedstock.argument,m.raw_parameters,g,H;ndrange=(V,3))
    Interior.device_cellular_launch!(p,Interior.dc_commit_resources!,h.feedstock.argument;ndrange=(V,3))
    Interior.device_cellular_launch!(p,Interior.dc_diffuse!,p.environment,m.free_parameters,g,H;ndrange=(V,3))
    Interior.device_cellular_launch!(p,Interior.dc_commit_resources!,p.environment;ndrange=(V,3))
    for j in 1:m.ratio
        m.hook.finish=j==m.ratio
        Interior.device_cellular_step!(m.world,p,m.fast_dt;anabolism=m.hook,check_status=false)
    end
    Interior.device_cellular_status!(p)
    m.complete_outer_steps+=1
    # Conservative bound for fresh catalyst's excess reaction exposure. It is
    # not a bound on the fully coupled nonlinear trajectory or all split error.
    fresh=max(0.,secretion_total(m.process)-before)
    m.fresh_catalyst_release_bound+=Float64(h.parameters.catalytic_rate)*Float64(H)*fresh
    m
end

function multirate_comparison(;size=32,founders=32,capacity=128,fast_steps=32,fast_dt=.001,ratios=(1,2,4,8),seed=911)
    fast_steps>0&&all(r->r>0&&fast_steps%r==0,ratios)||error("positive cadence-aligned duration required")
    # Target-shape compilation matters for CPU-KA workgroup specializations.
    ww,pp=fixture(CPU();size,founders,capacity,seed,support_strategy=:materialized)
    Interior.device_cellular_step!(ww,pp.cells,fast_dt;anabolism=pp.habitat)
    pp=nothing;ww=nothing
    for ratio in ratios
        warm=multirate_fixture(;size,founders,capacity,seed,ratio,fast_dt)
        multirate_step!(warm)
        warm=nothing
    end
    GC.gc()
    w,baseline=fixture(CPU();size,founders,capacity,seed,support_strategy=:materialized)
    initial=Interior.device_habitat_audit(Interior.device_habitat_snapshot(baseline))
    baseline_seconds=@elapsed for _ in 1:fast_steps
        Interior.device_cellular_step!(w,baseline.cells,fast_dt;anabolism=baseline.habitat)
    end
    reference=Interior.device_habitat_snapshot(baseline);base_final=Interior.device_habitat_audit(reference)
    rows=NamedTuple[]
    for ratio in ratios
        m=multirate_fixture(;size,founders,capacity,seed,ratio,fast_dt)
        elapsed=@elapsed for _ in 1:fast_steps÷ratio;multirate_step!(m);end
        state=Interior.device_habitat_snapshot(m.process);audit=Interior.device_habitat_audit(state)
        exact=state.cells.state==reference.cells.state&&state.cells.environment==reference.cells.environment&&
            state.arrays==reference.arrays&&state.feedstock==reference.feedstock
        volume=prod(m.process.cells.geometry.spacing)
        push!(rows,(ratio,seconds=elapsed,complete_step_speedup=baseline_seconds/elapsed,
            duration=fast_steps*Float64(m.fast_dt),fast_steps,outer_steps=m.complete_outer_steps,
            all_arrays_exact=exact,identities_equal=state.cells.state.ids==reference.cells.state.ids&&state.cells.state.active==reference.cells.state.active,
            audit,material_error=audit.material-initial.material,
            potential_error=audit.potential+audit.work+audit.heat-initial.potential,
            fresh_catalyst_release_bound=m.fresh_catalyst_release_bound,
            old_catalyst_relative_hold_bound=expm1(Float64(m.process.habitat.parameters.decay_rate)*Float64(m.hook.interval)),
            free_error=BlockReservoirEcology.field_error(state.cells.environment.fields,reference.cells.environment.fields,volume),
            enzyme_error=BlockReservoirEcology.field_error(reinterpret(Float32,state.arrays.enzyme_bits),reinterpret(Float32,reference.arrays.enzyme_bits),volume),
            bound_error=BlockReservoirEcology.field_error(state.arrays.bound,reference.arrays.bound,volume)))
    end
    (recipe=(size,founders,capacity,fast_steps,fast_dt,seed),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,
        baseline_seconds,initial,base_final,rows)
end
end
