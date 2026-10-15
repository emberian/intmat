#!/usr/bin/env julia
"""Recipe-derived scale estimates, not a simulation or a stability verdict."""
module GardenOperatingMap
using TOML, JSON3, SHA
# Pure parameter constructor: no simulation allocation or backend compilation.
include(joinpath(@__DIR__,"..","src","cellular_initialization.jl"))

function operating_map(cfg; catalytic_rate=cfg["catalysis"]["catalytic_rate"],
        transport_barrier=cfg["matrix"]["transport_barrier"])
    defaults=cellular_parameters()
    q=merge(defaults,(; (Symbol(k)=>v for (k,v) in cfg["cellular"] if Symbol(k) in keys(defaults))...))
    volume=prod(cfg["grid"]["extent"])
    L=q.exchange_radius
    h=cfg["catalysis"];m=cfg["matrix"];f=cfg["feedstock"]
    # Source rates are world totals. These estimates use spatial means before
    # depletion; the actual source masks, movement and uptake are heterogeneous.
    source=Float64.(f["source_rates"])./volume
    initial=Float64.(f["initial"])
    all(>(0),initial) && catalytic_rate>0 || error("positive raw inventory and catalytic rate required")
    enzyme=source.*(h["half_saturation"].+initial)./(catalytic_rate.*initial)
    total=sum(enzyme)
    permeability=inv(1+transport_barrier*total)
    # L²/D convention, not an RMS first-passage time; state this for comparison.
    diffusion=Float64.(cfg["resources"]["diffusion"])[1:2]
    clock(D)=D>0 ? L^2/D : nothing
    bare=clock.(diffusion);retained=clock.(diffusion.*permeability)
    decay=q.plasticity_decay;eta=q.plasticity_rate
    mark_fast=decay+eta>0 ? inv(decay+eta) : nothing
    mark_slow=decay>0 ? inv(decay) : nothing
    capacity=m["capacity_per_enzyme"].*enzyme
    free=Float64.(cfg["resources"]["initial"])[1:2]
    # Fixed TOTAL soluble+bound and capacity local linearization at B=0.
    binding_rate=m["adsorption_rate"].*(free.+capacity).+m["desorption_rate"]
    binding=map(x->x>0 ? inv(x) : nothing,binding_rate)
    sourcespan=[source[k]>0 ? initial[k]/source[k] : nothing for k in 1:2]
    budgetspan=[f["source_rates"][k]>0 ? f["source_budgets"][k]/f["source_rates"][k] : nothing for k in 1:2]
    (;catalytic_rate,transport_barrier,world_volume=volume,exchange_radius=L,
        mean_raw_source_rate=source,supply_matching_enzyme= enzyme,
        supply_matching_total_enzyme=total,permeability_at_supply_match=permeability,
        transport_resistance_at_supply_match=transport_barrier*total,
        matching_bound_capacity=capacity,
        capacity_over_initial_free=[free[k]>0 ? capacity[k]/free[k] : nothing for k in 1:2],
        clocks=(neural_response=q.response_time,adaptation=q.adaptation_time,
            plastic_mark_relaxation_range=(mark_fast,mark_slow),
            bare_exchange_diffusion=bare,matched_exchange_diffusion=retained,
            minimum_nominal_exchange_travel=q.speed>0 ? L/q.speed : nothing,
            catalyst_decay=h["decay_rate"]>0 ? inv(h["decay_rate"]) : nothing,
            free_desorption=m["desorption_rate"]>0 ? inv(m["desorption_rate"]) : nothing,
            binding_local_at_empty_bound=binding,
            initial_raw_over_mean_source=sourcespan,finite_source_duration=budgetspan),
        ratios=(adaptation_over_response=q.adaptation_time/q.response_time,
            diffusion_over_adaptation=[x===nothing ? nothing : x/q.adaptation_time for x in retained],
            catalyst_lifetime_over_adaptation=h["decay_rate"]>0 ? inv(h["decay_rate"]*q.adaptation_time) : nothing),
        exchange_lengths_per_adaptation=[sqrt(D*permeability*q.adaptation_time)/L for D in diffusion])
end

function report(profile,output)
    ispath(output) && error("fresh output required")
    cfg=TOML.parsefile(profile)
    baseline=operating_map(cfg)
    scenarios=[operating_map(cfg;catalytic_rate=k,transport_barrier=b)
        for k in (6.,.6,.3), b in (20.,100.,500.)]
    result=(kind="analytical operating map; no simulated trajectories",
        profile=relpath(profile,dirname(@__DIR__)),profile_sha256=bytes2hex(sha256(read(profile))),
        source_sha256=bytes2hex(sha256(read(@__FILE__))),
        defaults_sha256=bytes2hex(sha256(read(joinpath(dirname(@__DIR__),"src","cellular_initialization.jl")))),
        assumptions=["L^2/D characteristic diffusion convention; no factor of 2d or first-passage claim",
            "Spatial mean source and initial raw concentration frozen; no uptake or heterogeneous depletion",
            "Supply-matching enzyme is a balance estimate, not predicted expression or an optimal phenotype",
            "Bound-capacity estimate assumes zero cross_binding; this report requires it",
            "Binding clock linearizes fixed total soluble+bound at B=0; changing free concentration changes the clock",
            "Travel uses maximum nominal speed; actual paid motion may be slower",
            "Mark relaxation freezes hidden activity in [-1,1] and leaves clipping inactive",
            "Candidate contrasts are not validated viable regimes or criticality predictions"],
        baseline,scenarios=vec(scenarios))
    cfg["matrix"]["cross_binding"]==0 || error("cross-binding estimate not implemented")
    open(io->JSON3.pretty(io,result),output,"w")
    result
end
end
if abspath(PROGRAM_FILE)==(@__FILE__) || (isdefined(Main,:INTERIOR_RUNNING_STUDY) && Main.INTERIOR_RUNNING_STUDY=="garden_operating_map")
    length(ARGS)==2 || error("usage: garden_operating_map.jl PROFILE FRESH_OUTPUT.json")
    GardenOperatingMap.report(ARGS...)
end
