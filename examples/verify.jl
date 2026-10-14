"""Reproducible numerical checks and matched causal ablations of the supplied toy laws.

Run `julia --project=. examples/verify.jl`. Results describe these invented laws,
not physical calibration, learning, or intelligence. No renderer is involved.
"""
module SceneVerification
using Interior, TOML, JSON3, Statistics

function configuration(; dimension=3, resolution=16, particles=64)
    config=TOML.parsefile(joinpath(@__DIR__,"choir.toml"))
    config["simulation"]["dimension"]=dimension
    config["simulation"]["resolution"]=resolution
    config["particles"]["count"]=particles
    config
end

inventory(scene)=sum(Float64,scene.world.fields.nutrient)+
    sum(Float64,scene.world.fields.deposit)+
    (scene.particles===nothing ? 0.0 : sum(Float64,scene.particles.storage))

function ablate!(scene, variant)
    if variant == "no_neural_secretion"
        p=scene.particles
        p===nothing && error("secretion ablation requires particles")
        r=p.rule
        # Keep recurrence, sensing and motion; remove latent output to secretion.
        p.rule=ParticleRule(r.recurrence,r.sensing,r.motion,zero.(r.secretion),
                           r.uptake_rate,r.secretion_rate,r.speed,r.relaxation)
    end
    scene
end

function run_scenario(config; dt=0.02, duration=8.0, variant="coupled")
    config=deepcopy(config)
    if variant == "no_particles"
        config["particles"]["count"]=0
    elseif variant == "no_configured_material_feedback"
        config["excitable"]["material_feedback"]=0.0
        config["transport"]["material_feedback"]=0.0
    end
    scene=ablate!(build_scene(config),variant)
    initial=inventory(scene)
    positions=scene.particles===nothing ? nothing : copy(scene.particles.positions)
    steps=round(Int,duration/dt)
    isapprox(steps*dt,duration;atol=1e-10) || error("duration must be a multiple of dt")
    max_error=0.0
    elapsed=@elapsed for i in 1:steps
        step!(scene.world,scene.processes,dt)
        max_error=max(max_error,abs(inventory(scene)-initial)/initial)
    end
    f=scene.world.fields
    finite=all(a->all(isfinite,a),values(f))
    nonnegative=minimum(f.nutrient)>=0 && minimum(f.deposit)>=0
    motion=0.0; storage=0.0
    if scene.particles!==nothing
        p=scene.particles
        finite &= all(isfinite,p.positions) && all(isfinite,p.hidden) && all(isfinite,p.storage)
        nonnegative &= minimum(p.storage;init=0f0)>=0
        displacement=p.positions-positions
        for d in 1:3
            n=scene.world.grid.size[d]
            displacement[d,:] .= mod.(displacement[d,:].+n/2,n).-n/2
            displacement[d,:] .*= scene.world.grid.spacing[d]
        end
        motion=sqrt(sum(abs2,displacement)/max(1,length(p.ids)))
        storage=sum(Float64,p.storage)
    end
    metrics=Dict("variant"=>variant,"dimension"=>config["simulation"]["dimension"],
        "resolution"=>config["simulation"]["resolution"],"dt"=>dt,"duration"=>duration,
        "steps"=>steps,"initial_inventory"=>initial,"final_inventory"=>inventory(scene),
        "maximum_relative_mass_error"=>max_error,"finite"=>finite,"nonnegative_inventory"=>nonnegative,
        "deposit"=>sum(Float64,f.deposit),"deposit_mean"=>mean(Float64.(f.deposit)),
        "activity_rms"=>sqrt(mean(abs2,Float64.(f.activity))),
        "particle_storage"=>storage,"particle_rms_displacement"=>motion,"elapsed_seconds"=>elapsed)
    (scene=scene,metrics=metrics)
end

function comparison(a,b)
    fa,fb=a.scene.world.fields,b.scene.world.fields
    Dict("deposit_relative_change"=>(b.metrics["deposit"]-a.metrics["deposit"])/max(a.metrics["deposit"],eps()),
         "deposit_field_rms_difference"=>sqrt(mean(abs2,Float64.(fb.deposit).-Float64.(fa.deposit))),
         "activity_field_rms_difference"=>sqrt(mean(abs2,Float64.(fb.activity).-Float64.(fa.activity))))
end

function main(output=joinpath(@__DIR__,"..","outputs","verification.json"))
    reports=Any[]
    for dimension in (2,3)
        config=configuration(;dimension,resolution=dimension==2 ? 24 : 16)
        runs=[run_scenario(config;dt) for dt in (0.02,0.01,0.005)]
        coupled=runs[1]
        variants=[run_scenario(config;variant) for variant in
            ("no_particles","no_configured_material_feedback","no_neural_secretion")]
        all(r->r.metrics["finite"] && r.metrics["nonnegative_inventory"] &&
            r.metrics["maximum_relative_mass_error"]<1e-5,vcat(runs,variants)) ||
            error("scene conservation or validity check failed")
        push!(reports,Dict("dimension"=>dimension,"configuration"=>config,
            "time_refinement"=>[r.metrics for r in runs],
            "time_refinement_differences"=>[comparison(runs[1],runs[2]),comparison(runs[2],runs[3])],
            "ablations"=>[Dict("metrics"=>r.metrics,"difference_from_coupled"=>comparison(coupled,r)) for r in variants]))
        println("Validated ",dimension,"D, mass residual ",maximum(r.metrics["maximum_relative_mass_error"] for r in vcat(runs,variants)))
    end
    result=Dict("julia_version"=>string(VERSION),"backend"=>"CPU","seed"=>17,
        "model"=>"Fixed invented recurrent particle and excitable grid rules; no training or evolution",
        "scope"=>"Inventory conservation, finite state, temporal refinement and matched ablations; no physical calibration",
        "ablation_notes"=>Dict(
            "no_particles"=>"Remove carrier and its exchange; retain identical grid initial state and all grid rules.",
            "no_configured_material_feedback"=>"Zero excitable threshold and nutrient mobility feedback coefficients. Deposit-dependent excitable diffusion and wave attenuation remain.",
            "no_neural_secretion"=>"Zero latent secretion output weights, replacing state-dependent gate with 0.5. Recurrence, sensing, motion, uptake and secretion rate remain."),
        "reports"=>reports)
    mkpath(dirname(output))
    open(output,"w") do io; JSON3.pretty(io,result); end
    println("Saved ",abspath(output))
    result
end
end

if abspath(PROGRAM_FILE)==@__FILE__
    SceneVerification.main(isempty(ARGS) ? joinpath(@__DIR__,"..","outputs","verification.json") : abspath(ARGS[1]))
end
