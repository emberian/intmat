#!/usr/bin/env julia
using Interior,TOML,JSON3

"""Matched mechanism toggles in two independently founded gardens.
These are developmental/ecological trajectories, not adaptation assays.
"""
function sorptive_study(profile,output;steps=10000,capture_every=100,seeds=(211,307),
        conditions=("full","sorption_only","transport_only","catalytic_only"))
    ispath(output) && error("use a fresh output directory")
    mkpath(output);config=TOML.parsefile(profile)
    report=Dict{String,Any}("kind"=>"sorptive matrix ecological exploration; no adaptation claim",
        "loaded_sources"=>Interior.arena_sources(),"steps"=>steps,"capture_every"=>capture_every,"arms"=>Dict())
    for seed in seeds,condition in conditions
        cfg=deepcopy(config);cfg["seed"]=seed
        mechanism=replace(condition,r"^(residual_|shared_)"=>"")
        startswith(condition,"residual_") && (cfg["catalysis"]["shared_anabolism"]=false)
        startswith(condition,"shared_") && (cfg["catalysis"]["shared_anabolism"]=true)
        mechanism in ("full","sorption_only","transport_only","catalytic_only") || error("unknown mechanism $condition")
        mechanism in ("transport_only","catalytic_only") && (cfg["matrix"]["capacity_per_enzyme"]=0.)
        mechanism in ("sorption_only","catalytic_only") && (cfg["matrix"]["transport_barrier"]=0.)
        a=build_arena(cfg);name="$seed-$condition";directory=joinpath(output,name);mkpath(directory)
        frames=Any[];maxmass=0.;maxenergy=0.;peakbound=0.;peakpopulation=0;peakfile=""
        function capture()
            snapshot=arena_snapshot(a);file="frame-"*lpad(a.steps,7,'0')*".json"
            open(io->JSON3.write(io,snapshot),joinpath(directory,file),"w")
            open(io->JSON3.write(io,snapshot),joinpath(directory,"latest.json"),"w")
            m=snapshot["metrics"];bound=sum(m["matrix"]["bound_mass"])
            if bound>peakbound;peakbound=bound;peakfile=file;end
            peakpopulation=max(peakpopulation,m["population"])
            maxmass=max(maxmass,abs(m["material_error"]));maxenergy=max(maxenergy,abs(m["chemical_energy_error"]))
            push!(frames,Dict("file"=>file,"time"=>a.world.time,"step"=>a.steps))
            open(joinpath(directory,"metrics.jsonl"),"a") do io
                JSON3.write(io,Dict("step"=>a.steps,"time"=>a.world.time,"metrics"=>m));write(io,'\n')
            end
        end
        capture();checkpoint!(a,joinpath(directory,"initial.jls"));started=time()
        for i in 1:steps
            arena_step!(a)
            (i%capture_every==0 || i==steps) && capture()
            i==steps÷2 && checkpoint!(a,joinpath(directory,"middle.jls"))
        end
        elapsed=time()-started;checkpoint!(a,joinpath(directory,"final.jls"))
        open(io->JSON3.write(io,Dict("kind"=>"actual recorded trajectory","condition"=>condition,"seed"=>seed,"frames"=>frames)),joinpath(directory,"sequence.json"),"w")
        final=arena_snapshot(a)["metrics"]
        report["arms"][name]=Dict("config"=>cfg,"final"=>final,"peak_bound"=>peakbound,"peak_bound_frame"=>peakfile,
            "peak_population"=>peakpopulation,"max_abs_material_error"=>maxmass,"max_abs_chemical_error"=>maxenergy,"elapsed_seconds_including_captures"=>elapsed)
        open(io->JSON3.write(io,report),joinpath(output,"report.json"),"w")
        println(name,": pop=",final["population"]," births=",final["births"]," bonds=",final["bonds"]," retained=",sum(final["matrix"]["bound_mass"])," peak_retained=",peakbound," elapsed=",round(elapsed;digits=2));flush(stdout)
    end
    report
end

if abspath(PROGRAM_FILE)==(@__FILE__)
    length(ARGS)>=2 || error("usage: sorptive_garden.jl PROFILE OUTPUT [STEPS=10000] [CAPTURE=100]")
    sorptive_study(ARGS[1],abspath(ARGS[2]);steps=length(ARGS)>2 ? parse(Int,ARGS[3]) : 10000,
        capture_every=length(ARGS)>3 ? parse(Int,ARGS[4]) : 100)
end
