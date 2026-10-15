#!/usr/bin/env julia
# Three matched exploratory conditions; no live runs are read or changed.
using Interior, TOML, JSON3, SHA

function catalytic_study(profile,output;steps=3000,capture_every=25)
    ispath(output) && error("use a fresh study output directory")
    mkpath(output);config=TOML.parsefile(profile)
    report=Dict{String,Any}("kind"=>"constructed catalytic habitat exploration; not an adaptation assay",
        "loaded_sources"=>Interior.arena_sources(),"steps"=>steps,"capture_every"=>capture_every,"arms"=>Dict())
    for condition in ("regulated","catalysis_off","fixed_readout")
        cfg=deepcopy(config)
        condition=="catalysis_off" && (cfg["catalysis"]["catalytic_rate"]=0.)
        condition=="fixed_readout" && (cfg["catalysis"]["regulation"]=false)
        a=build_arena(cfg);directory=joinpath(output,condition);mkpath(directory)
        frames=Any[];peak_matrix=-Inf;peak_file="";max_mass_error=0.;max_energy_error=0.
        function capture()
            snapshot=arena_snapshot(a);name="frame-"*lpad(a.steps,7,'0')*".json"
            open(io->JSON3.write(io,snapshot),joinpath(directory,name),"w")
            open(io->JSON3.write(io,snapshot),joinpath(directory,"latest.json"),"w")
            metrics=snapshot["metrics"];matrix=sum(metrics["catalytic"]["enzyme_mass"])
            if matrix>peak_matrix && metrics["population"]>0;peak_matrix=matrix;peak_file=name;end
            max_mass_error=max(max_mass_error,abs(metrics["material_error"]))
            max_energy_error=max(max_energy_error,abs(metrics["chemical_energy_error"]))
            push!(frames,Dict("file"=>name,"step"=>a.steps,"time"=>a.world.time))
            open(joinpath(directory,"metrics.jsonl"),"a") do io
                JSON3.write(io,Dict("step"=>a.steps,"time"=>a.world.time,"metrics"=>metrics));write(io,'\n')
            end
        end
        capture();checkpoint!(a,joinpath(directory,"initial.jls"))
        for i in 1:steps
            arena_step!(a)
            (i%capture_every==0 || i==steps) && capture()
            i==steps÷2 && checkpoint!(a,joinpath(directory,"middle.jls"))
        end
        checkpoint!(a,joinpath(directory,"final.jls"))
        open(io->JSON3.write(io,Dict("kind"=>"actual recorded trajectory","condition"=>condition,"frames"=>frames)),joinpath(directory,"sequence.json"),"w")
        final=arena_snapshot(a)["metrics"]
        report["arms"][condition]=Dict("config"=>cfg,"final"=>final,"peak_matrix_frame"=>peak_file,
            "peak_matrix_with_living_cells"=>peak_matrix,"max_abs_material_error"=>max_mass_error,"max_abs_chemical_error"=>max_energy_error)
        println(condition,": population=",final["population"]," births=",final["births"]," matrix=",sum(final["catalytic"]["enzyme_mass"]),
            " unlocked=",sum(final["catalytic"]["solubilized"])," peak_frame=",peak_file)
        open(io->JSON3.write(io,report),joinpath(output,"report.json"),"w")
    end
    report
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    length(ARGS)>=2 || error("usage: catalytic_habitat.jl PROFILE OUTPUT [STEPS=3000]")
    catalytic_study(ARGS[1],abspath(ARGS[2]);steps=length(ARGS)>2 ? parse(Int,ARGS[3]) : 3000)
end
