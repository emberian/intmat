#!/usr/bin/env julia
using Interior,TOML,JSON3,SHA

"""Actual 3D cells and chemistry, with physical lengths shared across resolutions.
Full raw fields and spherical-body geometry accompany each sampled state. No
extrusion or invented connection surface is exported as simulated anatomy.
"""
function constructive_volume(profile,output;resolution=24,steps=3000,capture_every=50,indexed=true,
        resume=nothing,wall_seconds=Inf,artifact_bytes=typemax(Int))
    resolution isa Integer && resolution>0 || throw(ArgumentError("resolution must be a positive integer"))
    steps isa Integer && steps>=0 || throw(ArgumentError("steps must be a nonnegative integer"))
    capture_every isa Integer && capture_every>0 || throw(ArgumentError("capture interval must be a positive integer"))
    wall_seconds>0 || throw(ArgumentError("wall time must be positive"))
    artifact_bytes isa Integer && artifact_bytes>0 || throw(ArgumentError("artifact budget must be positive"))
    ispath(output) && error("use a fresh output directory")
    cfg=TOML.parsefile(profile);cfg["grid"]["shape"]=[resolution,resolution,resolution]
    a=resume===nothing ? build_arena(cfg) : restore_arena(resume)
    # Arena owns this derived initial-energy ledger entry; it is not a recipe setting.
    saved_recipe=Dict(k=>v for (k,v) in a.config if k!="_initial_energy")
    saved_recipe==cfg || throw(ArgumentError("continuation profile must equal the saved configuration after resolution selection"))
    w=a.world;owner=only(filter(p->p isa ConstructiveEcology,a.processes))
    c=owner.cells;h=owner.habitat;m=only(filter(p->p isa SorptiveMatrix,a.processes));r=c.resources
    if resume===nothing
        indexed && Interior._cellular_enable_indexed_storage!(c)
    else
        Interior._cellular_storage_statistics(c).mode==(indexed ? "indexed" : "dense") ||
            throw(ArgumentError("continuation preserves the saved genome representation"))
    end
    volume=c.stencil.cell_volume;mkpath(output)
    frames=Any[];maxmass=0.;maxenergy=0.
    manifest=Dict{String,Any}("kind"=>"actual three-dimensional constructive ecology",
        "config"=>cfg,"loaded_sources"=>Interior.arena_sources(),"frames"=>frames,
        "study_source"=>Dict("path"=>abspath(@__FILE__),"sha256"=>bytes2hex(SHA.sha256(read(@__FILE__)))),
        "profile_source"=>Dict("path"=>abspath(profile),"sha256"=>bytes2hex(SHA.sha256(read(profile)))),
        "parent_checkpoint"=>(resume===nothing ? nothing : Dict("path"=>abspath(resume),"sha256"=>bytes2hex(SHA.sha256(read(resume))))),
        "trajectory_kind"=>(resume===nothing ? "fresh founders" : "strict checkpoint continuation"),
        "initial_step"=>a.steps,"initial_time"=>w.time,
        "requested_additional_steps"=>steps,"capture_every"=>capture_every,
        "wall_seconds_after_initial_export"=>(isfinite(wall_seconds) ? wall_seconds : nothing),"artifact_budget_bytes"=>artifact_bytes,
        "arena_dt"=>a.dt,
        "resource_diffusion"=>[Dict("species"=>string(r.names[k]),"coefficient"=>D,
            "explicit_limit"=>(D>0 ? Interior.diffusion_limit(w,D) : nothing),
            "internal_substeps"=>(D>0 ? max(1,ceil(Int,Float64(a.dt)/Interior.diffusion_limit(w,D))) : 0))
            for (k,D) in enumerate(r.diffusion)],
        "genome_storage"=>(indexed ? "indexed" : "dense"),
        "body_geometry"=>"physical positions and spherical radii; no extrusion",
        "field_units"=>"all chemistry is mass per physical volume; permeability dimensionless",
        "format"=>"little-endian Float32 raw arrays in Julia/Fortran axis order")
    function capture()
        directory="frame-"*lpad(a.steps,7,'0');dest=joinpath(output,directory);mkpath(dest)
        snapshot=arena_snapshot(a;max_agents=c.capacity);metrics=snapshot["metrics"]
        # Keep the existing viewer's honest central slice as well as full volumes.
        open(io->JSON3.write(io,snapshot),joinpath(dest,"snapshot.json"),"w")
        geometry=Dict("grid"=>snapshot["grid"],"agents"=>snapshot["agents"],"bonds"=>metrics["cell_bonds"],
            "time"=>w.time,"step"=>a.steps,"metrics"=>metrics)
        open(io->JSON3.write(io,geometry),joinpath(dest,"geometry.json"),"w")
        fields=(enzyme_a=h.enzymes[1],enzyme_b=h.enzymes[2],retained_a=m.bound[1],retained_b=m.bound[2],
            precursor_a=r.fields[1],precursor_b=r.fields[2],permeability=r.permeability)
        hashes=Dict{String,String}()
        for (name,field) in pairs(fields)
            data=name===:permeability ? Float32.(field) : Float32.(field./volume)
            path=joinpath(dest,string(name)*".f32");open(io->write(io,data),path,"w")
            hashes[string(name)]=bytes2hex(SHA.sha256(read(path)))
        end
        open(io->JSON3.write(io,Dict("shape"=>collect(w.grid.size),"extent"=>collect(w.grid.extent),
            "time"=>w.time,"step"=>a.steps,"dtype"=>"float32","order"=>"F","sha256"=>hashes)),joinpath(dest,"volume.json"),"w")
        maxmass=max(maxmass,abs(metrics["material_error"]));maxenergy=max(maxenergy,abs(metrics["chemical_energy_error"]))
        push!(frames,Dict("directory"=>directory,"step"=>a.steps,"time"=>w.time,"population"=>length(c.ids),
            "enzyme_mass"=>sum(sum,h.enzymes),"retained_mass"=>sum(sum,m.bound)))
        open(io->JSON3.write(io,manifest),joinpath(output,"run.json"),"w")
        open(joinpath(output,"metrics.jsonl"),"a") do io
            JSON3.write(io,Dict("step"=>a.steps,"time"=>w.time,"metrics"=>metrics));write(io,'\n')
        end
        println("resolution=",resolution," step=",a.steps," pop=",length(c.ids)," enzyme=",sum(sum,h.enzymes));flush(stdout)
    end
    capture();checkpoint!(a,joinpath(output,"initial.jls"));started=time();reason="steps complete"
    # Reserve several initial checkpoint sizes for the final restart state. The
    # trajectory stops before this reserve, never deleting earlier artifacts.
    reserve=4*filesize(joinpath(output,"initial.jls"))
    bytes_written()=sum(filesize(joinpath(root,file)) for (root,_,files) in walkdir(output) for file in files)
    bytes_written()+reserve<artifact_bytes || error("artifact budget cannot hold initial state and restart reserve")
    for i in 1:steps
        arena_step!(a)
        (i%capture_every==0 || i==steps) && capture()
        i==steps÷2 && checkpoint!(a,joinpath(output,"middle.jls"))
        if time()-started>=wall_seconds
            reason="wall budget";last(frames)["step"]==a.steps || capture();break
        elseif i%capture_every==0 && bytes_written()+reserve>=artifact_bytes
            reason="artifact budget";break
        end
    end
    checkpoint!(a,joinpath(output,"final.jls"))
    manifest["elapsed_seconds_including_exports"]=time()-started
    manifest["stop_reason"]=reason;manifest["artifact_bytes_before_manifest"]=bytes_written()
    manifest["max_abs_material_error"]=maxmass;manifest["max_abs_chemical_error"]=maxenergy
    manifest["final"]=arena_snapshot(a)["metrics"]
    open(io->JSON3.write(io,manifest),joinpath(output,"run.json"),"w")
    manifest
end

if abspath(PROGRAM_FILE)==(@__FILE__)
    2<=length(ARGS)<=8 || error("usage: constructive_volume.jl PROFILE OUTPUT [RESOLUTION=24] [STEPS=3000] [CAPTURE=50] [CHECKPOINT] [WALL_SECONDS] [ARTIFACT_BYTES]")
    constructive_volume(ARGS[1],abspath(ARGS[2]);resolution=length(ARGS)>2 ? parse(Int,ARGS[3]) : 24,
        steps=length(ARGS)>3 ? parse(Int,ARGS[4]) : 3000,capture_every=length(ARGS)>4 ? parse(Int,ARGS[5]) : 50,
        resume=length(ARGS)>5 ? ARGS[6] : nothing,wall_seconds=length(ARGS)>6 ? parse(Float64,ARGS[7]) : Inf,
        artifact_bytes=length(ARGS)>7 ? parse(Int,ARGS[8]) : typemax(Int))
end
