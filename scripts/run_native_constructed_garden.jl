#!/usr/bin/env julia
include("device_cellular_setup.jl")
include("native_constructed_garden.jl")
include("device_habitat_capture.jl")
using Dates
plain_config(x::JSON3.Object)=Dict(string(k)=>plain_config(v) for (k,v) in pairs(x))
plain_config(x::JSON3.Array)=plain_config.(collect(x))
plain_config(x)=x
hashfile(path)=bytes2hex(Interior.SHA.sha256(read(path)))
function garden_frame!(path,w,p,buffer,side)
    c=p.cells;h=p.habitat;Interior.device_cellular_status!(c)
    habitat_coarse!(c.backend,128)(buffer,c.environment,h.arrays,h.feedstock.arrays,c.geometry,Int32(side);ndrange=side^3)
    fields=Array(buffer);active=findall(==(1),Array(c.state.active))
    select(a)=Array(a)[active,:]
    frame=(schema="interior-native-constructive-frame",version=1,time=w.time,step=c.attempted_steps,
        completed_steps=c.attempted_steps,trace_ids=Array(c.state.ids)[active],
        energy_storage_type=string(eltype(c.state.energy)),
        extent=c.geometry.extent,shape=(side,side,side),native_shape=c.geometry.dims,
        field_names=(:free_a,:free_b,:raw_a,:raw_b,:enzyme_a,:enzyme_b,:bound_a,:bound_b,:permeability),
        field_units="concentration except dimensionless permeability; exact native voxel values or explicit block means",
        fields=reshape(fields,side,side,side,9),ids=Array(c.state.ids)[active],parents=Array(c.state.parents)[active],
        lineages=Array(c.state.lineages)[active],generations=Array(c.state.generations)[active],
        position=select(c.state.position),radii=Array(c.state.radii)[active],mass=Array(c.state.mass)[active],
        energy=Array(c.state.energy)[active],metabolites=select(c.state.metabolites),bond_ids=select(c.state.bond_ids),
        hidden=select(c.state.hidden),adaptation=select(c.state.adaptation),actions=select(c.state.actions),
        catalytic_expression=select(h.arrays.expression),last_step_secretion=select(h.arrays.last_secretion))
    serialize(path,frame);(file=basename(path),time=w.time,step=c.attempted_steps,side=side,
        population=length(active),bytes=filesize(path),sha256=hashfile(path))
end
function main(profile,output;resume=nothing)
    ispath(output) && error("fresh output required")
    cfg=TOML.parsefile(profile);run=cfg["run"];q=cfg["native"]
    all(k->run[k]>0,("target_time","wall_seconds","metrics_seconds","checkpoint_seconds","capture_time","native_capture_time","artifact_budget_bytes")) || error("positive run controls required")
    dt=Float32(q["dt"]);batch=Int(q["check_every"]);batch>0 || error("positive check batch required")
    coarse=Int(run["coarse_side"]);native=Int(cfg["grid"]["shape"][1]);coarse>0&&native%coarse==0 || error("coarse side must divide native side")
    deadline=min(time()+Float64(run["wall_seconds"]),Dates.datetime2unix(DateTime(run["deadline_utc"])))
    deadline>time()+30 || error("deadline too close or elapsed")
    maxsteps=parse(Int,get(ENV,"NATIVE_GARDEN_MAX_STEPS",string(typemax(Int))))
    maxsteps>=0 || error("nonnegative optional step bound required")
    mkpath(output);mkpath(joinpath(output,"frames","coarse"));mkpath(joinpath(output,"frames","native"));mkdir(joinpath(output,"checkpoints"))
    started=time();construction=@elapsed begin
        if resume===nothing
            w,p,funding=native_constructed_garden(backend,cfg)
        else
            payload=deserialize(resume)
            grid=Grid(Tuple(Int.(cfg["grid"]["shape"]))...;extent=Tuple(Float64.(cfg["grid"]["extent"])))
            w=World(grid,backend,NamedTuple(),NamedTuple(),0.,Int(cfg["seed"]))
            p=Interior.device_habitat_restore(w,payload)
            prior=JSON3.read(read(joinpath(dirname(dirname(resume)),"run.json"),String))
            plain_config(JSON3.read(JSON3.write(cfg)))==plain_config(prior.config) || error("resume config must match the original run exactly")
            funding=prior.funding
        end
    end
    initial=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
    dt<=Interior.stability_limit(w,p) || error("dt exceeds coupled limit")
    bytes=resident_bytes(p);budget=Int(run["artifact_budget_bytes"])
    5bytes+2native^3*9*4+2_000_000<budget || error("budget cannot retain initial/final/two rolling checkpoints and capture reserve")
    buffers=Dict(side=>Interior.KernelAbstractions.zeros(backend,Float32,side^3,9) for side in unique((coarse,native)))
    frames=Any[];archives=Dict{String,Any}();audit_failure=false;initialstep=p.cells.attempted_steps;initialtime=w.time
    manifest=Dict{String,Any}("kind"=>resume===nothing ? "fresh native constructive ecology" : "strict native checkpoint continuation",
        "complete"=>false,"config"=>cfg,"funding"=>funding,"source"=>Interior.DEVICE_CELLULAR_SOURCE_PIN,
        "habitat_source"=>Interior.DEVICE_HABITAT_SOURCE_PIN,"runtime"=>device_cellular_runtime(),
        "energy_storage_type"=>string(eltype(p.cells.state.energy)),
        "loaded_project_sources"=>Interior.arena_sources(),"source_commit"=>get(ENV,"DEVICE_SOURCE_COMMIT","unspecified"),
        "driver_sha256"=>hashfile(@__FILE__),"fixture_sha256"=>hashfile(joinpath(@__DIR__,"native_constructed_garden.jl")),
        "capture_helper_sha256"=>hashfile(joinpath(@__DIR__,"device_habitat_capture.jl")),"profile_sha256"=>hashfile(profile),
        "parent_checkpoint"=>resume===nothing ? nothing : Dict("sha256"=>hashfile(resume),"file"=>basename(resume)),
        "initial"=>initial,"initial_step"=>initialstep,"initial_time"=>initialtime,"resident_array_bytes"=>bytes,
        "deadline_unix"=>deadline,"construction_seconds"=>construction,"frames"=>frames,"checkpoints"=>archives,
        "limits"=>"explicit native Float32 fields, geometry and controller with separately declared cellular energy precision; fixed checked contact substeps and native junction/bond rules; no CPU MatrixDrag; snapshots are actual states, no interpolation; sampled per-cell secretion is last physical step")
    write_manifest()=open(io->JSON3.pretty(io,manifest),joinpath(output,"run.json"),"w")
    diskbytes()=sum(filesize(joinpath(root,name)) for (root,_,files) in walkdir(output) for name in files;init=0)
    function accounting(audit)
        material=audit.material-initial.material-audit.supplied+initial.supplied+audit.removed-initial.removed
        chemical=audit.potential-initial.potential+audit.work-initial.work+audit.heat-initial.heat-audit.supplied_potential+initial.supplied_potential+audit.removed_potential-initial.removed_potential
        rtol=Float64(run["ledger_relative_tolerance"])
        material_tolerance=max(1e-4,rtol*(initial.material+max(0.,audit.supplied-initial.supplied)))
        chemical_tolerance=max(1e-4,rtol*(initial.potential+max(0.,audit.supplied_potential-initial.supplied_potential)))
        (material_error=material,chemical_error=chemical,material_tolerance=material_tolerance,
            chemical_tolerance=chemical_tolerance,passed=audit.finite&&abs(material)<=material_tolerance&&abs(chemical)<=chemical_tolerance)
    end
    function checkpoint(label)
        path=joinpath(output,"checkpoints",label*".jls");temporary=path*".pending"
        payload=Interior.device_habitat_checkpoint(w,p);serialize(temporary,payload)
        mv(temporary,path;force=true) # only this run's documented rolling slots are replaced
        observed=Interior.device_habitat_audit((cells=payload.cells.snapshot,arrays=payload.arrays,feedstock=payload.feedstock,
            parameters=payload.parameters,feed_parameters=payload.feed_parameters))
        ledger=accounting(observed);audit_failure|=!ledger.passed
        entry=(file=relpath(path,output),step=p.cells.attempted_steps,time=w.time,bytes=filesize(path),sha256=hashfile(path),audit=observed,accounting=ledger)
        archives[label]=entry;manifest["latest_complete_checkpoint"]=entry
        open(joinpath(output,"checkpoint-audits.jsonl"),"a") do io;JSON3.write(io,entry);write(io,'\n');end
        write_manifest()
        payload
    end
    function capture(side)
        any(f->f.step==p.cells.attempted_steps&&f.side==side,frames) && return
        file="frame-"*lpad(p.cells.attempted_steps,9,'0')*".jls"
        path=joinpath(output,"frames",side==native ? "native" : "coarse",file)
        entry=garden_frame!(path,w,p,buffers[side],side)
        push!(frames,merge(entry,(file=relpath(path,output),)));write_manifest()
    end
    # Strict restore is checked on an independent CPU carrier to avoid retaining
    # a second device world throughout the production trajectory.
    payload=checkpoint("initial")
    cw=World(w.grid,CPU(),NamedTuple(),NamedTuple(),0.,w.seed);restored=Interior.device_habitat_restore(cw,payload)
    rs=Interior.device_habitat_snapshot(restored)
    rs.cells==payload.cells.snapshot&&rs.arrays==payload.arrays&&rs.feedstock==payload.feedstock || error("strict initial restore differs")
    manifest["initial_checkpoint_strict_restore_exact"]=true
    payload=nothing;restored=nothing;rs=nothing;GC.gc()
    capture(native);coarse!=native&&capture(coarse)
    nextcapture=w.time+run["capture_time"];nextnative=w.time+run["native_capture_time"]
    lastmetrics=time();lastcheckpoint=time();lastbudgetcheck=0.;slot=0;reason="target time";warmbatches=Any[]
    manifest["warm_batches"]=warmbatches
    while w.time<Float64(run["target_time"]) && p.cells.attempted_steps-initialstep<maxsteps
        if audit_failure;reason="completed-state accounting guard";break;end
        if time()>=deadline;reason="wall deadline";break;end
        if time()-lastbudgetcheck>=30
            if diskbytes()+2bytes+2native^3*9*4>budget;reason="artifact budget reserve";break;end
            lastbudgetcheck=time()
        end
        count=min(batch,maxsteps-(p.cells.attempted_steps-initialstep),ceil(Int,(Float64(run["target_time"])-w.time)/Float64(dt)))
        batch_seconds=@elapsed begin
            for _ in 1:count
                Interior.device_cellular_step!(w,p.cells,dt;anabolism=p.habitat,check_status=false)
            end
            Interior.device_cellular_status!(p.cells)
        end
        p.cells.attempted_steps-initialstep>=128&&length(warmbatches)<64&&push!(warmbatches,(steps=count,seconds=batch_seconds))
        w.time=initialtime+(p.cells.attempted_steps-initialstep)*Float64(dt)
        if w.time>=nextcapture;capture(coarse);nextcapture+=run["capture_time"];end
        if w.time>=nextnative;coarse!=native&&capture(native);nextnative+=run["native_capture_time"];end
        if time()-lastmetrics>=run["metrics_seconds"]
            stats=Interior.device_cellular_statistics(p.cells)
            open(joinpath(output,"metrics.jsonl"),"a") do io
                JSON3.write(io,(wall_unix=time(),time=w.time,cells=stats));write(io,'\n')
            end
            stats.blocked_births==0 || error("population guard reached; do not interpret capped selection")
            lastmetrics=time()
        end
        if time()-lastcheckpoint>=run["checkpoint_seconds"]
            slot=1-slot;checkpoint("rolling-"*string(slot));lastcheckpoint=time()
        end
    end
    p.cells.attempted_steps-initialstep>=maxsteps&&(reason="explicit step bound")
    capture(native);checkpoint("final")
    final=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
    manifest["final"]=final;manifest["final_accounting"]=accounting(final);manifest["material_error"]=final.material-initial.material-final.supplied+initial.supplied+final.removed-initial.removed
    manifest["chemical_error"]=final.potential-initial.potential+final.work-initial.work+final.heat-initial.heat-final.supplied_potential+initial.supplied_potential+final.removed_potential-initial.removed_potential
    manifest["elapsed_seconds"]=time()-started;manifest["completed_steps"]=p.cells.attempted_steps-initialstep
    audit_failure&&(reason="completed-state accounting guard");manifest["stop_reason"]=reason;manifest["complete"]=true;manifest["stored_bytes"]=diskbytes();write_manifest()
    println((stop=reason,time=w.time,population=final.population,steps=manifest["completed_steps"],seconds=manifest["elapsed_seconds"]));flush(stdout)
end
length(ARGS)>=2 || error("PROFILE FRESH_OUTPUT [PARENT_CHECKPOINT] required")
try
    main(ARGS[1],ARGS[2];resume=length(ARGS)>2 ? ARGS[3] : nothing)
catch err
    isdir(ARGS[2])&&open(io->JSON3.pretty(io,(error=sprint(showerror,err),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,
        habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,partial_state_resumable=false)),joinpath(ARGS[2],"failure.json"),"w")
    rethrow()
end
