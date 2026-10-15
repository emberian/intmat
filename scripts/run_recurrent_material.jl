#!/usr/bin/env julia
using Interior,Statistics
using Interior: JSON3
const loaded_model_source=read(joinpath(@__DIR__,"..","src","recurrent_material.jl"),String)
const loaded_core_source=read(joinpath(@__DIR__,"..","src","core.jl"),String)
isdefined(Interior,:RecurrentMaterial) || Base.include(Interior,joinpath(@__DIR__,"..","src","recurrent_material.jl"))
using Interior: RecurrentMaterial
include(joinpath(@__DIR__,"recurrent_runtime.jl"))
const recurrent_profiles=(("braided",(coupling=.12,feedback=1.2,writing=.25,erosion=.06,contrast=.7)),
    ("slow-memory",(coupling=.28,feedback=2.,writing=.5,erosion=.035,contrast=1.1)),
    ("threshold-feedback-off",(coupling=.12,feedback=0.,writing=.25,erosion=.06,contrast=.7)))
function parse_recurrent_args(args)
    positional=String[]; flags=Dict{String,String}();i=1
    allowed=("backend","profile","captures","capture-mode","states","budget-mb","resident-mb","batch")
    while i<=length(args)
        if startswith(args[i],"--")
            flag=args[i][3:end]
            flag in allowed || throw(ArgumentError("unknown option --$flag"))
            i<length(args) || throw(ArgumentError("--$flag requires a value"))
            haskey(flags,flag) && throw(ArgumentError("duplicate --$flag"))
            flags[flag]=args[i+1];i+=2
        else
            push!(positional,args[i]);i+=1
        end
    end
    length(positional)<=3 || throw(ArgumentError("side steps output-directory [options]"))
    side=length(positional)>0 ? parse(Int,positional[1]) : 64
    steps=length(positional)>1 ? parse(Int,positional[2]) : 3000
    root=length(positional)>2 ? positional[3] : "outputs/recurrent-material"
    backend=get(flags,"backend","metal");profile=get(flags,"profile",get(ENV,"RECURRENT_PROFILE","all"))
    captures=parse(Int,get(flags,"captures","12"));states=parse(Int,get(flags,"states","1"))
    mode=get(flags,"capture-mode","slices");budget=parse(Int,get(flags,"budget-mb","1536"))*1_000_000
    resident_budget=parse(Int,get(flags,"resident-mb","8192"))*1_000_000
    # Three world fields, two world scratch buffers, two eight-channel latent buffers.
    resident_bytes=big(side)^3*84
    resident_budget>0 && resident_bytes<=resident_budget || throw(ArgumentError("resident arrays need $(Float64(resident_bytes)/1e6) MB; --resident-mb budget is $(resident_budget/1e6)"))
    batch=parse(Int,get(flags,"batch","64"))
    side>0 && steps>0 && batch>0 && budget>0 || throw(ArgumentError("side, steps, batch and budget must be positive"))
    backend in ("cpu","metal","amd") || throw(ArgumentError("backend must be cpu, metal or amd"))
    profile in ("all",first.(recurrent_profiles)...) || throw(ArgumentError("unknown profile"))
    mode in ("slices","volumes","none") || throw(ArgumentError("capture mode must be slices, volumes or none"))
    captures>=0 && states>=0 || throw(ArgumentError("capture counts must be nonnegative"))
    selected=filter(p->profile=="all" || first(p)==profile,recurrent_profiles)
    frames=recurrent_capture_steps(steps,captures);checkpoints=recurrent_capture_steps(steps,states)
    # Raw payload plus a conservative metadata allowance, over the entire run.
    perframe=mode=="volumes" ? big(side)^3*8 : mode=="slices" ? big(side)^2*24 : 0
    estimated=length(selected)*(length(frames)*perframe+length(checkpoints)*big(side)^3*44+1_000_000)
    estimated<=budget || throw(ArgumentError("capture plan needs $(round(Float64(estimated)/1e6;digits=1)) MB, budget $(budget/1e6) MB; reduce --captures/--states, select one --profile, or raise --budget-mb"))
    (side=side,steps=steps,root=root,backend=backend,profile=profile,captures=captures,states=states,
        capture_mode=mode,budget_bytes=budget,resident_array_bytes=Int(resident_bytes),resident_budget_bytes=resident_budget,batch=batch,estimated_bytes=Int(estimated),selected=selected,
        frames=frames,checkpoints=checkpoints)
end
function write_recurrent_capture(dir,field,array,k,mode)
    if mode=="volumes"
        open(joinpath(dir,"$field-$k.f32"),"w") do io;write(io,array);end
    elseif mode=="slices"
        for axis in 1:3
            plane=Array(selectdim(array,axis,cld(size(array,axis),2)))
            open(joinpath(dir,"$field-$(axis)-$k.f32"),"w") do io;write(io,plane);end
        end
    end
end
function warm_recurrent_material(backend)
    w=recurrent_world(Grid(4,4,4;extent=(4,4,4)),backend)
    p=RecurrentMaterial(w)
    advance!(w,p,.05)
    Interior.KernelAbstractions.synchronize(backend)
    nothing
end
recurrent_world(grid,backend)=World(grid;backend,seed=31,
    fields=(:nutrient,:deposit,:activity),scratch_fields=(:nutrient,:deposit))
function run_recurrent_material(options)
    o=options;root=o.root;side=o.side;steps=o.steps
    ispath(root) && throw(ArgumentError("output directory exists; use a fresh destination"))
    backend=recurrent_backend(o.backend)
    mkpath(dirname(abspath(root)));mkdir(root) # Atomic exclusive output claim.
    write(joinpath(root,"loaded-model-source.jl"),loaded_model_source)
    write(joinpath(root,"loaded-core-source.jl"),loaded_core_source)
    cp(@__FILE__,joinpath(root,"runner-source.jl"))
    cp(joinpath(@__DIR__,"recurrent_runtime.jl"),joinpath(root,"runtime-source.jl"))
    open(joinpath(root,"run-config.json"),"w") do io
        JSON3.pretty(io,(backend=o.backend,side=side,steps=steps,profile=o.profile,capture_mode=o.capture_mode,
            capture_steps=sort!(collect(o.frames)),state_steps=sort!(collect(o.checkpoints)),
            estimated_capture_bytes=o.estimated_bytes,capture_budget_bytes=o.budget_bytes,resident_array_bytes=o.resident_array_bytes,
            resident_budget_bytes=o.resident_budget_bytes,batch=o.batch,
            model_source_sha256=bytes2hex(Interior.SHA.sha256(loaded_model_source)),
            core_source_sha256=bytes2hex(Interior.SHA.sha256(loaded_core_source))))
    end
    # Small-grid compilation avoids retaining two full-size worlds during reset.
    warm_recurrent_material(backend)
    for (name,kwargs) in o.selected
        w=recurrent_world(Grid(side,side,side;extent=Float32.((side,side,side))),backend)
        p=RecurrentMaterial(w;kwargs...);dir=joinpath(root,name);mkdir(dir)
        before=Float64(side)^3;samples=Any[]
        open(joinpath(dir,"volume.json"),"w") do io
            JSON3.pretty(io,(shape=w.grid.size,extent=w.grid.extent,dtype="float32",order="Julia-column-major",dt=.05,
                capture_mode=o.capture_mode,capture_steps=sort!(collect(o.frames)),
                slice_axes="axis 1=yz, 2=xz, 3=xy; middle cell cld(side,2)",latent_channels=8))
        end
        started=time()
        for k in 1:steps
            advance!(w,p,.05);w.time+=Float64(Float32(.05))
            # Bound queued launches even across arbitrarily long sparse captures.
            k%o.batch==0 && Interior.KernelAbstractions.synchronize(backend)
            if k in o.frames || k==steps
                Interior.KernelAbstractions.synchronize(backend)
                activity=Array(w.fields.activity);solid=Array(w.fields.deposit);latent=Array(p.state)
                push!(samples,(step=k,time=w.time,activity_std=std(activity),deposit_mean=mean(solid),deposit_std=std(solid),
                    latent_std=std(latent),material_error=sum(Float64,Array(w.fields.nutrient))+sum(Float64,solid)-before,
                    finite=all(isfinite,latent)))
                if k in o.frames
                    write_recurrent_capture(dir,"activity",activity,k,o.capture_mode)
                    write_recurrent_capture(dir,"deposit",solid,k,o.capture_mode)
                end
            end
            k in o.checkpoints && record_state(w,p,joinpath(dir,"state-$k");step=k)
        end
        elapsed=time()-started
        result=(side=side,steps=steps,dt=.05,backend=o.backend,seconds=elapsed,
            voxel_steps_per_second=side^3*steps/elapsed,settings=kwargs,samples=samples)
        open(joinpath(dir,"metrics.json"),"w") do io;JSON3.pretty(io,result);end
        println(name," ",side,"³ ",round(result.voxel_steps_per_second/1e6;digits=2)," M voxelsteps/s final ",last(samples));flush(stdout)
        # End the previous profile's ownership before constructing the next one.
        # In particular p retains both device latent/resource buffer pairs.
        Interior.KernelAbstractions.synchronize(backend)
        w=nothing;p=nothing;activity=nothing;solid=nothing;latent=nothing
        samples=nothing;result=nothing
        GC.gc(true)
    end
    root
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    if "--help" in ARGS
        println("side steps output-directory [--backend cpu|metal|amd] [--profile all|braided|slow-memory|threshold-feedback-off] [--captures 12] [--capture-mode slices|volumes|none] [--states 1] [--budget-mb 1536] [--resident-mb 8192] [--batch 64]")
    else
        options=parse_recurrent_args(ARGS)
        prepare_recurrent_backend!(options.backend)
        Base.invokelatest(run_recurrent_material,options)
    end
end
