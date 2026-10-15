# Small bounded CPU laboratory for the separate fixed-voxel material law.
using Interior,TOML,SHA,Statistics
using Interior:JSON3,KernelAbstractions
include("lattice_ecology.jl")
using .LatticeEcologyStudy
const LE=LatticeEcologyStudy
const RUNNER_SOURCE=read(@__FILE__,String)
function lattice_topology(mask,dims)
    seen=falses(length(mask));sizes=Int[]
    for site in eachindex(mask)
        if mask[site]&&!seen[site]
            queue=[site];seen[site]=true;cursor=1
            while cursor<=length(queue)
                current=queue[cursor];cursor+=1
                for d in 1:3
                    if dims[d]>1
                        for direction in (-1,1)
                            other=LE.neighbor(current,dims,d,direction)
                            if mask[other]&&!seen[other];seen[other]=true;push!(queue,other);end
                        end
                    end
                end
            end
            push!(sizes,length(queue))
        end
    end
    (components=length(sizes),largest=isempty(sizes) ? 0 : maximum(sizes),sites=sum(sizes))
end
function lattice_run_limits(cfg;coordinated_volume=false)
    dims=Tuple(cfg["dimensions"])
    0<cfg["steps"]<=10000&&cfg["capture_every"]>0||error("bounded schedule required")
    frames=1+sum(step==get(cfg,"injury_step",0) ? 2 : (step%cfg["capture_every"]==0||step==cfg["steps"] ? 1 : 0) for step in 1:cfg["steps"])
    if coordinated_volume
        dims==(32,32,32)&&Tuple(cfg["extent"])==(32.,32.,32.)||error("coordinated opt-in is limited to spacing-one 32 cubed")
        length(cfg["seeds"])==1&&cfg["arms"]==["full"]&&cfg["steps"]<=3200&&frames<=42||error("coordinated volume schedule exceeds approved bounds")
        Threads.nthreads()<=2||error("coordinated volume permits at most two Julia threads")
    else
        prod(dims)<=4096||error("ordinary CPU pilot is bounded to 4096 voxels; approved 32 cubed run requires --coordinated-cpu-volume")
    end
    (coordinated_volume=coordinated_volume,frames=frames,native_frame_bytes=frames*prod(dims)*(15sizeof(Float32)+3sizeof(Int32)))
end
function lattice_run_main(recipe,out;coordinated_volume=false)
    cfg=TOML.parsefile(recipe);dims=Tuple(cfg["dimensions"])
    limits=lattice_run_limits(cfg;coordinated_volume)
    ispath(out)&&error("choose a new output directory");mkpath(out);cp(recipe,joinpath(out,"recipe.toml"))
    write(joinpath(out,"lattice_ecology.jl"),LE.MODEL_SOURCE);write(joinpath(out,"lattice_ecology_run.jl"),RUNNER_SOURCE)
    q=lattice_parameters(;NamedTuple(Symbol(k)=>v for (k,v) in get(cfg,"parameters",Dict()))...)
    warm=lattice_world(Grid(dims...;extent=Tuple(cfg["extent"])));warm_p=LatticeEcology(warm;parameters=q,loci=cfg["loci"],founder_genotypes=cfg["founder_genotypes"],occupancy=cfg["occupancy"])
    for _ in 1:32;step!(warm,warm_p,Float32(cfg["dt"]));end;KernelAbstractions.synchronize(warm.backend)
    results=Any[]
    for seed in cfg["seeds"],arm in cfg["arms"]
        patch=arm=="no_mutation" ? (mutation_probability=0f0,) : arm=="no_transport_barrier" ? (matrix_barrier=0f0,) : arm=="no_recurrence" ? (recurrence=0f0,) : arm=="full" ? (;) : error("unknown arm")
        parameters=lattice_parameters(;merge(q,patch)...)
        w=lattice_world(Grid(dims...;extent=Tuple(cfg["extent"]));seed)
        p=LatticeEcology(w;parameters,loci=cfg["loci"],occupancy=cfg["occupancy"],founder_genotypes=cfg["founder_genotypes"],
            reserve=cfg["reserve"],precursor=cfg["precursor"],fuel=cfg["fuel"],niches=cfg["spatial_niches"])
        Interior.check_dt(w,p,cfg["dt"]);dir=joinpath(out,"seed-$seed-$arm");mkpath(dir);mkdir(joinpath(dir,"frames"))
        initial=lattice_statistics(w,p);history=Any[];step_seconds=0.;injury=get(cfg,"injury_step",0);radius=Float32(get(cfg,"injury_radius",3.))
        region=[sum((Float32(LE.coords(i,dims)[d])-.5f0-Float32(dims[d])/2)^2 for d in 1:3 if dims[d]>1)<radius^2 for i in 1:prod(dims)]
        function capture(step,label)
            s=lattice_statistics(w,p);e=Array(p.environment);c=Array(p.cells);g=Array(p.genotype)
            live=g.>0;matrix=e[:,4].>.05f0
            sample=merge(s,(step=step,label=label,material_residual=s.material-initial.material,energy_residual=s.energy_accounted-initial.energy_accounted,
                injury_region_population=count(live[region]),injury_region_matrix=sum(Float64,e[region,4]),
                living_topology=lattice_topology(live,dims),matrix_topology=lattice_topology(matrix,dims)))
            sample.finite&&sample.minimum_physical>=0&&sample.capacity_excess<=2e-6&&sample.hidden_maxabs<=1.000001||error("invalid lattice physical state")
            abs(sample.material_residual)<1e-4initial.material&&abs(sample.energy_residual)<1e-4initial.energy_accounted||error("lattice budget drift exceeded pilot tolerance")
            sample.population==initial.population+sample.births-sample.deaths||error("lifecycle accounting failed")
            push!(history,sample)
            open(joinpath(dir,"history.jsonl"),"a") do io;JSON3.write(io,sample);write(io,'\n');end
            prefix=lpad(step,7,'0')*"-$label"
            floats=hcat(e,c,Array(p.hidden));ints=hcat(g,Array(p.lineage),Array(p.generation))
            open(io->write(io,floats),joinpath(dir,"frames",prefix*".f32"),"w")
            open(io->write(io,ints),joinpath(dir,"frames",prefix*".i32"),"w")
            open(joinpath(dir,"frames",prefix*".json"),"w") do io
                JSON3.write(io,(step=step,time=w.time,shape=dims,array_order="F",float_shape=size(floats),integer_shape=size(ints),
                    float_channels=(LE.ENV_CHANNELS...,:biomass,:stored_work,:age,:hidden1,:hidden2,:hidden3,:hidden4),integer_channels=(:genotype,:lineage,:generation)))
            end
        end
        capture(0,"initial")
        wall=@elapsed for step in 1:cfg["steps"]
            step_seconds+=@elapsed step!(w,p,Float32(cfg["dt"]))
            if step==injury
                capture(step,"pre-injury");lattice_damage!(w,p;radius);capture(step,"post-injury")
            elseif step%cfg["capture_every"]==0||step==cfg["steps"]
                capture(step,"sample")
            end
        end
        lattice_checkpoint(w,p,joinpath(dir,"final.jls"));record_state(w,p,joinpath(dir,"native-final");step=p.tick)
        result=(seed=seed,arm=arm,dimensions=dims,extent=w.grid.extent,parameters=parameters,initial=initial,final=last(history),history=history,
            run_limits=limits,step_seconds=step_seconds,wall_seconds_including_capture=wall,voxel_steps_per_second=prod(dims)*cfg["steps"]/step_seconds,
            policy="Four hidden channels; seven funded flux/event gates; nine local observations, a bias and four recurrent channels. Immutable finite-locus shared rule bank; bit mutation only at funded birth.",
            law_source_sha256=LE.MODEL_SHA,runner_sha256=bytes2hex(sha256(RUNNER_SOURCE)),shared_sources=Interior.arena_sources())
        push!(results,result);open(joinpath(dir,"result.json"),"w") do io;JSON3.write(io,result);end
        open(joinpath(out,"results.json"),"w") do io;JSON3.write(io,results);end
        println((seed=seed,arm=arm,population=result.final.population,births=result.final.births,deaths=result.final.deaths,
            mutations=result.final.mutations,matrix=result.final.matrix,material_residual=result.final.material_residual,
            energy_residual=result.final.energy_residual,voxel_steps_per_second=result.voxel_steps_per_second));flush(stdout)
    end
    println("READY ",out)
end
if abspath(PROGRAM_FILE)==(@__FILE__) || (isdefined(Main,:INTERIOR_RUNNING_STUDY)&&Main.INTERIOR_RUNNING_STUDY=="lattice_ecology_run")
    length(ARGS) in (2,3)||error("recipe.toml fresh-output-directory [--coordinated-cpu-volume]")
    length(ARGS)==3&&ARGS[3]!="--coordinated-cpu-volume"&&error("unknown option")
    lattice_run_main(ARGS[1],ARGS[2];coordinated_volume=length(ARGS)==3)
end
