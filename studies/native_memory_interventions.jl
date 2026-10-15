#!/usr/bin/env julia
# Study-only interventions; every ordinary advance uses the native production law.
using Interior, Random, TOML, JSON3, SHA, Statistics
include(joinpath(@__DIR__,"../scripts/native_constructed_garden.jl"))

function assay_config(profile;seed=211,dt=.02)
    cfg=TOML.parsefile(profile)
    # Preserve physical density, spacing, concentrations, per-volume supply and
    # supply duration, in a 1^3 periodic subworld; not a production-world replica.
    cfg["seed"]=seed;cfg["grid"]["shape"]=[16,16,16];cfg["grid"]["extent"]=[1.,1.,1.]
    cfg["cellular"]["founders"]=4;cfg["cellular"]["capacity"]=256
    for section in ("resources","feedstock"),key in ("source_rates","source_budgets","sink_rates","sink_budgets")
        cfg[section][key]./=64
    end
    cfg["native"]["dt"]=dt
    cfg
end

function advance_assay!(w,p,n,dt)
    for _ in 1:n
        Interior.device_cellular_step!(w,p.cells,dt;anabolism=p.habitat,check_status=false)
    end
    Interior.device_cellular_status!(p.cells)
    w.time=Float64(p.cells.attempted_steps)*Float64(dt)
end

audit(p)=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
function intervention!(p,erase,pulse)
    before=audit(p);s=p.cells.state;h=p.habitat;a=h.arrays;e=p.cells.environment
    erase=="neural" && (fill!(s.hidden,0);fill!(s.adaptation,0))
    erase=="marks" && fill!(s.marks,0)
    # Release retained food once, retaining the scaffold and all subsequent laws.
    if erase=="retention"
        for k in 1:2,i in axes(a.bound,1)
            e.fields[i,k]=eltype(e.fields)(Float64(e.fields[i,k])+Float64(a.bound[i,k]));a.bound[i,k]=0
        end
    end
    # Accessible-food pulse: convert 10% of A/B raw inventory in x<1/2 to
    # same-potential soluble food. This is redistribution, never rescue funding.
    if pulse
        dims=Int.(p.cells.geometry.dims)
        for k in 1:2,i in axes(e.fields,1)
            mod(i-1,dims[1])<dims[1]÷2 || continue
            old=h.feedstock.arrays.fields[i,k];new=typeof(old)(.9*Float64(old))
            h.feedstock.arrays.fields[i,k]=new
            e.fields[i,k]=eltype(e.fields)(Float64(e.fields[i,k])+Float64(old)-Float64(new))
        end
    end
    after=audit(p)
    abs(after.material-before.material)<1e-6 || error("intervention material rounding excessive")
    abs(after.potential-before.potential)<4e-6 || error("intervention potential rounding excessive")
    (;material_jump=after.material-before.material,potential_jump=after.potential-before.potential)
end

function observe(w,p,origin)
    s=p.cells.state;a=p.habitat.arrays;e=p.cells.environment
    live=findall(==(1),s.active);x=audit(p)
    cells=[(;id=Int(s.ids[i]),parent=Int(s.parents[i]),lineage=Int(s.lineages[i]),
        hidden=collect(s.hidden[i,:]),adaptation=collect(s.adaptation[i,:]),actions=collect(s.actions[i,:]),
        mass=s.mass[i],energy=s.energy[i],metabolites=collect(s.metabolites[i,:]),
        secretion=collect(a.last_secretion[i,:])) for i in live]
    (;time=w.time,population=length(live),births=Int(e.metadata[3]),deaths=Int(e.metadata[4]),cells,
        free=sum(Float64,e.fields[:,1:2]),raw=sum(Float64,p.habitat.feedstock.arrays.fields),
        bound=sum(Float64,a.bound),catalyst=sum(Float64,reinterpret(Float32,a.enzyme_bits)),
        permeability_mean=mean(a.permeability),mark_rms=isempty(live) ? 0. : sqrt(mean(abs2,s.marks[live,:,:])),
        material_residual=x.material-origin.material-x.supplied+origin.supplied+x.removed-origin.removed,
        potential_residual=x.potential-origin.potential+x.work-origin.work+x.heat-origin.heat-x.supplied_potential+origin.supplied_potential+x.removed_potential-origin.removed_potential,
        finite=x.finite)
end

function spatial_frame(w,p)
    s=p.cells.state;a=p.habitat.arrays;live=findall(==(1),s.active)
    volume=Float64(p.cells.geometry.volume)
    (;time=w.time,shape=Int.(p.cells.geometry.dims),extent=p.cells.geometry.extent,
        field_order=["free_a","free_b","raw_a","raw_b","enzyme_a","enzyme_b","bound_a","bound_b"],
        field_units="concentration except dimensionless scaffold permeability; Julia column-major voxel ordering",
        fields=[collect(Float64.(v)./volume) for v in (eachcol(p.cells.environment.fields[:,1:2])...,eachcol(p.habitat.feedstock.arrays.fields[:,1:2])...,eachcol(reinterpret(Float32,a.enzyme_bits))...,eachcol(a.bound)...)],
        permeability=collect(a.permeability),
        cells=[(;id=Int(s.ids[i]),position=collect(s.position[i,:]),radius=s.radii[i],mass=s.mass[i],energy=s.energy[i],hidden=collect(s.hidden[i,:])) for i in live])
end

function run_assay(profile,out;warm=40.,horizon=64.,dt=.02,cadence=.2,seed=211)
    all(isfinite,(warm,horizon,dt,cadence)) && warm>=0 && horizon>0 && dt>0 && cadence>=dt || error("invalid assay durations")
    for duration in (warm,horizon,cadence)
        isapprox(duration/dt,round(duration/dt);atol=1e-8) || error("durations must be integer multiples of dt")
    end
    ispath(out) && error("fresh output required");mkpath(out)
    driver_text=read(@__FILE__,String);driver_hash=bytes2hex(sha256(driver_text))
    write(joinpath(out,"executed-driver.jl"),driver_text)
    factory_hash=bytes2hex(sha256(read(joinpath(@__DIR__,"../scripts/native_constructed_garden.jl"))))
    cfg=assay_config(profile;seed,dt);delta=Float32(dt)
    w,p,funding=native_constructed_garden(CPU(),cfg)
    advance_assay!(w,p,round(Int,warm/dt),delta)
    checkpoint=Interior.device_habitat_checkpoint(w,p)
    expected=Interior.device_habitat_snapshot(p)
    open(io->JSON3.write(io,spatial_frame(w,p)),joinpath(out,"before.json"),"w")
    arms=Dict{String,Any}();start=time()
    for erase in ("none","neural","marks","retention"),pulse in (false,true)
        name=erase*(pulse ? "_pulse" : "_sham")
        wr=World(w.grid,CPU(),NamedTuple(),NamedTuple(),0.,w.seed)
        branch=Interior.device_habitat_restore(wr,checkpoint)
        Interior.device_habitat_snapshot(branch)==expected || error("branch differs before intervention")
        jump=intervention!(branch,erase,pulse);origin=audit(branch)
        open(io->JSON3.write(io,spatial_frame(wr,branch)),joinpath(out,name*"-initial-field.json"),"w")
        trace=[observe(wr,branch,origin)];stride=round(Int,cadence/dt);steps=round(Int,horizon/dt)
        for i in 1:steps
            advance_assay!(wr,branch,1,delta)
            (i%stride==0 || i==steps) && push!(trace,observe(wr,branch,origin))
        end
        open(io->JSON3.write(io,spatial_frame(wr,branch)),joinpath(out,name*"-final-field.json"),"w")
        all(t->t.finite,trace) || error("nonfinite branch")
        arms[name]=(;erase,pulse,jump,trace)
        open(io->JSON3.write(io,arms[name]),joinpath(out,name*".json"),"w")
        println(name," population=",last(trace).population," bound=",last(trace).bound," elapsed=",round(time()-start;digits=1));flush(stdout)
    end
    # Same-ID difference in pulse responses between intact/erased memory, not
    # population-average variance or a proxy for learning.
    interactions=Dict{String,Any}()
    for erase in ("neural","marks","retention")
        rows=Any[]
        for j in eachindex(arms["none_sham"].trace)
            ts=[arms[k].trace[j] for k in ("none_sham","none_pulse",erase*"_sham",erase*"_pulse")]
            maps=[Dict(c.id=>c for c in t.cells) for t in ts]
            ids=sort!(collect(intersect((Set(keys(m)) for m in maps)...)))
            rms(field)=isempty(ids) ? nothing : sqrt(mean((Float64(getproperty(maps[4][id],field)[k])-Float64(getproperty(maps[3][id],field)[k])-Float64(getproperty(maps[2][id],field)[k])+Float64(getproperty(maps[1][id],field)[k]))^2 for id in ids for k in eachindex(getproperty(maps[1][id],field))))
            push!(rows,(;time=ts[1].time,common_ids=ids,hidden_interaction_rms=rms(:hidden),action_interaction_rms=rms(:actions),
                bound_interaction=ts[4].bound-ts[3].bound-ts[2].bound+ts[1].bound))
        end
        interactions[erase]=rows
    end
    result=(;schema="native-memory-interventions-v1",seed,config=cfg,warm,horizon,dt=Float64(delta),cadence,
        source_pin=checkpoint.source,cellular_source_pin=checkpoint.cells.source,driver_sha256=driver_hash,factory_sha256=factory_hash,
        runtime=(julia=string(VERSION),resource_inventory=string(eltype(p.cells.environment.fields)),feedstock_inventory=string(eltype(p.habitat.feedstock.arrays.fields)),battery=string(eltype(p.cells.state.energy))),
        profile_sha256=bytes2hex(sha256(read(profile))),scope="fresh CPU native 16^3 periodic subworld, 4 founders; matched common-state causal interventions; no adaptation/criticality claim",
        elapsed_seconds=time()-start,interactions,
        final=Dict(k=>(;population=last(v.trace).population,material_residual=last(v.trace).material_residual,potential_residual=last(v.trace).potential_residual,jump=v.jump) for (k,v) in arms))
    open(io->JSON3.pretty(io,result),joinpath(out,"report.json"),"w")
    result
end
if abspath(PROGRAM_FILE)==(@__FILE__) || (isdefined(Main,:INTERIOR_RUNNING_STUDY) && Main.INTERIOR_RUNNING_STUDY=="native_memory_interventions")
    length(ARGS)>=2 || error("PROFILE FRESH_OUTPUT [WARM=40] [HORIZON=64] [DT=.02]")
    run_assay(ARGS[1],ARGS[2];warm=length(ARGS)>2 ? parse(Float64,ARGS[3]) : 40.,horizon=length(ARGS)>3 ? parse(Float64,ARGS[4]) : 64.,dt=length(ARGS)>4 ? parse(Float64,ARGS[5]) : .02)
end
