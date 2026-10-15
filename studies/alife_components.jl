#!/usr/bin/env julia
# Exploratory ancestry-matched component transplants; never alters primary outputs.
using Interior, Random, JSON3, Serialization, TOML, SHA
include(joinpath(@__DIR__,"alife_adaptation.jl"))
function component_genotype(descendant,founder,indices,ancestors;controller=true,metabolism=true)
    (;weights=copy(controller ? descendant.weights[:,:,indices] : founder.weights[:,:,ancestors]),
      allocation=copy(metabolism ? descendant.allocation[:,indices] : founder.allocation[:,ancestors]),
      (k=>copy(metabolism ? getfield(descendant,k)[indices] : getfield(founder,k)[ancestors]) for k in TRAITS)...)
end
function components_main()
    protocol=TOML.parsefile(ARGS[2]);root=abspath(protocol["primary_directory"])
    cfg=TOML.parsefile(joinpath(root,"protocol.toml"))
    cfg["assay_capacity"]=protocol["assay_capacity"]; cfg["assay_time"]=protocol["assay_time"]
    out=abspath(ARGS[1]);ispath(out)&&!isempty(readdir(out))&&error("Use a fresh directory");mkpath(out)
    cp(ARGS[2],joinpath(out,"protocol.toml"));cp(@__FILE__,joinpath(out,"study-source.jl"))
    cp(joinpath(@__DIR__,"alife_adaptation.jl"),joinpath(out,"included-study-source.jl"))
    results=Any[]
    for seed in protocol["seeds"]
        primary=joinpath(root,"seed-$seed-$(protocol["selection_environment"])-mutation")
        bank=deserialize(joinpath(primary,"founders.jls"))
        saved=deserialize(joinpath(primary,"selection","terminal.jls"))
        saved.sources==Interior.arena_sources() || error("Physics hashes differ from primary study")
        selected=saved.state[4];desc=genotype(selected); ancestry=selected.controller.lineages
        @assert !isempty(ancestry) && all(i->1<=i<=size(bank.weights,3),ancestry)
        for rep in 1:protocol["position_replicates"]
            n=protocol["assay_per_group"];rng=MersenneTwister(seed+100000rep)
            indices=rand(rng,axes(desc.weights,3),n); ancestors=ancestry[indices]
            slots=randperm(rng,2n);ds=slots[1:n];fs=slots[n+1:end]
            baseline=component_genotype(desc,bank,indices,ancestors;controller=false,metabolism=false)
            for component in protocol["components"]
                test=component_genotype(desc,bank,indices,ancestors;
                    controller=component in ("evolved_both","ancestral_metabolism"),
                    metabolism=component in ("evolved_both","ancestral_controller"))
                state=fresh(cfg,seed+100000rep,protocol["assay_environment"],2n;mutation=false,plasticity=protocol["plasticity"])
                p=state[4];c=p.controller;c.capacity=cfg["assay_capacity"]
                install!(p,test,1:n,ds);install!(p,baseline,1:n,fs)
                c.lineages[ds].=1;c.lineages[fs].=2
                @assert all(iszero,c.hidden)&&all(iszero,c.plastic_marks)&&all(iszero,c.age)
                @assert all(==(.4f0),c.storage)&&all(==(1f0),c.energy)
                dir=joinpath(out,"seed-$seed-rep-$rep-$component")
                run=evolve!(state,cfg,cfg["assay_time"],dir;groups=true)
                counts=[count(==(i),c.lineages) for i in 1:2]
                r=(seed=seed,replicate=rep,component=component,descendant_indices=indices,
                    ancestral_founder_indices=ancestors,living=counts,births=run.group_births,deaths=run.group_deaths,
                    log_frequency_change=log((counts[1]+.5)/(counts[2]+.5)),both_extinct=sum(counts)==0,
                    test_extinct=counts[1]==0,ancestor_extinct=counts[2]==0,final=last(run.samples))
                push!(results,r);open(joinpath(dir,"result.json"),"w") do io;JSON3.write(io,r);end
                open(joinpath(out,"results.json"),"w") do io;JSON3.write(io,results);end
                println((seed=seed,replicate=rep,component=component,living=counts));flush(stdout)
            end
        end
    end
    open(joinpath(out,"provenance.json"),"w") do io
        JSON3.write(io,(sources=Interior.arena_sources(),julia=string(VERSION),script_sha256=bytes2hex(sha256(read(@__FILE__))),primary_directory=root))
    end
end
if abspath(PROGRAM_FILE)==(@__FILE__) || (isdefined(Main,:INTERIOR_RUNNING_STUDY) && Main.INTERIOR_RUNNING_STUDY=="alife_components")
    components_main()
end
