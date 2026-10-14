#!/usr/bin/env julia
# Bounded causal study; all membrane and signaling actions remain funded.
include(joinpath(@__DIR__,"arena_experiments.jl"))
function boundary_metrics(phi,previous,h)
    crossings=0
    for I in CartesianIndices(phi),d in 1:2
        J=CartesianIndex(ntuple(k->k==d ? mod1(I[k]+1,size(phi,k)) : I[k],3))
        crossings+=(phi[I]>0)!=(phi[J]>0)
    end
    negative=morphology_metrics(-phi)
    (interface_length=crossings*h,phase_sign_changes=sum((phi.>0).!=(previous.>0)),
     negative_components=negative.components,negative_largest_component_share=negative.largest_component_share)
end
function boundary_fixture(name,condition,params)
    dims=(48,48,1); seed=11
    w=World(Grid(dims...;extent=(24,24,1));seed)
    phase=ImplicitPhaseSeparation(w;mean=-.12,noise=0,mineral_pinning=0,nutrient_diffusion=0,interface_activity=0,activity_relaxation=0)
    copyto!(phase.phi,mode_initial_state(dims;seed)); Interior.phase_observe!(w,phase)
    medium=ActiveMedium(w,phase;remodeling_gain=params.gain,signal_bias=params.bias,membrane_pinning=params.pinning)
    founders=condition=="no_inhabitants" ? 0 : 96
    rng=MersenneTwister(seed+1000)
    positions=rand(rng,Float32,3,founders).*reshape(Float32[48,48,0],3,1)
    feedback=condition=="full"
    cells=NeuralEcology(positions;medium,seed,capacity=700,lifespan=80,uptake_rate=.7,
        secretion_rate=feedback ? params.secretion : 0,remodeling_rate=feedback ? params.rate : 0,
        signal_rate=feedback ? .1 : 0,catabolic_fraction=.5,energy_yield=2)
    flux=Chemostat(;source_rate=.023,sink_rate=.01)
    arena=Arena(w,(flux,cells,medium,phase);dt=.1,config=Dict("initial_energy"=>cells.energy_input))
    arena,phase,medium,cells,flux
end
function sweep_run(out,name,condition,params;duration=160)
    a,p,m,c,f=boundary_fixture(name,condition,params)
    dir=joinpath(out,name*"-"*condition); mkpath(dir)
    previous=copy(p.phi); samples=Any[]; started=time()
    function sample!()
        push!(samples,merge(arena_measure(a,p,m,c,f),boundary_metrics(p.phi,previous,.5),
            (membrane_mass=sum(a.world.fields.deposit),forcing_rms=sqrt(sum(abs2,m.forcing)/length(m.forcing)))))
        copyto!(previous,p.phi)
    end
    sample!(); record_state(a.world,a.processes,joinpath(dir,"snapshot-0");step=0)
    steps=round(Int,duration/a.dt)
    for k in 1:steps
        arena_step!(a)
        if k%200==0 || k==steps
            sample!()
            println(name,"/",condition," t=",round(a.world.time;digits=1)," population=",length(c.ids)," components=",last(samples).components," interface=",last(samples).interface_length);flush(stdout)
        end
        if k==steps÷2 || k==steps
            record_state(a.world,a.processes,joinpath(dir,"snapshot-"*string(k));step=k)
            open(joinpath(dir,"replay-"*string(k)*".json"),"w") do io; JSON3.pretty(io,arena_snapshot(a)); end
        end
    end
    result=(name=name,condition=condition,settings=params,model_version="catabolic-fuel-fairshare-v2",sources=Interior.arena_sources(),wall_seconds=time()-started,samples=samples)
    open(joinpath(dir,"metrics.json"),"w") do io; JSON3.pretty(io,result);end
    result
end
function sweep_main()
    out=isempty(ARGS) ? "outputs/compartment-sweep" : ARGS[1]; mkpath(out)
    configurations=[("default",(gain=.3,bias=.1,pinning=.5,rate=.03,secretion=.025)),
        ("moderate",(gain=1.,bias=.5,pinning=2.,rate=.15,secretion=.08)),
        ("strong",(gain=4.,bias=.5,pinning=2.,rate=.15,secretion=.08))]
    results=Any[]
    for (name,params) in configurations
        push!(results,sweep_run(out,name,"full",params))
    end
    for condition in ("no_feedback","no_inhabitants")
        push!(results,sweep_run(out,"strong",condition,last(configurations)[2]))
    end
    open(joinpath(out,"summary.json"),"w") do io; JSON3.pretty(io,results);end
end
if abspath(PROGRAM_FILE)==@__FILE__
    sweep_main()
end
