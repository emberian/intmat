#!/usr/bin/env julia
"""Batched matched-state lag assay; no continuing run is reset or migrated."""
module ConstructionLagMap
using Interior, KernelAbstractions, SHA, Serialization
using Interior: JSON3
include("native_construction_lag.jl")
include("construction_lag_world_batches.jl")
const P=NativeConstructionLag
const J=P.J
const C=P.C
const L=P.L
const START_TIMES=(7.8,8.2,10.,11.8,12.2,14.)
const CAPTURE_TIMES=(0.,.2,.4,.6,1.,2.,4.,8.)
function checkpoints(profile;seeds=(211,217),starts=START_TIMES,dt=.02)
    records=[]
    for seed in seeds
        f=C.fixture(profile;seed,reserve_fraction=.2)
        w=World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed);p=Interior.device_habitat_restore(w,f.checkpoint)
        h=J.sidecar(p;dt,seed=6201);previous=0.
        for t in starts
            L.advance!(w,p,h,t-previous;dt);previous=t
            resident_valid=[p.cells.state.active[i]==1 && p.cells.state.ids[i]==i for i in 1:4]
            s=p.cells.state;q=p.cells.parameters;hidden=Int(p.cells.geometry.hidden)
            features=[J.feature(s,q,i,hidden,c,true) for i in 1:4,c in 1:size(h.state.learned,2)]
            feature_norm_squared=hcat([sum(abs2,features[i,1:7]) for i in 1:4],[sum(abs2,features[i,:]) for i in 1:4])
            logit_displacement=[sum(h.state.learned[i,c,a]*features[i,c] for c in axes(features,2)) for i in 1:4,a in 1:3]
            push!(records,(;seed,start=t,resident_valid,feature_norm_squared,logit_displacement,grid=f.grid,checkpoint=Interior.device_habitat_checkpoint(w,p),
                learned=copy(h.state.learned),original_masks=copy(h.original_masks),steps=h.steps))
        end
    end
    records
end
function member(record,axis,amplitude;dt=.02)
    w=World(record.grid,CPU(),NamedTuple(),NamedTuple(),record.start,record.seed)
    p=Interior.device_habitat_restore(w,record.checkpoint)
    h=J.sidecar(p;dt,rate=0.,sigma=.3,learning=false,critic=false,seed=6201)
    h.steps=record.steps;copyto!(h.original_masks,record.original_masks);copyto!(h.state.learned,record.learned)
    copyto!(h.state.ids,p.cells.state.ids)
    for i in eachindex(h.state.previous);h.state.previous[i]=L.local_value(p.cells.state,p.cells.parameters,i);end
    (w,p,P.pulse(h;axis,amplitude))
end
function observation(x,origin,offset)
    a=Interior.device_habitat_audit(x.physical);s=x.physical.cells.state;q=x.physical.cells.parameters
    material=a.material-origin.material-a.supplied+origin.supplied+a.removed-origin.removed
    potential=a.potential-origin.potential+a.work-origin.work+a.heat-origin.heat-a.supplied_potential+origin.supplied_potential+a.removed_potential-origin.removed_potential
    a.finite && abs(material)<1e-8 && abs(potential)<1e-8 || error("lag-map conservation failure")
    valid=[s.active[i]==1 && s.ids[i]==i for i in 1:4]
    usable=[valid[i] ? L.local_value(s,q,i) : nothing for i in 1:4]
    cl=Float64.(s.ledger).+Float64.(s.ledger_low)
    hl=Float64.(x.physical.arrays.cell_ledger).+Float64.(x.physical.arrays.cell_ledger_low)
    (;offset,usable,resident_valid=valid,energy=s.energy[1:4],stores=s.metabolites[1:4,:],
        uptake_a=cl[1:4,3],uptake_b=cl[1:4,4],oxidized=cl[1:4,7],
        secretion_a=hl[1:4,7],secretion_b=hl[1:4,8],requested=x.state.requested[1:4,:],
        reserve_fraction=(s.actions[1:4,11].+1f0).*.5f0,expression=x.physical.arrays.expression[1:4,:],
        enzyme=a.enzyme,bound=a.bound,solubilized=a.solubilized,work=a.work-origin.work,heat=a.heat-origin.heat,
        supplied=a.supplied-origin.supplied,supplied_potential=a.supplied_potential-origin.supplied_potential,
        removed=a.removed-origin.removed,removed_potential=a.removed_potential-origin.removed_potential,material,potential)
end
function cohort(records;backend=CPU(),dt=.02,captures=CAPTURE_TIMES)
    specs=[];worlds=World[];physical=[];learners=[]
    for r in records,(axis,amplitude) in ((1,0.),(1,-.15),(1,.15),(2,-.15),(2,.15),(3,-.15),(3,.15))
        w,p,h=member(r,axis,amplitude;dt)
        push!(worlds,w);push!(physical,p);push!(learners,h)
        push!(specs,(;seed=r.seed,start=r.start,initial_resident_valid=r.resident_valid,axis,amplitude,producer=1))
    end
    batch=learning_world_batch(physical,learners;backend,identities=["$(s.seed)-$(s.start)-$(s.axis)-$(s.amplitude)" for s in specs])
    worlds=[World(w.grid,backend,NamedTuple(),NamedTuple(),w.time,w.seed) for w in worlds]
    initial=learning_world_batch_snapshot(batch);origins=[Interior.device_habitat_audit(x.physical) for x in initial]
    traces=[Any[observation(x,o,0.)] for (x,o) in zip(initial,origins)]
    previous=0.
    for t in captures[2:end]
        steps=round(Int,(t-previous)/dt);isapprox(steps*dt,t-previous;atol=1e-8)||error("capture clock")
        learning_world_batch_advance!(worlds,batch,steps,dt)
        snapshots=learning_world_batch_snapshot(batch)
        for (trace,x,o) in zip(traces,snapshots,origins);push!(trace,observation(x,o,t));end
        previous=t
    end
    [(;specs[i]...,trace=traces[i]) for i in eachindex(specs)]
end
function run(profile,out;backend=CPU(),seeds=(211,217),starts=START_TIMES,captures=CAPTURE_TIMES,dt=.02)
    ispath(out)&&error("fresh output required");mkpath(out)
    maximum(starts)+maximum(captures)<=28 || error("assay exceeds original finite source duration")
    write(joinpath(out,"executed-driver.jl"),read(@__FILE__))
    write(joinpath(out,"executed-pulse.jl"),read(joinpath(@__DIR__,"native_construction_lag.jl")))
    haskey(WORLD_SIDECAR_POLICIES,:construction_lag) || register_construction_lag_world_batches!(P)
    records=checkpoints(profile;seeds,starts,dt)
    serialize(joinpath(out,"private-checkpoints.bin"),records)
    results=[]
    for t in starts
        members=filter(r->r.start==t,records)
        result=cohort(members;backend,dt,captures);append!(results,result)
        open(io->JSON3.write(io,result),joinpath(out,"start-$(t).json"),"w")
        println((;start=t,worlds=length(result),horizon=last(captures)));flush(stdout)
    end
    report=(;schema="native-construction-lag-v1",seeds,starts,captures,dt,amplitude=.15,pulse_time=.2,noise_seed=6201,sigma=.3,
        source=Interior.DEVICE_HABITAT_SOURCE_PIN,cellular_source=Interior.DEVICE_CELLULAR_SOURCE_PIN,
        driver_sha256=bytes2hex(sha256(read(joinpath(out,"executed-driver.jl")))),pulse_sha256=bytes2hex(sha256(read(joinpath(out,"executed-pulse.jl")))),
        contexts=[(;seed=r.seed,start=r.start,resident_valid=r.resident_valid,feature_norm_squared=r.feature_norm_squared,logit_displacement=r.logit_displacement) for r in records],
        scope="one held-logit perturbation of producer1 from learned physical checkpoints, then frozen policy and common future innovations; original remaining budgets and source clocks retained",results)
    open(io->JSON3.pretty(io,report),joinpath(out,"report.json"),"w");report
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    ConstructionLagMap.run(ARGS[1],ARGS[2])
end
