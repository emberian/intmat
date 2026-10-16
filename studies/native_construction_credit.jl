#!/usr/bin/env julia
"""Paid local policy for habitat investment/allocation under accounted alternating supply.
The environmental clock is deliberately absent from actor and critic observations.
"""
module NativeConstructionCredit
using Interior, KernelAbstractions, SHA
using Interior: JSON3
include("native_lifetime_credit.jl")
const L=NativeLifetimeCredit
const A=L.NativeActionCredit
@inline feature(s,i,h,c)=c==1 ? 1f0 : s.inputs[i,2h+14+c]
@kernel function choose!(s,status,q,z,h,decision,sigma,seed,clock,window_start,critic_mode)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1
        changed=z.ids[i]!=s.ids[i]
        if changed
            z.ids[i]=s.ids[i];z.ready[i]=0;z.baseline[i]=0.;z.previous[i]=L.local_value(s,q,i)
            z.paid[i]=0.;z.last_reward[i]=0.;z.older_reward[i]=0.;z.prediction[i]=0.
            for c in 1:8;z.critic[i,c]=0f0;z.features[i,c]=0f0;end
            for a in 1:2,c in 1:4;z.learned[i,c,a]=0f0;z.trace[i,c,a]=0f0;end
        end
        if critic_mode && (window_start || changed)
            f1=feature(s,i,h,2);f2=feature(s,i,h,3)
            z.features[i,7]=f1-z.features[i,2];z.features[i,8]=f2-z.features[i,3]
            z.features[i,1]=1f0;z.features[i,2]=f1;z.features[i,3]=f2
            z.features[i,4]=Float32(L.local_value(s,q,i)/2)
            z.features[i,5]=Float32(z.last_reward[i]/.1)
            z.features[i,6]=Float32((z.last_reward[i]-z.older_reward[i])/.1)
            prediction=0f0
            for c in 1:8;prediction+=z.critic[i,c]*z.features[i,c];end
            z.prediction[i]=Float64(prediction)
        end
        if decision || changed
            for a in 1:2
                mu=a==1 ? -.8f0 : 0f0
                for c in 1:4;mu+=z.learned[i,c,a]*feature(s,i,h,c);end
                noise=sigma*Interior.dc_normal(seed,s.ids[i],clock,200+a)
                z.drive[i,a]=mu+noise
                if sigma>0
                    for c in 1:4;z.trace[i,c,a]+=feature(s,i,h,c)*noise/(sigma*sigma);end
                end
            end
        end
    end
end
@kernel function actuate!(a,s,status,q,z,dt,enabled)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1
        investment=enabled ? inv(1f0+exp(-clamp(z.drive[i,1],-40f0,40f0))) : 0f0
        allocation=inv(1f0+exp(-clamp(z.drive[i,2],-40f0,40f0)))
        Interior.dh_set_requests!(a,s,i,q,dt,investment,allocation)
        for k in 1:2;z.requested[i,k]+=Float64(a.requests[i,k]);end
    end
end
@kernel function feed_phase!(masks,original,swapped)
    i,k=@index(Global,NTuple)
    masks[i,k]=original[i,swapped && k<=2 ? 3-k : k]
end
mutable struct ConstructionHabitat{H,Z,B,Q,M,P}
    habitat::H
    state::Z
    backend::B
    options::Q
    original_masks::M
    steps::Int
    launcher::P
end
function sidecar(p;dt=.02,decision_time=.2,rate=.003,sigma=.3,seed=6201,learning=true,critic=true,enabled=true,phase_offset=0,launcher=L.NativePolicyLaunch())
    old=L.sidecar(p;dt,decision_time,rate,sigma,seed,learning,critic,critic_learning=learning || rate==0,launcher)
    C=Int(p.cells.geometry.capacity)
    alloc(dims...)=fill!(similar(p.cells.state.inputs,Float32,dims...),0f0)
    z=merge(old.state,(learned=alloc(C,4,2),trace=alloc(C,4,2),drive=alloc(C,2),requested=fill!(similar(p.cells.state.inputs,Float64,C,2),0.)))
    opts=merge(old.options,(phase_stride=round(Int,4/dt),phase_offset,enabled))
    masks=p.habitat.feedstock.arrays.source_masks
    totals=p.habitat.feedstock.parameters.source_totals
    totals[1]==totals[2] || error("mask permutation requires equal declared totals")
    ConstructionHabitat(p.habitat,z,p.cells.backend,opts,copy(masks),0,launcher)
end
Interior.device_anabolic_argument(h::ConstructionHabitat)=Interior.device_anabolic_argument(h.habitat)
function Interior.device_anabolic_begin!(h::ConstructionHabitat,w,p,dt)
    if h.steps%h.options.phase_stride==0
        swapped=isodd(h.steps÷h.options.phase_stride+h.options.phase_offset)
        masks=h.habitat.feedstock.arrays.source_masks
        L.policy_launch!(h,feed_phase!,masks,h.original_masks,swapped;ndrange=size(masks))
    end
    Interior.device_anabolic_begin!(h.habitat,w,p,dt)
end
function Interior.device_anabolic_sense!(h::ConstructionHabitat,w,p,dt)
    dt==h.options.dt||error("construction physics clock changed")
    Interior.device_anabolic_sense!(h.habitat,w,p,dt)
    s=p.state;q=h.options
    v=(active=s.active,ids=s.ids,energy=s.energy,metabolites=s.metabolites,inputs=s.inputs)
    chemistry=(efficiency=p.parameters.efficiency,potentials=p.parameters.potentials)
    L.policy_launch!(h,choose!,v,p.environment.status,chemistry,h.state,q.hidden,h.steps%q.decision_stride==0,
        q.sigma,q.seed,UInt32(h.steps÷q.decision_stride+1),h.steps%q.reward_stride==0,q.critic;ndrange=Int(p.geometry.capacity))
end
function Interior.device_anabolic_prepare!(h::ConstructionHabitat,w,p,dt)
    Interior.device_anabolic_prepare!(h.habitat,w,p,dt)
    a=h.habitat.arrays;q=h.habitat.parameters;s=p.state
    L.policy_launch!(h,actuate!,(expression=a.expression,requests=a.requests),(active=s.active,mass=s.mass),p.environment.status,
        (secretion_rate=q.secretion_rate,allocation_exponent=q.allocation_exponent),h.state,dt,h.options.enabled;ndrange=Int(p.geometry.capacity))
end
Interior.device_anabolic_inherit!(h::ConstructionHabitat,w,p,step)=Interior.device_anabolic_inherit!(h.habitat,w,p,step)
function Interior.device_anabolic_finish!(h::ConstructionHabitat,w,p,dt)
    Interior.device_anabolic_finish!(h.habitat,w,p,dt);h.steps+=1
    q=h.options
    if h.steps%q.reward_stride==0
        s=p.state;v=(active=s.active,ids=s.ids,energy=s.energy,metabolites=s.metabolites,ledger=s.ledger,ledger_low=s.ledger_low)
        chemistry=(efficiency=p.parameters.efficiency,potentials=p.parameters.potentials)
        L.policy_launch!(h,L.reward!,v,p.environment.status,chemistry,h.state,q.rate,.3f0,q.kappa,q.learning,q.critic,q.critic_learning;ndrange=Int(p.geometry.capacity))
    end
end
function fixture(profile;seed=211,feed=true,reserve_fraction=nothing)
    cfg=A.fixture(profile,(.03,.03);seed).config
    cfg["cellular"]["speed"]=0.;cfg["cellular"]["capacity"]=256
    cfg["feedstock"]["initial"]=[0.,0.];cfg["feedstock"]["niche_contrast"]=.95
    cfg["feedstock"]["source_rates"]=feed ? [.5,.5] : [0.,0.]
    cfg["feedstock"]["source_budgets"]=feed ? [14.,14.] : [0.,0.]
    cfg["catalysis"]["decay_rate"]=.2
    cfg["material_sensing"]=Dict("channels"=>["raw_balance","raw_total","catalyst_total"],"scales"=>[1.,1.,.5])
    w,p,_=A.native_constructed_garden(CPU(),cfg)
    if reserve_fraction!==nothing
        0<reserve_fraction<1 || error("reserve fraction must lie in (0,1)")
        row=Int(p.cells.geometry.hidden)+11
        p.cells.state.weights[:,:,row].=0f0;p.cells.state.marks[:,:,row].=0f0
        p.cells.state.weights[:,end,row].=Float32(atanh(2reserve_fraction-1))
    end
    (;grid=w.grid,seed=w.seed,checkpoint=Interior.device_habitat_checkpoint(w,p),config=cfg)
end
function observation(w,p,h,origin)
    snap=Interior.device_habitat_snapshot(p);a=Interior.device_habitat_audit(snap);s=snap.cells.state
    material=a.material-origin.material-a.supplied+origin.supplied+a.removed-origin.removed
    potential=a.potential-origin.potential+a.work-origin.work+a.heat-origin.heat-a.supplied_potential+origin.supplied_potential+a.removed_potential-origin.removed_potential
    a.finite && abs(material)<1e-8 && abs(potential)<1e-8 || error("construction accounting failure")
    q=p.cells.parameters;z=map(Array,h.state)
    ledger=Float64.(snap.arrays.cell_ledger).+Float64.(snap.arrays.cell_ledger_low)
    resident_valid=[s.active[i]==1 && s.ids[i]==i for i in 1:4]
    usable=[resident_valid[i] ? Float64(s.energy[i])+Float64(q.efficiency)*sum(Float64(s.metabolites[i,k])*Float64(q.potentials[k]) for k in 1:2) : nothing for i in 1:4]
    (;time=w.time,ids=s.ids[s.active.==1],usable,resident_valid,all_residents_valid=all(resident_valid),
        energy=s.energy[1:4],stores=s.metabolites[1:4,:],reserve_fraction=(s.actions[1:4,11].+1f0).*.5f0,
        learning_work=z.paid[1:4],work=a.work-origin.work,heat=a.heat-origin.heat,
        uptake_a=[Float64(s.ledger[i,3])+Float64(s.ledger_low[i,3]) for i in 1:4],
        uptake_b=[Float64(s.ledger[i,4])+Float64(s.ledger_low[i,4]) for i in 1:4],
        decisions=h.steps÷h.options.decision_stride,reward_windows=h.steps÷h.options.reward_stride,
        baseline=z.baseline[1:4],last_return=z.last_reward[1:4],prediction_error=z.prediction_error[1:4],ema_error=z.ema_error[1:4],
        credit_moment=z.credit_moment[1:4],clipped=z.clipped[1:4],critic_max=maximum(abs,z.critic),mark_max=maximum(abs,z.learned),material,potential,
        supplied=a.supplied-origin.supplied,supplied_potential=a.supplied_potential-origin.supplied_potential,
        removed=a.removed-origin.removed,removed_potential=a.removed_potential-origin.removed_potential,
        secretion_a=ledger[1:4,7],secretion_b=ledger[1:4,8],expression=snap.arrays.expression[1:4,:],
        requests=snap.arrays.requests[1:4,:],cumulative_requests=z.requested[1:4,:],enzyme=a.enzyme,bound=a.bound,solubilized=a.solubilized)
end

function run(profile,out;seed=211,noise_seed=6201,training=24.,test_time=4.,dt=.02,reserve_fraction=nothing,make_sidecar=sidecar,driver_file=@__FILE__,
        policy_metadata=(law="two_action",actor_features=("bias","raw_balance","raw_total","catalyst_total")))
    training>=4 && isapprox(training/4,round(training/4);atol=1e-8) || error("training requires whole four-time feed phases")
    test_time>0 && isapprox(test_time/4,round(test_time/4);atol=1e-8) || error("probe requires whole four-time feed phases")
    ispath(out)&&error("fresh output required");mkpath(out);write(joinpath(out,"executed-driver.jl"),read(driver_file))
    protocol_sha256=bytes2hex(sha256(read(@__FILE__)))
    driver_file!=(@__FILE__) && write(joinpath(out,"executed-protocol.jl"),read(@__FILE__))
    results=Dict{String,Any}()
    for (arm,rate,learning,enabled,feed) in (("credit",.003,true,true,true),("no_credit",0.,false,true,true),
            ("reversed_credit",-.003,true,true,true),("no_secretion",0.,false,false,true),("no_raw",.003,true,true,false))
        f=fixture(profile;seed,feed,reserve_fraction);w=World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed)
        p=Interior.device_habitat_restore(w,f.checkpoint);h=make_sidecar(p;dt,rate,learning,enabled,seed=noise_seed)
        origin=A.audit(p);trace=Any[]
        for _ in 1:round(Int,training/2)
            L.advance!(w,p,h,2.;dt);push!(trace,observation(w,p,h,origin))
            println((;arm,time=w.time,paid=sum(h.state.paid),mark=maximum(abs,h.state.learned)));flush(stdout)
        end
        cp=Interior.device_habitat_checkpoint(w,p);branches=Dict{String,Any}()
        for erase in (false,true),shift in (0,1)
            wr=World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed);pr=Interior.device_habitat_restore(wr,cp)
            hr=make_sidecar(pr;dt,rate=0.,sigma=0.,learning=false,critic=false,enabled,phase_offset=shift)
            # Continue the known external forcing clock, but expose none of it to the policy.
            copyto!(hr.original_masks,h.original_masks);hr.steps=h.steps
            hr.state.ids.=pr.cells.state.ids
            !erase && copyto!(hr.state.learned,h.state.learned)
            for i in eachindex(hr.state.previous);hr.state.previous[i]=L.local_value(pr.cells.state,pr.cells.parameters,i);end
            base=A.audit(pr);u0=A.usable(pr)
            valid0=[pr.cells.state.active[i]==1 && pr.cells.state.ids[i]==i for i in 1:4]
            L.advance!(wr,pr,hr,test_time;dt)
            observed=observation(wr,pr,hr,base);u1=A.usable(pr)
            returns=[valid0[i] && observed.resident_valid[i] ? u1[i]-u0[i] : nothing for i in 1:4]
            branches[(erase ? "erased" : "learned")*(shift==0 ? "_continued" : "_shifted")]=
                (;return_value=returns,observed...)
        end
        results[arm]=(;trace,branches,learned=Array(h.state.learned)[1:4,:,:])
        open(io->JSON3.write(io,results[arm]),joinpath(out,arm*".json"),"w")
    end
    report=(;schema="native-construction-credit-v1",seed,noise_seed,training,test_time,dt,decision_time=.2,reward_time=2.,feed_phase_time=4.,reserve_fraction,
        policy_metadata...,protocol_sha256,material_scales=(1.,1.,.5),
        actor_rate=.003,sigma=.3,kappa=.001,decay_rate=.2,source_rates=(.5,.5),source_budgets=(14.,14.),
        scope="study habitat policy; paid actor and predictor; no phase observations; matched frozen/erased physical states; no cooperation claim",
        driver_sha256=bytes2hex(sha256(read(joinpath(out,"executed-driver.jl")))),source=Interior.DEVICE_HABITAT_SOURCE_PIN,
        cellular_source=Interior.DEVICE_CELLULAR_SOURCE_PIN,results)
    open(io->JSON3.pretty(io,report),joinpath(out,"report.json"),"w");report
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    NativeConstructionCredit.run(ARGS[1],ARGS[2];training=length(ARGS)>2 ? parse(Float64,ARGS[3]) : 24.)
end
