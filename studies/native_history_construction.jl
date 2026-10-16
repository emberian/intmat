#!/usr/bin/env julia
"""Engineered delayed-material-cue allocation task on paid native chemistry.

Independent physical episodes, explicit policy-sensor occlusion, one scored
allocation decision. This is training a useful material skill, not a claim of
autonomous lifetime or evolutionary learning. Native sensing remains untouched.
"""
module NativeHistoryConstruction
using Interior, KernelAbstractions, Random, Statistics, SHA
using Interior: JSON3
include("native_construction_credit.jl")
include("local_temporal_features.jl")
const C=NativeConstructionCredit
const L=C.L
const A=C.A
const T=LocalTemporalFeatures

@kernel function observe_decide!(s,status,z,m,q,step,noise,mode,forced,swap)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1
        if step%q.stride==0
            cue=step<q.cue_end ? C.feature(s,i,q.hidden,2) : 0f0
            T.observe!(m,i,s.ids[i],cue,0f0,0f0,q.af,q.as)
        end
        if step==q.decision
            trend=T.contrast(m,i,1)
            feature=trend/(.02f0+abs(trend))
            feature=swap ? -feature : feature
            signed_feature=feature
            feature=mode==1 ? abs(feature) : feature
            mu=z.learned[i,1,1]+z.learned[i,2,1]*feature
            u=isfinite(forced) ? forced : mu+noise[i]
            z.drive[i,1]=u;z.action_feature[i]=feature;z.signed_feature[i]=signed_feature
            z.trace[i,1,1]=noise[i]/(q.sigma*q.sigma)
            z.trace[i,2,1]=feature*noise[i]/(q.sigma*q.sigma)
            # The same critic receives pre-action local cue history in all arms;
            # its prediction is frozen before the allocation innovation.
            z.features[i,1]=1f0;z.features[i,2]=signed_feature
            for k in 3:8;z.features[i,k]=0f0;end
            prediction=0f0
            for k in 1:8;prediction+=z.critic[i,k]*z.features[i,k];end
            z.prediction[i]=Float64(prediction)
            z.previous[i]=L.local_value(s,q,i)
            z.ids[i]=s.ids[i]
        end
    end
end
@kernel function actuate!(a,s,status,z,q,dt,pulse)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1
        allocation=inv(1f0+exp(-clamp(z.drive[i,1],-40f0,40f0)))
        Interior.dh_set_requests!(a,s,i,q,dt,pulse ? .4f0 : 0f0,allocation)
        s.actions[i,11]=-.6f0 # same reserve fraction .2 in every arm
    end
end
mutable struct HistoryHabitat{H,Z,M,Q,B,V}
    habitat::H
    state::Z
    memory::M
    options::Q
    steps::Int
    mode::Int
    forced::Float32
    swap::Bool
    backend::B
    noise::V
end
Interior.device_anabolic_argument(h::HistoryHabitat)=Interior.device_anabolic_argument(h.habitat)
Interior.device_anabolic_begin!(h::HistoryHabitat,w,p,dt)=Interior.device_anabolic_begin!(h.habitat,w,p,dt)
Interior.device_anabolic_inherit!(h::HistoryHabitat,w,p,step)=Interior.device_anabolic_inherit!(h.habitat,w,p,step)
function Interior.device_anabolic_sense!(h::HistoryHabitat,w,p,dt)
    Interior.device_anabolic_sense!(h.habitat,w,p,dt)
    observe_decide!(h.backend,64)(p.state,p.environment.status,h.state,h.memory,h.options,h.steps,h.noise,h.mode,h.forced,h.swap;ndrange=Int(p.geometry.capacity))
end
function Interior.device_anabolic_prepare!(h::HistoryHabitat,w,p,dt)
    Interior.device_anabolic_prepare!(h.habitat,w,p,dt)
    q=h.habitat.parameters
    actuate!(h.backend,64)(h.habitat.arrays,p.state,p.environment.status,h.state,q,dt,
        h.options.decision<=h.steps<h.options.pulse_end;ndrange=Int(p.geometry.capacity))
end
function Interior.device_anabolic_finish!(h::HistoryHabitat,w,p,dt)
    Interior.device_anabolic_finish!(h.habitat,w,p,dt);h.steps+=1
end
function fixture(profile,context;seed=211)
    cfg=A.fixture(profile,(.03,.03);seed).config
    cfg["cellular"]["speed"]=0.;cfg["cellular"]["capacity"]=64
    cfg["catalysis"]["allocation_exponent"]=1. # constant total requested construction across allocation
    cfg["feedstock"]["initial"]=collect(context);cfg["feedstock"]["niche_contrast"]=0.
    w,p,_=A.native_constructed_garden(CPU(),cfg)
    (;grid=w.grid,seed=w.seed,checkpoint=Interior.device_habitat_checkpoint(w,p),config=cfg)
end
function episode(f,learned,critic;mode=0,noise_seed=1,sigma=.5,forced=NaN32,swap=false,horizon=8.,decision=1.,cue_time=.4,pulse_time=.4,rate=0.,ready=true,train=false,explore=true)
    dt=.02
    w=World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed);p=Interior.device_habitat_restore(w,f.checkpoint)
    base=L.sidecar(p;dt,rate=0.,sigma,learning=false,critic=true)
    n=Int(p.cells.geometry.capacity);m=T.state(p.cells.state.inputs,n);tc=T.coefficients(.2)
    z=merge(base.state,(learned=copy(learned),critic=copy(critic),trace=zeros(Float32,n,2,1),drive=zeros(Float32,n,1),
        action_feature=zeros(Float32,n),signed_feature=zeros(Float32,n)))
    z.ready.=Int32(ready)
    options=merge(p.cells.parameters,(hidden=Int(p.cells.geometry.hidden),stride=10,cue_end=round(Int,cue_time/dt),decision=round(Int,decision/dt),
        pulse_end=round(Int,(decision+pulse_time)/dt),af=tc.fast,as=tc.slow,sigma=Float32(sigma)))
    rng=MersenneTwister(noise_seed);noise=Float32(sigma).*randn(rng,Float32,n)
    !explore && fill!(noise,0f0)
    h=HistoryHabitat(p.habitat,z,m,options,0,mode,Float32(forced),swap,CPU(),noise)
    origin=A.audit(p)
    for step in 1:round(Int,horizon/dt)
        Interior.device_cellular_step!(w,p.cells,dt;anabolism=h,check_status=false)
    end
    reward=[L.local_value(p.cells.state,p.cells.parameters,i)-z.previous[i] for i in 1:4]
    advantage=reward.-z.prediction[1:4];before=copy(z.learned)
    if train
        q=merge(p.cells.parameters,(critic_gain=.3f0,))
        L.reward!(CPU(),64)(p.cells.state,p.cells.environment.status,q,z,Float32(rate),2f0,.001f0,true,true,true;ndrange=n)
    end
    a=A.audit(p)
    material=a.material-origin.material-a.supplied+origin.supplied+a.removed-origin.removed
    potential=a.potential-origin.potential+a.work-origin.work+a.heat-origin.heat-a.supplied_potential+origin.supplied_potential+a.removed_potential-origin.removed_potential
    a.finite && abs(material)<1e-8 && abs(potential)<1e-8 || error("history physical closure")
    all(p.cells.state.active[1:4].==1) || error("fixed resident died")
    (;learned=z.learned,critic=z.critic,reward,advantage,postpaid_usable=A.usable(p),paid=z.paid[1:4],
        feature=z.action_feature[1:4],signed_feature=z.signed_feature[1:4],drive=z.drive[1:4,1],
        update_l1=sum(abs,z.learned.-before),funded=a.secreted,converted=a.solubilized,material,potential)
end
compact(x)=NamedTuple{Tuple(filter(k->!(k in (:learned,:critic)),collect(keys(x))))}(x)
function authority(profile,out)
    ispath(out)&&error("fresh output required");mkpath(out)
    write(joinpath(out,"executed-driver.jl"),read(@__FILE__))
    fixtures=[fixture(profile,c) for c in ((.8,.08),(.08,.8))]
    learned=zeros(Float32,64,2,1);critic=zeros(Float32,64,8)
    results=[(;context=i,logit=u,compact(episode(f,learned,critic;forced=Float32(u)))...) for (i,f) in enumerate(fixtures) for u in (-1.,0.,1.)]
    report=(;schema="history-construction-authority-v1",contexts=((.8,.08),(.08,.8)),results)
    open(io->JSON3.write(io,report),joinpath(out,"report.json"),"w");report
end

function run(profile,out;episodes=128,rate=.5,seed=211,noise_seed=9101)
    iseven(episodes)&&episodes>=4 || error("balanced even episode count required")
    ispath(out)&&error("fresh output required");mkpath(out)
    write(joinpath(out,"executed-driver.jl"),read(@__FILE__))
    contexts=((.8,.08),(.08,.8),(.6,.12),(.12,.6))
    fixtures=[fixture(profile,c;seed) for c in contexts]
    initial=zeros(Float32,64,2,1);critic0=zeros(Float32,64,8)
    before=[compact(episode(f,initial,critic0;explore=false)) for f in fixtures]
    order=vcat([shuffle(MersenneTwister(noise_seed+trial),[1,2]) for trial in 1:episodes÷2]...)
    arms=Dict{String,Any}();started=time()
    for (name,mode,sign) in (("history",0,1.),("unsigned_history",1,1.),("no_credit",0,0.),("reverse_credit",0,-1.))
        learned=copy(initial);critic=copy(critic0);history=Any[]
        for trial in 1:episodes
            result=episode(fixtures[order[trial]],learned,critic;mode,noise_seed=noise_seed+trial,rate=rate*sign,train=true,ready=trial>2)
            learned=result.learned;critic=result.critic
            push!(history,(;trial,context=order[trial],compact(result)...))
            if trial%32==0
                println((;arm=name,trial,mean_reward=mean(result.reward),weights=vec(learned[1:4,2,1]),elapsed=time()-started));flush(stdout)
            end
        end
        probes=Any[]
        for (context,f) in enumerate(fixtures),swap in (false,true),decision in (1.,1.4)
            result=episode(f,learned,critic;mode,explore=false,swap,decision)
            push!(probes,(;context,swap,decision,compact(result)...))
        end
        arms[name]=(;history,probes,learned=learned[1:4,:,1],critic=critic[1:4,:],
            total_learning_work=sum(sum(r.paid) for r in history),total_update_l1=sum(r.update_l1 for r in history))
        open(io->JSON3.write(io,arms[name]),joinpath(out,name*".json"),"w")
    end
    report=(;schema="native-history-construction-v1",seed,noise_seed,episodes,rate,contexts,before,arms,
        scope="episodic external training with fresh physical states; one local Gaussian allocation score and pre-action baseline; no autonomous lifetime claim",
        architecture="two actor features: bias and signed or absolute fast-minus-slow local raw-balance history; policy sensor occluded after .4; decision1; pulse.4; return8",
        costs="native fixed-total paid construction, shared .2 reserve, actor/critic turnover .001; physical episode reset and reinstallation are external training apparatus",
        config=fixtures[1].config,source=fixtures[1].checkpoint.source,cellular_source=fixtures[1].checkpoint.cells.source,
        driver_sha256=bytes2hex(sha256(read(joinpath(out,"executed-driver.jl")))),elapsed_seconds=time()-started)
    open(io->JSON3.write(io,report),joinpath(out,"report.json"),"w");report
end

function transfer_probe(profile,trained,out;seed=211)
    ispath(out)&&error("fresh output required");mkpath(out)
    saved=JSON3.read(read(joinpath(trained,"history.json"),String))
    learned=zeros(Float32,64,2,1);learned[1:4,:,1].=reshape(Float32.(saved.learned),4,2)
    critic=zeros(Float32,64,8);results=Any[]
    for (context,c) in enumerate(((.8,.08),(.08,.8)))
        f=fixture(profile,c;seed)
        for (name,weights,cue,pulse) in (("learned_occluded",learned,.4,.4),("learned_unoccluded",learned,20.,.4),
                ("zero_readout",zero(learned),20.,.4),("no_construction",zero(learned),20.,0.))
            r=episode(f,weights,critic;explore=false,cue_time=cue,pulse_time=pulse)
            push!(results,(;context,name,compact(r)...))
        end
    end
    report=(;schema="native-history-transfer-v1",results,seed,
        purpose="unoccluded sensing in the same physical assay, not a continuous mobile deployment; no-construction opportunity control",
        trained_sha256=bytes2hex(sha256(read(joinpath(trained,"history.json")))))
    open(io->JSON3.write(io,report),joinpath(out,"report.json"),"w");report
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    NativeHistoryConstruction.authority(ARGS[1],ARGS[2])
end
