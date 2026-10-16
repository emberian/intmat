#!/usr/bin/env julia
"""Continuously observed, finite physical supply ramps for paid construction.

At a raw-composition crossing, recent change predicts the forthcoming substrate
imbalance. Episodes remain externally initialized; no sensor is hidden or reset
inside an episode. The experiment is a control component, not a living ecology.
"""
module NativeRampConstruction
include("native_history_construction.jl")
using .NativeHistoryConstruction: Interior, KernelAbstractions, Random, Statistics, SHA
using KernelAbstractions, Random, Statistics, SHA
using Interior: JSON3, World
const H=NativeHistoryConstruction
const C=H.C
const L=H.L
const A=H.A
const T=H.T

# Both banks have squared norm exactly 2 in real arithmetic. At zero trend use
# +1, explicitly: absence of a trend is not an extra zero-amplitude feature.
@inline function actor_features(balance,trend,reactive,swap)
    x=tanh(balance/.2f0)
    direction=trend<0f0 ? -1f0 : 1f0
    direction=swap ? -direction : direction
    y=sqrt(max(0f0,1f0-x*x))
    (1f0,x,(reactive ? 1f0 : direction)*y)
end
@kernel function observe_choose!(s,status,z,m,q,step,noise,reactive,forced,swap)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1
        if step%q.stride==0
            T.observe!(m,i,s.ids[i],C.feature(s,i,q.hidden,2),0f0,0f0,q.af,q.as)
        end
        if step==q.decision
            balance=C.feature(s,i,q.hidden,2);trend=T.contrast(m,i,1)
            x=actor_features(balance,trend,reactive,swap)
            cx=actor_features(balance,trend,false,false)
            mu=0f0
            for k in 1:3
                mu+=z.learned[i,k,1]*x[k]
                z.trace[i,k,1]=x[k]*noise[i]/(q.sigma*q.sigma)
                z.actor_features[i,k]=x[k];z.features[i,k]=q.critic_scale*cx[k]
            end
            for k in 4:8;z.features[i,k]=0f0;end
            z.drive[i,1]=isfinite(forced) ? forced : mu+noise[i]
            prediction=0f0
            for k in 1:8;prediction+=z.critic[i,k]*z.features[i,k];end
            z.prediction[i]=Float64(prediction)
            z.previous[i]=L.local_value(s,q,i);z.ids[i]=s.ids[i]
            z.balance[i]=balance;z.trend[i]=trend
            z.raw_total[i]=C.feature(s,i,q.hidden,3)
        end
    end
end
struct RampHabitat{H}
    inner::H
end
Interior.device_anabolic_argument(h::RampHabitat)=Interior.device_anabolic_argument(h.inner)
Interior.device_anabolic_begin!(h::RampHabitat,w,p,dt)=Interior.device_anabolic_begin!(h.inner,w,p,dt)
Interior.device_anabolic_prepare!(h::RampHabitat,w,p,dt)=Interior.device_anabolic_prepare!(h.inner,w,p,dt)
Interior.device_anabolic_finish!(h::RampHabitat,w,p,dt)=Interior.device_anabolic_finish!(h.inner,w,p,dt)
Interior.device_anabolic_inherit!(h::RampHabitat,w,p,step)=Interior.device_anabolic_inherit!(h.inner,w,p,step)
function Interior.device_anabolic_sense!(h::RampHabitat,w,p,dt)
    a=h.inner
    Interior.device_anabolic_sense!(a.habitat,w,p,dt)
    observe_choose!(a.backend,64)(p.state,p.environment.status,a.state,a.memory,a.options,a.steps,a.noise,
        a.mode==1,a.forced,a.swap;ndrange=Int(p.geometry.capacity))
end
function fixture(profile,direction;seed=211,rate=.72,decision=1.,supply_duration=8.)
    direction in (1,2) || error("direction must name supplied species")
    cfg=A.fixture(profile,(.03,.03);seed).config
    cfg["cellular"]["speed"]=0.;cfg["catalysis"]["allocation_exponent"]=1.
    # Native feed runs before sensing: decision step 50 observes 51 grants.
    low=.8-rate*(decision+.02);low>=0 || error("negative initial stock")
    cfg["feedstock"]["niche_contrast"]=0.
    cfg["feedstock"]["initial"]=direction==1 ? [low,.8] : [.8,low]
    cfg["feedstock"]["source_rates"]=direction==1 ? [rate,0.] : [0.,rate]
    cfg["feedstock"]["source_budgets"]=direction==1 ? [supply_duration*rate,0.] : [0.,supply_duration*rate]
    w,p,_=A.native_constructed_garden(CPU(),cfg)
    (;grid=w.grid,seed=w.seed,checkpoint=Interior.device_habitat_checkpoint(w,p),config=cfg)
end
function episode(f,learned,critic;reactive=false,noise_seed=1,sigma=.5,forced=NaN32,swap=false,
        horizon=8.,decision=1.,pulse_time=.4,rate=0.,ready=true,train=false,explore=true,critic_scale=4.)
    dt=.02;w=World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed)
    p=Interior.device_habitat_restore(w,f.checkpoint)
    base=L.sidecar(p;dt,rate=0.,sigma,learning=false,critic=true)
    n=Int(p.cells.geometry.capacity);m=T.state(p.cells.state.inputs,n);tc=T.coefficients(.2)
    z=merge(base.state,(learned=copy(learned),critic=copy(critic),trace=zeros(Float32,n,3,1),drive=zeros(Float32,n,1),
        actor_features=zeros(Float32,n,3),balance=zeros(Float32,n),trend=zeros(Float32,n),raw_total=zeros(Float32,n)))
    z.ready.=Int32(ready)
    options=merge(p.cells.parameters,(hidden=Int(p.cells.geometry.hidden),stride=10,decision=round(Int,decision/dt),
        pulse_end=round(Int,(decision+pulse_time)/dt),af=tc.fast,as=tc.slow,sigma=Float32(sigma),critic_scale=Float32(critic_scale)))
    noise=Float32(sigma).*randn(MersenneTwister(noise_seed),Float32,n);!explore && fill!(noise,0f0)
    h=RampHabitat(H.HistoryHabitat(p.habitat,z,m,options,0,Int(reactive),Float32(forced),swap,CPU(),noise))
    origin=A.audit(p)
    for step in 1:round(Int,horizon/dt)
        Interior.device_cellular_step!(w,p.cells,dt;anabolism=h,check_status=false)
    end
    reward=[L.local_value(p.cells.state,p.cells.parameters,i)-z.previous[i] for i in 1:4]
    before=copy(z.learned)
    if train
        q=merge(p.cells.parameters,(critic_gain=.3f0,))
        L.reward!(CPU(),64)(p.cells.state,p.cells.environment.status,q,z,Float32(rate),2f0,.001f0,true,true,true;ndrange=n)
    end
    a=A.audit(p)
    material=a.material-origin.material-a.supplied+origin.supplied+a.removed-origin.removed
    potential=a.potential-origin.potential+a.work-origin.work+a.heat-origin.heat-a.supplied_potential+origin.supplied_potential+a.removed_potential-origin.removed_potential
    a.finite && abs(material)<1e-8 && abs(potential)<1e-8 || error("ramp physical closure")
    all(p.cells.state.active[1:4].==1) || error("fixed resident died")
    (;learned=z.learned,critic=z.critic,reward,prediction=z.prediction[1:4],paid=z.paid[1:4],balance=z.balance[1:4],trend=z.trend[1:4],
        raw_total=z.raw_total[1:4],features=z.actor_features[1:4,:],drive=z.drive[1:4,1],
        update_l1=sum(abs,z.learned.-before),funded=a.secreted,converted=a.solubilized,supplied=a.supplied-origin.supplied,material,potential)
end
compact(x)=H.compact(x)
function authority(profile,out)
    ispath(out)&&error("fresh output required");mkpath(out)
    write(joinpath(out,"executed-driver.jl"),read(@__FILE__))
    learned=zeros(Float32,64,3,1);critic=zeros(Float32,64,8)
    results=[(;direction,logit=u,supply_duration,compact(episode(fixture(profile,direction;supply_duration),learned,critic;forced=u,explore=false))...)
        for supply_duration in (8.,1.02) for direction in (1,2) for u in (-1.,0.,1.)]
    report=(;schema="native-ramp-authority-v1",results)
    open(io->JSON3.write(io,report),joinpath(out,"report.json"),"w");report
end
function run(profile,out;episodes=64,rate=.5,seed=211,noise_seed=9301,critic_scale=4.)
    iseven(episodes)&&episodes>=4 || error("balanced even episode count required")
    ispath(out)&&error("fresh output required");mkpath(out)
    write(joinpath(out,"executed-driver.jl"),read(@__FILE__))
    contexts=[(;direction,source_rate) for source_rate in (.72,.5) for direction in (1,2)]
    fixtures=[fixture(profile,c.direction;seed,rate=c.source_rate) for c in contexts]
    initial=zeros(Float32,64,3,1);critic0=zeros(Float32,64,8)
    order=vcat([shuffle(MersenneTwister(noise_seed+trial),[1,2]) for trial in 1:episodes÷2]...)
    arms=Dict{String,Any}();started=time()
    for (name,reactive,credit) in (("history",false,1.),("instantaneous",true,1.),("no_credit",false,0.),("reverse_credit",false,-1.))
        learned=copy(initial);critic=copy(critic0);history=Any[]
        for trial in 1:episodes
            r=episode(fixtures[order[trial]],learned,critic;reactive,noise_seed=noise_seed+trial,rate=rate*credit,train=true,ready=trial>2,critic_scale)
            learned=r.learned;critic=r.critic
            push!(history,(;trial,context=order[trial],compact(r)...))
        end
        probes=[(;context,swap,compact(episode(f,learned,critic;reactive,explore=false,swap,critic_scale))...)
            for (context,f) in enumerate(fixtures) for swap in (false,true)]
        arms[name]=(;history,probes,learned=learned[1:4,:,1],critic=critic[1:4,:],
            total_learning_work=sum(sum(r.paid) for r in history),total_update_l1=sum(r.update_l1 for r in history))
        open(io->JSON3.write(io,arms[name]),joinpath(out,name*".json"),"w")
        println((;arm=name,elapsed=time()-started));flush(stdout)
    end
    report=(;schema="native-ramp-construction-v1",seed,noise_seed,episodes,rate,critic_scale,contexts,arms,
        scope="external fresh physical episodes, continuous native observations, no within-episode masking; not autonomous living learning",
        features="history [1,x,sign_plus(trend)*sqrt(1-x²)] versus instantaneous [1,x,sqrt(1-x²)], x=tanh(balance/.2), norm²=2",
        schedule="real finite raw source continues past matched raw crossing at decision1; pulse.4; payoff8; training source.72, held-out.5",
        config=fixtures[1].config,source=fixtures[1].checkpoint.source,cellular_source=fixtures[1].checkpoint.cells.source,
        driver_sha256=bytes2hex(sha256(read(joinpath(out,"executed-driver.jl")))),elapsed_seconds=time()-started)
    open(io->JSON3.write(io,report),joinpath(out,"report.json"),"w");report
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    NativeRampConstruction.authority(ARGS[1],ARGS[2])
end
