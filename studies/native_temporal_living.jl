#!/usr/bin/env julia
"""Continuing living allocation policy; native payment, birth and cohort mechanics."""
module NativeTemporalLiving
include("native_living_construction.jl")
include("local_temporal_features.jl")
using .NativeLivingConstruction: Interior, KernelAbstractions, TOML, SHA
using KernelAbstractions, Statistics
using Interior: JSON3
const N=NativeLivingConstruction
const D=N.D
const J=N.J
const C=N.C
const L=N.L
const T=LocalTemporalFeatures
@inline function feature_bank(balance,trend,reactive,swap=false)
    x=tanh(balance/.2f0);y=sqrt(max(0f0,1f0-x*x));sign=trend<0f0 ? -1f0 : 1f0
    (1f0,x,(reactive ? (balance<0f0 ? -1f0 : 1f0) : (swap ? -sign : sign))*y)
end
@kernel function choose!(s,status,q,z,m,opts,slot,clock,step,observe,decision,acting,reactive,swap,forced)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1 && (observe || decision || z.birth_pending[i]==1 || z.ids[i]!=s.ids[i])
        newcomer=z.birth_pending[i]==1 || z.ids[i]!=s.ids[i]
        if z.ids[i]!=s.ids[i]
            z.ids[i]=s.ids[i];z.ready[i]=0;z.last_reward[i]=0.;z.older_reward[i]=0.;z.paid[i]=0.
            for c in 1:8;z.critic[i,c]=0f0;z.features[i,c]=0f0;end
            for a in 1:3,c in 1:3;z.learned[i,c,a]=0f0;z.trace[i,c,a]=0f0;end
            for k in axes(z.queue_ids,2);z.invalidated[i]+=z.queue_ids[i,k]!=0;z.queue_ids[i,k]=0;end
        end
        balance=C.feature(s,i,opts.hidden,2);total=C.feature(s,i,opts.hidden,3);catalyst=C.feature(s,i,opts.hidden,4)
        interval=Float32(step-m.last_step[i])*opts.dt
        af=-expm1(-interval/.4f0);as=-expm1(-interval/2f0)
        T.observe!(m,i,s.ids[i],balance,total,catalyst,af,as);m.last_step[i]=step
        if decision || newcomer
        trend=T.contrast(m,i,1);x=feature_bank(balance,trend,reactive,swap);cx=feature_bank(balance,trend,false)
        # Shared predictor information in both arms. Explicit bounded physical basis.
        z.current_features[i,1]=opts.critic_scale
        z.current_features[i,2]=opts.critic_scale*cx[2];z.current_features[i,3]=opts.critic_scale*cx[3]
        z.current_features[i,4]=Float32(s.energy[i]/max(s.mass[i]*q.energy_capacity,eps(Float32)))
        z.current_features[i,5]=Float32((s.metabolites[i,1]-s.metabolites[i,2])/max(s.mass[i]*q.store_capacity,eps(Float32)))
        z.current_features[i,6]=total;z.current_features[i,7]=catalyst;z.current_features[i,8]=tanh(Float32(z.last_reward[i]))
        prediction=0f0
        for c in 1:8;prediction+=z.critic[i,c]*z.current_features[i,c];end
        record=decision && acting && (opts.learning || opts.critic)
        if record
            z.queue_ids[i,slot]=s.ids[i];z.queue_ready[i,slot]=z.ready[i]
            z.queue_value[i,slot]=L.local_value(s,q,i);z.queue_prediction[i,slot]=Float64(prediction)
            for c in 1:8;z.queue_critic[i,c,slot]=z.current_features[i,c];end
            for c in 1:3;z.queue_features[i,c,slot]=x[c];end
            z.enqueued[i]+=1
        end
        mu=0f0
        for c in 1:3;mu+=z.learned[i,c,2]*x[c];end
        sigma=decision && acting ? opts.sigma : 0f0
        noise=sigma*Interior.dc_normal(opts.seed,s.ids[i],clock,202)
        z.drive[i,1]=opts.investment_logit;z.drive[i,2]=isfinite(forced) ? forced : mu+noise;z.drive[i,3]=opts.reserve_logit
        if record
            z.queue_score[i,1,slot]=0f0;z.queue_score[i,3,slot]=0f0
            z.queue_score[i,2,slot]=sigma>0 ? noise/(sigma*sigma) : 0f0
        end
        z.birth_pending[i]=0
        end
    end
end
struct TemporalHabitat{H,M}
    base::H
    temporal::M
    reactive::Bool
    swap::Bool
    forced::Float32
end
Base.getproperty(h::TemporalHabitat,k::Symbol)=k in (:base,:temporal,:reactive,:swap,:forced) ? getfield(h,k) : getproperty(getfield(h,:base),k)
Base.setproperty!(h::TemporalHabitat,k::Symbol,v)=setproperty!(getfield(h,:base),k,v)
function sidecar(p;reactive=false,swap=false,forced=NaN32,critic_scale=4.,critic_gain=.3,investment=inv(1+exp(.8))/2,reserve=.2,phase_time=4.,observation_time=.2,kwargs...)
    p.habitat.parameters.allocation_exponent==1 || error("fixed-total allocation requires exponent1 recipe")
    b=N.sidecar(p;kwargs...);n=Int(p.cells.geometry.capacity);K=b.options.queue_length
    alloc(d...)=fill!(similar(p.cells.state.inputs,Float32,d...),0f0)
    z=merge(b.state,(learned=alloc(n,3,3),trace=alloc(n,3,3),queue_features=alloc(n,3,K)))
    observation_stride=round(Int,observation_time/b.options.dt)
    observation_stride>0 && isapprox(observation_stride*b.options.dt,observation_time;atol=1e-6) || error("observation clock must divide physics clock")
    b.options.decision_stride%observation_stride==0 || error("observation clock must divide decision clock")
    tc=T.coefficients(Float64(b.options.dt)*b.options.decision_stride)
    opts=merge(b.options,(observation_stride=observation_stride,critic_gain=Float32(critic_gain),phase_stride=round(Int,phase_time/b.options.dt),af=tc.fast,as=tc.slow,critic_scale=Float32(critic_scale),investment_logit=Float32(log(investment/(1-investment))),reserve_logit=Float32(atanh(2reserve-1))))
    c=C.ConstructionHabitat(b.habitat,z,b.backend,opts,b.original_masks,b.steps,b.launcher)
    living=N.LivingHabitat(D.CohortHabitat(J.JointHabitat(c,false)),b.events,b.inheritance_fraction,b.birth_transfer_adjusted)
    temporal=merge(T.state(p.cells.state.inputs,n),(last_step=fill!(similar(p.cells.state.inputs,Int32,n),0),))
    TemporalHabitat(living,temporal,reactive,swap,Float32(forced))
end
Interior.device_anabolic_argument(h::TemporalHabitat)=Interior.device_anabolic_argument(h.base)
Interior.device_anabolic_begin!(h::TemporalHabitat,w,p,dt)=Interior.device_anabolic_begin!(h.base,w,p,dt)
Interior.device_anabolic_prepare!(h::TemporalHabitat,w,p,dt)=Interior.device_anabolic_prepare!(h.base,w,p,dt)
Interior.device_anabolic_finish!(h::TemporalHabitat,w,p,dt)=Interior.device_anabolic_finish!(h.base,w,p,dt)
Interior.device_anabolic_inherit!(h::TemporalHabitat,w,p,step)=Interior.device_anabolic_inherit!(h.base,w,p,step)
function Interior.device_anabolic_sense!(h::TemporalHabitat,w,p,dt)
    dt==h.options.dt || error("temporal living clock changed")
    Interior.device_anabolic_sense!(h.habitat,w,p,dt)
    q=h.options;s=p.state
    v=(active=s.active,ids=s.ids,energy=s.energy,metabolites=s.metabolites,mass=s.mass,inputs=s.inputs)
    chemistry=(efficiency=p.parameters.efficiency,potentials=p.parameters.potentials,energy_capacity=p.parameters.energy_capacity,store_capacity=p.parameters.store_capacity)
    slot=mod(h.steps÷q.decision_stride,q.queue_length)+1
    L.policy_launch!(h,choose!,v,p.environment.status,chemistry,h.state,h.temporal,q,slot,UInt32(h.steps÷q.decision_stride+1),h.steps,
        h.steps%q.observation_stride==0,h.steps%q.decision_stride==0,h.steps<q.act_stride,h.reactive,h.swap,h.forced;ndrange=Int(p.geometry.capacity))
end
function observation(w,p,h,origin)
    f=N.frame(w,p,h,origin);s=map(Array,p.cells.state);z=map(Array,h.state);m=map(Array,h.temporal)
    cells=[(;id=s.ids[i],balance=C.feature(s,i,h.options.hidden,2),total=C.feature(s,i,h.options.hidden,3),trend=T.contrast(m,i,1),
        mass=s.mass[i],energy=s.energy[i],stores=vec(s.metabolites[i,:]),generation=s.generations[i],parent_id=s.parents[i],
        usable=L.local_value(s,p.cells.parameters,i),allocation=z.drive[i,2],predictor=vec(z.critic[i,:]),last_return=z.last_reward[i],
        prediction_error=z.prediction_error[i],actor_updates=z.actor_updates[i],critic_updates=z.critic_updates[i],clipped=z.clip_total[i],
        learned=vec(z.learned[i,:,2])) for i in eachindex(s.active) if s.active[i]==1]
    (;time=w.time,cells,ledger=f.ledger,funded=f.funded_secretion,requested=f.cumulative_requests,births=f.births,deaths=f.deaths)
end
function authority(profile,out;checkpoints=(2.,4.,6.),horizon=4.)
    ispath(out)&&error("fresh output required");mkpath(out)
    write(joinpath(out,"executed-driver.jl"),read(@__FILE__))
    f=N.fixture(profile);w,p=f.w,f.p;h=sidecar(p;rate=0.,learning=false,critic=false,sigma=0.,act_time=40.,inheritance_fraction=0.,horizon,
        (Symbol(k)=>v for (k,v) in get(f.cfg,"temporal_policy",Dict()))...)
    origin=C.A.audit(p);rows=Any[];last=0.
    for t in checkpoints
        L.advance!(w,p,h,t-last;dt=.02);last=t
        initial=observation(w,p,h,origin)
        for force in (-1f0,0f0,1f0)
            branch=deepcopy((w=w,p=p,h=h));a=branch.h
            pulse=TemporalHabitat(a.base,a.temporal,a.reactive,a.swap,force)
            samples=Any[];L.advance!(branch.w,branch.p,pulse,.4;dt=.02)
            # Continue ordinary deterministic decisions after the one paid pulse.
            for lag in (1.,2.,horizon)
                elapsed=branch.w.time-t;L.advance!(branch.w,branch.p,a,lag-elapsed;dt=.02)
                obs=observation(branch.w,branch.p,a,origin)
                gifts=N.birth_events(a)
                targets=[(;id=c.id,value=c.usable-old.usable+sum(e.child_usable_before_copy for e in gifts if e.parent_id==c.id && e.time>t;init=0.))
                    for c in obs.cells for old in initial.cells if c.id==old.id]
                push!(samples,(;lag,observation=obs,transfer_adjusted_returns=targets))
            end
            push!(rows,(;checkpoint=t,force,initial,samples))
        end
    end
    report=(;schema="temporal-living-authority-v1",rows,recipe=TOML.parsefile(profile),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
        driver_sha256=bytes2hex(SHA.sha256(read(joinpath(out,"executed-driver.jl")))))
    open(io->JSON3.write(io,report),joinpath(out,"report.json"),"w");report
end
function returns(initial,final,events,start)
    byid=Dict(c.id=>c for c in final.cells)
    [let alive=haskey(byid,c.id)
        gift=sum(e.child_usable_before_copy for e in events if e.parent_id==c.id && e.time>start;init=0.)
        value=(alive ? byid[c.id].usable : 0.)-c.usable+gift
        (;id=c.id,alive,value=alive ? value : nothing,all_initial_value=value)
    end for c in initial.cells]
end
function frozen_probe(w,p,h,origin,duration,mode)
    h.backend isa KernelAbstractions.CPU || error("frozen_probe requires CPU: GPU forks need explicit physical/private/temporal restoration")
    branch=deepcopy((w=w,p=p,h=h));a=branch.h
    a.options=merge(a.options,(sigma=0f0,learning=false,critic=false,critic_learning=false))
    fill!(a.state.queue_ids,0) # explicit diagnostic freeze: no pending training transactions
    mode==:erased && fill!(a.state.learned,0f0)
    probe=TemporalHabitat(a.base,a.temporal,a.reactive,mode==:swapped,a.forced)
    initial=observation(branch.w,branch.p,probe,origin);start=branch.w.time
    L.advance!(branch.w,branch.p,probe,duration;dt=.02)
    final=observation(branch.w,branch.p,probe,origin)
    (;mode,initial,final,returns=returns(initial,final,N.birth_events(probe),start))
end
function run(profile,out;duration=24.,probe_time=4.,horizon=4.,rate=.1,noise_seed=6201)
    ispath(out)&&error("fresh output required");mkpath(out)
    write(joinpath(out,"executed-driver.jl"),read(@__FILE__))
    arms=Dict{String,Any}()
    for (name,reactive,credit) in (("history",false,1.),("instantaneous",true,1.),("no_credit",false,0.),("reversed_credit",false,-1.))
        f=N.fixture(profile);w,p=f.w,f.p
        h=sidecar(p;reactive,rate=rate*credit,learning=credit!=0,seed=noise_seed,horizon,act_time=duration,
            inheritance_fraction=0.,(Symbol(k)=>v for (k,v) in get(f.cfg,"temporal_policy",Dict()))...)
        origin=C.A.audit(p);samples=Any[observation(w,p,h,origin)]
        for _ in 1:round(Int,duration/2)
            L.advance!(w,p,h,2.;dt=.02);push!(samples,observation(w,p,h,origin))
        end
        probes=[frozen_probe(w,p,h,origin,probe_time,mode) for mode in (:learned,:erased,:swapped)]
        z=h.state
        arms[name]=(;samples,probes,events=N.birth_events(h),options=h.options,
            learning_work=sum(z.paid),actor_updates=sum(z.actor_updates),critic_updates=sum(z.critic_updates),
            pending=count(!iszero,z.queue_ids),matured=sum(z.matured),invalidated=sum(z.invalidated),calibration_skips=sum(z.skipped),
            clipping=sum(z.clip_total),critic_max=maximum(abs,z.critic),learned_max=maximum(abs,z.learned))
        open(io->JSON3.write(io,arms[name]),joinpath(out,name*".json"),"w")
        println((;arm=name,population=length(last(samples).cells),learning_work=sum(z.paid)));flush(stdout)
    end
    report=(;schema="temporal-living-learning-v1",duration,probe_time,horizon,rate,noise_seed,arms,recipe=TOML.parsefile(profile),
        source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
        driver_sha256=bytes2hex(SHA.sha256(read(joinpath(out,"executed-driver.jl")))),
        scope="uninterrupted moving/dividing lives; fixed-total paid allocation; finite external sources/sinks; frozen probes share each arm's exact physical state")
    open(io->JSON3.write(io,report),joinpath(out,"report.json"),"w");report
end

end
