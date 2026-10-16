#!/usr/bin/env julia
"""Decision-aligned delayed returns, using the original paid update transaction."""
module NativeCohortConstructionCredit
using Interior, KernelAbstractions
include("native_joint_construction_credit.jl")
const J=NativeJointConstructionCredit
const C=J.C
const L=J.L

@kernel function choose_cohort!(s,status,q,z,h,slot,sigma,seed,clock,acting,learning,critic,memory)
    i=@index(Global,Linear)
    record=acting && (learning || critic)
    if status[1]==0 && s.active[i]==1
        if z.ids[i]!=s.ids[i]
            z.ids[i]=s.ids[i];z.ready[i]=0;z.baseline[i]=0.;z.last_reward[i]=0.;z.older_reward[i]=0.
            for c in 1:8;z.critic[i,c]=0f0;z.features[i,c]=0f0;z.current_features[i,c]=0f0;end
            z.current_features[i,2]=C.feature(s,i,h,2);z.current_features[i,3]=C.feature(s,i,h,3)
            for a in 1:3,c in 1:size(z.learned,2);z.learned[i,c,a]=0f0;z.trace[i,c,a]=0f0;end
            for k in 1:size(z.queue_ids,2)
                z.invalidated[i]+=z.queue_ids[i,k]!=0;z.queue_ids[i,k]=0
            end
        end
        # Features are sampled before this decision's innovation. Old predictor
        # features live separately from whichever cohort just matured.
        f1=C.feature(s,i,h,2);f2=C.feature(s,i,h,3)
        z.current_features[i,7]=f1-z.current_features[i,2]
        z.current_features[i,8]=f2-z.current_features[i,3]
        z.current_features[i,1]=1f0;z.current_features[i,2]=f1;z.current_features[i,3]=f2
        z.current_features[i,4]=Float32(L.local_value(s,q,i)/2)
        z.current_features[i,5]=Float32(z.last_reward[i]/.1)
        z.current_features[i,6]=Float32((z.last_reward[i]-z.older_reward[i])/.1)
        prediction=0f0
        for c in 1:8;prediction+=z.critic[i,c]*z.current_features[i,c];end
        if record
            z.queue_ids[i,slot]=s.ids[i];z.queue_ready[i,slot]=z.ready[i]
            z.queue_value[i,slot]=L.local_value(s,q,i);z.queue_prediction[i,slot]=Float64(prediction)
            for c in 1:8;z.queue_critic[i,c,slot]=z.current_features[i,c];end
            for c in 1:size(z.learned,2);z.queue_features[i,c,slot]=J.feature(s,q,i,h,c,memory);end
            z.enqueued[i]+=1
        end
        for a in 1:3
            mu=a==1 ? -.8f0 : a==2 ? 0f0 : -0.6931472f0
            for c in 1:size(z.learned,2);mu+=z.learned[i,c,a]*J.feature(s,q,i,h,c,memory);end
            noise=sigma*Interior.dc_normal(seed,s.ids[i],clock,200+a)
            z.drive[i,a]=mu+noise
            if record;z.queue_score[i,a,slot]=sigma>0 ? noise/(sigma*sigma) : 0f0;end
        end
    end
end
@kernel function mature!(s,status,q,z,slot,rate)
    i=@index(Global,Linear)
    if status[1]==0
        z.transaction[i]=0;z.before_paid[i]=z.paid[i];z.proposal[i]=0.;z.critic_proposal[i]=0.
        if z.queue_ids[i,slot]!=0
            if s.active[i]==1 && z.queue_ids[i,slot]==s.ids[i] && z.ids[i]==s.ids[i]
                z.transaction[i]=1;z.matured[i]+=1
                z.ready[i]=z.queue_ready[i,slot];z.eligible[i]=z.ready[i];z.skipped[i]+=z.ready[i]==0
                z.previous[i]=z.queue_value[i,slot];z.prediction[i]=z.queue_prediction[i,slot]
                for c in 1:8;z.features[i,c]=z.queue_critic[i,c,slot];end
                advantage=L.local_value(s,q,i)-z.previous[i]-z.prediction[i]
                for a in 1:3,c in 1:size(z.learned,2)
                    g=z.queue_features[i,c,slot]*z.queue_score[i,a,slot];z.trace[i,c,a]=g
                    if z.ready[i]==1
                        old=z.learned[i,c,a]
                        z.raw_increment[i]+=abs(Float64(rate)*advantage*g)
                        z.proposal[i]+=abs(Float64(L.actor_proposal(old,g,advantage,rate,.3f0))-old)
                    end
                end
                norm=1f-4
                for c in 1:8;norm+=z.features[i,c]^2;end
                current=0f0
                for c in 1:8;current+=z.critic[i,c]*z.features[i,c];end
                criticerror=L.local_value(s,q,i)-z.previous[i]-Float64(current)
                step=q.critic_gain*Float32(criticerror)/norm
                for c in 1:8
                    old=z.critic[i,c]
                    z.critic_proposal[i]+=abs(Float64(L.critic_proposal(old,z.features[i,c],step))-old)
                end
            else
                z.invalidated[i]+=1
            end
            z.queue_ids[i,slot]=0
        end
    end
end
@kernel function transaction_diagnostics!(status,z,learning,critic,critic_learning)
    i=@index(Global,Linear)
    if status[1]==0 && z.transaction[i]==1
        paid=z.paid[i]>z.before_paid[i]
        z.actor_updates[i]+=paid && learning && z.proposal[i]>0
        z.critic_updates[i]+=paid && critic && critic_learning && z.critic_proposal[i]>0
        z.funded_increment[i]+=paid && learning ? z.proposal[i] : 0.
        z.proposed_increment[i]+=learning ? z.proposal[i] : 0.
        z.clip_total[i]+=learning && z.eligible[i]==1 ? z.clipped[i] : 0
    end
end
struct CohortHabitat{H}
    base::H
end
Base.getproperty(h::CohortHabitat,k::Symbol)=k===:base ? getfield(h,k) : getproperty(getfield(h,:base),k)
Base.setproperty!(h::CohortHabitat,k::Symbol,v)=setproperty!(getfield(h,:base),k,v)
function sidecar(p;horizon=2.,act_time=16.,memory=true,kwargs...)
    b=J.sidecar(p;memory,kwargs...);q=b.options
    K=round(Int,horizon/(q.dt*q.decision_stride));K>0 && isapprox(K*q.dt*q.decision_stride,horizon;atol=1e-6)||error("horizon must divide decisions")
    n=Int(p.cells.geometry.capacity);F=size(b.state.learned,2)
    alloc(T,d...)=fill!(similar(p.cells.state.inputs,T,d...),zero(T))
    z=merge(b.state,(queue_ids=alloc(Int32,n,K),queue_ready=alloc(Int32,n,K),queue_value=alloc(Float64,n,K),queue_prediction=alloc(Float64,n,K),
        queue_features=alloc(Float32,n,F,K),queue_score=alloc(Float32,n,3,K),queue_critic=alloc(Float32,n,8,K),current_features=alloc(Float32,n,8),
        transaction=alloc(Int32,n),eligible=alloc(Int32,n),before_paid=alloc(Float64,n),proposal=alloc(Float64,n),critic_proposal=alloc(Float64,n),raw_increment=alloc(Float64,n),enqueued=alloc(Int32,n),matured=alloc(Int32,n),
        skipped=alloc(Int32,n),invalidated=alloc(Int32,n),actor_updates=alloc(Int32,n),critic_updates=alloc(Int32,n),
        funded_increment=alloc(Float64,n),proposed_increment=alloc(Float64,n),clip_total=alloc(Int32,n)))
    options=merge(q,(horizon_stride=K*q.decision_stride,queue_length=K,act_stride=round(Int,act_time/q.dt),critic_gain=.03f0))
    base=C.ConstructionHabitat(b.habitat,z,b.backend,options,b.original_masks,b.steps,b.launcher)
    CohortHabitat(J.JointHabitat(base,memory))
end
Interior.device_anabolic_argument(h::CohortHabitat)=Interior.device_anabolic_argument(h.base)
Interior.device_anabolic_begin!(h::CohortHabitat,w,p,dt)=Interior.device_anabolic_begin!(h.base,w,p,dt)
Interior.device_anabolic_prepare!(h::CohortHabitat,w,p,dt)=Interior.device_anabolic_prepare!(h.base,w,p,dt)
Interior.device_anabolic_inherit!(h::CohortHabitat,w,p,step)=Interior.device_anabolic_inherit!(h.base,w,p,step)
function Interior.device_anabolic_sense!(h::CohortHabitat,w,p,dt)
    dt==h.options.dt||error("cohort clock changed")
    Interior.device_anabolic_sense!(h.habitat,w,p,dt)
    choose_step!(h,p)
end
function choose_step!(h,p)
    q=h.options
    if h.steps%q.decision_stride==0
        s=p.state
        v=(active=s.active,ids=s.ids,energy=s.energy,metabolites=s.metabolites,mass=s.mass,inputs=s.inputs,hidden=s.hidden,adaptation=s.adaptation)
        chemistry=(efficiency=p.parameters.efficiency,potentials=p.parameters.potentials,energy_capacity=p.parameters.energy_capacity,store_capacity=p.parameters.store_capacity)
        acting=h.steps<q.act_stride
        slot=mod(h.steps÷q.decision_stride,q.queue_length)+1
        L.policy_launch!(h,choose_cohort!,v,p.environment.status,chemistry,h.state,q.hidden,slot,acting ? q.sigma : 0f0,q.seed,
            UInt32(h.steps÷q.decision_stride+1),acting,q.learning,q.critic,h.memory;ndrange=Int(p.geometry.capacity))
    end
end
function Interior.device_anabolic_finish!(h::CohortHabitat,w,p,dt)
    Interior.device_anabolic_finish!(h.habitat,w,p,dt);h.steps+=1
    q=h.options
    if h.steps>=q.horizon_stride && h.steps%q.decision_stride==0
        s=p.state;chemistry=(efficiency=p.parameters.efficiency,potentials=p.parameters.potentials,critic_gain=q.critic_gain,critic_current=true)
        hasproperty(q,:development_coefficient) && (chemistry=merge(chemistry,(development_coefficient=q.development_coefficient,)))
        slot=mod((h.steps-q.horizon_stride)÷q.decision_stride,q.queue_length)+1
        L.policy_launch!(h,mature!,(active=s.active,ids=s.ids,energy=s.energy,metabolites=s.metabolites,mass=s.mass),p.environment.status,chemistry,h.state,slot,q.rate;ndrange=Int(p.geometry.capacity))
        # Transaction mask prevents the shared updater from inventing a return
        # for a missing, dead or replaced resident's queue entry.
        v=(active=h.state.transaction,ids=s.ids,energy=s.energy,metabolites=s.metabolites,mass=s.mass,ledger=s.ledger,ledger_low=s.ledger_low)
        L.policy_launch!(h,L.reward!,v,p.environment.status,chemistry,h.state,q.rate,.3f0,q.kappa,q.learning,q.critic,q.critic_learning;ndrange=Int(p.geometry.capacity))
        L.policy_launch!(h,transaction_diagnostics!,p.environment.status,h.state,q.learning,q.critic,q.critic_learning;ndrange=Int(p.geometry.capacity))
    end
end
function C.observation(w,p,h::CohortHabitat,origin)
    observed=invoke(C.observation,Tuple{Any,Any,Any,Any},w,p,h,origin)
    z=map(Array,h.state)
    merge(observed,(cohorts_enqueued=z.enqueued[1:4],cohorts_matured=z.matured[1:4],
        cohorts_pending=[count(!iszero,z.queue_ids[i,:]) for i in 1:4],cohorts_invalidated=z.invalidated[1:4],
        calibration_skips=z.skipped[1:4],actor_updates=z.actor_updates[1:4],critic_updates=z.critic_updates[1:4],
        raw_increment=z.raw_increment[1:4],proposed_increment=z.proposed_increment[1:4],funded_increment=z.funded_increment[1:4],clip_total=z.clip_total[1:4]))
end
function run(profile,out;horizon=2.,memory=true,kwargs...)
    factory=(p;options...)->sidecar(p;horizon,memory,options...)
    C.run(profile,out;reserve_fraction=.2,make_sidecar=factory,driver_file=@__FILE__,
        policy_metadata=(law="decision_aligned_cohorts",revision=2,critic_fit="current coefficients on stored pre-action features",horizon,act_time=16.,drain_end=24.,memory,
            calibration="pre-action queue ready flag; first horizon of decisions calibrates only; earliest actor update at twice horizon",critic_gain=.03),kwargs...)
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    NativeCohortConstructionCredit.run(ARGS[1],ARGS[2];horizon=length(ARGS)>2 ? parse(Float64,ARGS[3]) : 2.)
end
