#!/usr/bin/env julia
"""Joint paid reserve/investment/allocation policy with optional recurrent features."""
module NativeJointConstructionCredit
using Interior, KernelAbstractions
include("native_construction_credit.jl")
const C=NativeConstructionCredit
const L=C.L
@inline function feature(s,q,i,h,c,memory)
    c<=4 && return C.feature(s,i,h,c)
    c==5 && return Float32(s.energy[i]/max(s.mass[i]*q.energy_capacity,eps(Float32)))
    c<=7 && return Float32(s.metabolites[i,c-5]/max(s.mass[i]*q.store_capacity,eps(Float32)))
    !memory && return 0f0
    c<=7+h ? s.hidden[i,c-7] : s.adaptation[i,c-7-h]
end
@kernel function choose_joint!(s,status,q,z,h,decision,sigma,seed,clock,window_start,critic_mode,memory)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1
        changed=z.ids[i]!=s.ids[i]
        if changed
            z.ids[i]=s.ids[i];z.ready[i]=0;z.baseline[i]=0.;z.previous[i]=L.local_value(s,q,i)
            z.paid[i]=0.;z.last_reward[i]=0.;z.older_reward[i]=0.;z.prediction[i]=0.
            for c in 1:8;z.critic[i,c]=0f0;z.features[i,c]=0f0;end
            for a in 1:3,c in 1:size(z.learned,2);z.learned[i,c,a]=0f0;z.trace[i,c,a]=0f0;end
        end
        if critic_mode && (window_start || changed)
            f1=C.feature(s,i,h,2);f2=C.feature(s,i,h,3)
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
            for a in 1:3
                # Native reserve action tanh(atanh(-.6)) selects .2 of capacity.
                mu=a==1 ? -.8f0 : a==2 ? 0f0 : -0.6931472f0
                for c in 1:size(z.learned,2);mu+=z.learned[i,c,a]*feature(s,q,i,h,c,memory);end
                noise=sigma*Interior.dc_normal(seed,s.ids[i],clock,200+a)
                z.drive[i,a]=mu+noise
                if sigma>0
                    for c in 1:size(z.learned,2);z.trace[i,c,a]+=feature(s,q,i,h,c,memory)*noise/(sigma*sigma);end
                end
            end
        end
    end
end
@kernel function reserve!(s,status,z)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1
        s.actions[i,11]=tanh(z.drive[i,3])
    end
end
struct JointHabitat{H,M}
    base::H
    memory::M
end
# The study observer and common run protocol see the original state/clock API.
function Base.getproperty(h::JointHabitat,key::Symbol)
    key===:base || key===:memory ? getfield(h,key) : getproperty(getfield(h,:base),key)
end
Base.setproperty!(h::JointHabitat,key::Symbol,value)=setproperty!(getfield(h,:base),key,value)
function sidecar(p;memory=true,kwargs...)
    base=C.sidecar(p;kwargs...);n=Int(p.cells.geometry.capacity);f=7+2Int(p.cells.geometry.hidden)
    alloc(dims...)=fill!(similar(p.cells.state.inputs,Float32,dims...),0f0)
    z=merge(base.state,(learned=alloc(n,f,3),trace=alloc(n,f,3),drive=alloc(n,3)))
    JointHabitat(C.ConstructionHabitat(base.habitat,z,base.backend,base.options,base.original_masks,base.steps,base.launcher),memory)
end
Interior.device_anabolic_argument(h::JointHabitat)=Interior.device_anabolic_argument(h.base)
Interior.device_anabolic_begin!(h::JointHabitat,w,p,dt)=Interior.device_anabolic_begin!(h.base,w,p,dt)
Interior.device_anabolic_inherit!(h::JointHabitat,w,p,step)=Interior.device_anabolic_inherit!(h.base,w,p,step)
Interior.device_anabolic_finish!(h::JointHabitat,w,p,dt)=Interior.device_anabolic_finish!(h.base,w,p,dt)
function Interior.device_anabolic_sense!(h::JointHabitat,w,p,dt)
    dt==h.options.dt||error("joint construction clock changed")
    Interior.device_anabolic_sense!(h.habitat,w,p,dt)
    s=p.state;q=h.options
    v=(active=s.active,ids=s.ids,energy=s.energy,metabolites=s.metabolites,mass=s.mass,inputs=s.inputs,hidden=s.hidden,adaptation=s.adaptation)
    chemistry=(efficiency=p.parameters.efficiency,potentials=p.parameters.potentials,energy_capacity=p.parameters.energy_capacity,store_capacity=p.parameters.store_capacity)
    L.policy_launch!(h,choose_joint!,v,p.environment.status,chemistry,h.state,q.hidden,h.steps%q.decision_stride==0,
        q.sigma,q.seed,UInt32(h.steps÷q.decision_stride+1),h.steps%q.reward_stride==0,q.critic,h.memory;ndrange=Int(p.geometry.capacity))
end
function Interior.device_anabolic_prepare!(h::JointHabitat,w,p,dt)
    # Existing prepare has already produced native requests and the two held
    # material commands. Only reserve priority is replaced before metabolism.
    Interior.device_anabolic_prepare!(h.base,w,p,dt)
    L.policy_launch!(h,reserve!,(active=p.state.active,actions=p.state.actions),p.environment.status,h.state;ndrange=Int(p.geometry.capacity))
end
function run(profile,out;memory=true,kwargs...)
    factory=(p;options...)->sidecar(p;memory,options...)
    metadata=(law="joint_reserve_material",memory,actor_features=("bias","raw_balance","raw_total","catalyst_total",
        "energy_capacity_fraction","store_a_capacity_fraction","store_b_capacity_fraction","hidden","adaptation"),
        recurrent_features="2h existing hidden/adaptation features; zeroed when memory=false")
    C.run(profile,out;reserve_fraction=.2,make_sidecar=factory,driver_file=@__FILE__,policy_metadata=metadata,kwargs...)
end

end

if abspath(PROGRAM_FILE)==(@__FILE__)
    NativeJointConstructionCredit.run(ARGS[1],ARGS[2])
end
