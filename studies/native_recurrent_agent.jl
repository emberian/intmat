#!/usr/bin/env julia
"""Native observation/action/reward adapter for externally trained recurrent agents.

Shared trained parameters belong to the external trainer. Per-ID recurrent state
belongs to its runtime. This adapter never calls central training biological work.
"""
module NativeRecurrentAgent
include("native_predictive_programs.jl")
using .NativePredictivePrograms: Interior,KernelAbstractions
using KernelAbstractions
const P=NativePredictivePrograms
const N=P.N
const L=N.L
const INPUTS=25
const TARGETS=12
const ACTIVE_CHANNELS=(12,13,11)
const LOWER=(0f0,0f0,.02f0)
const UPPER=(.3f0,1f0,.25f0)
const INITIAL_BIAS=(0f0,0f0,Float32(log((.035-.02)/(.25-.035))))
@inline decode(u,lower,upper)=lower+(upper-lower)/(1f0+exp(-u))
@inline value(s,q,i)=L.local_value(s,q,i)+(Float64(q.efficiency)*(q.potentials[1]+q.potentials[2])/2+q.growth_work)*s.mass[i]
@kernel function pack!(s,status,q,z,a,step)
    i=@index(Global,Linear)
    if status[1]==0
        live=s.active[i]==1
        if a.observed_step[i]!=step
            a.incoming[i]=live && a.ids[i]==s.ids[i] ? 1f0 : 0f0
            a.ids[i]=live ? s.ids[i] : 0;a.observed_step[i]=step
        end
        a.valid[i]=live ? 1f0 : 0f0
        for d in 1:25;a.x[d,i]=0f0;end
        if live
            for k in 1:6;a.x[k,i]=a.obs[i,k];end
            a.x[7,i]=Float32(s.energy[i]/max(s.mass[i]*q.energy_capacity,eps(Float32)))
            a.x[8,i]=Float32(s.metabolites[i,1]/max(s.mass[i]*q.store_capacity,eps(Float32)))
            a.x[9,i]=Float32(s.metabolites[i,2]/max(s.mass[i]*q.store_capacity,eps(Float32)))
            a.x[10,i]=Float32(s.mass[i]/(q.division_mass+s.mass[i]))
            a.x[11,i]=s.age[i]/(20f0+s.age[i])
            a.x[12,i]=Float32(a.funded[i]/max(Float64(s.mass[i]),1e-6))
            for k in 1:11;a.x[12+k,i]=s.actions[i,k];end
            a.x[24,i]=a.command_ids[i]==s.ids[i] ? a.commands[i,12] : 0f0
            a.x[25,i]=a.command_ids[i]==s.ids[i] ? a.commands[i,13] : 0f0
        end
    end
end
@kernel function install!(s,status,q,a,raw,channels,lower,upper)
    i=@index(Global,Linear)
    if status[1]==0
        live=s.active[i]==1
        a.start_ids[i]=live ? s.ids[i] : 0;a.start_value[i]=live ? value(s,q,i) : 0.
        a.gifts[i]=0.;a.funded[i]=0.;a.command_ids[i]=live ? s.ids[i] : 0
        if live
            for j in eachindex(channels)
                u=raw[j,i]
                if isfinite(u)
                    a.commands[i,channels[j]]=decode(u,lower[j],upper[j])
                else
                    KernelAbstractions.@atomic status[1]|=Int32(8)
                end
            end
        end
    end
end
@kernel function override!(s,status,g,a,channels)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1 && a.command_ids[i]==s.ids[i]
        for j in eachindex(channels)
            channel=channels[j]
            if channel<=11
                # Initial reserve is encoded as a physical fraction; other native channels use [-1,1].
                s.outputs[i,Int(g.hidden)+channel]=channel==11 ? 2a.commands[i,channel]-1 : a.commands[i,channel]
            end
        end
    end
end
@kernel function requests!(s,status,habitat,a,hq,z)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1 && a.command_ids[i]==s.ids[i]
        Interior.dh_set_requests!(habitat,s,i,hq,z,a.commands[i,12],a.commands[i,13])
    end
end
@kernel function track_funded!(s,status,habitat,a)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1
        if a.funded_ids[i]!=s.ids[i];a.funded[i]=0.;a.funded_ids[i]=s.ids[i];end
        a.funded[i]+=Float64(habitat.last_secretion[i,1])+Float64(habitat.last_secretion[i,2])
    end
end
@kernel function birth_gift!(s,status,q,a)
    child=@index(Global,Linear)
    if status[1]==0 && s.active[child]==1
        parent=s.new_parent[child]
        id=s.ids[child]
        if parent>0 && id<=length(a.gift_recorded) && a.gift_recorded[id]==0
            a.gift_recorded[id]=id
            if s.active[parent]==1 && a.start_ids[parent]==s.ids[parent]
                # Native division emits at most one child per parent in a step.
                a.gifts[parent]+=value(s,q,child)
            end
        end
    end
end
@kernel function outcome!(s,status,q,a)
    i=@index(Global,Linear)
    if status[1]==0
        started=a.start_ids[i]!=0;same=started && s.active[i]==1 && s.ids[i]==a.start_ids[i]
        a.reward[i]=started ? Float32((same ? value(s,q,i) : 0.)-a.start_value[i]+a.gifts[i]) : 0f0
        a.outgoing[i]=same ? 1f0 : 0f0;a.predvalid[i]=same ? 1f0 : 0f0
    end
end
@kernel function targets!(a)
    i=@index(Global,Linear)
    for k in 1:12;a.target[k,i]=a.x[k,i];end
end
struct AgentHabitat{H,A,C}
    base::H
    agent::A
    channels::C
end
Base.getproperty(h::AgentHabitat,k::Symbol)=k in (:base,:agent,:channels) ? getfield(h,k) : getproperty(getfield(h,:base),k)
Base.setproperty!(h::AgentHabitat,k::Symbol,v)=setproperty!(getfield(h,:base),k,v)
Interior.device_anabolic_argument(h::AgentHabitat)=Interior.device_anabolic_argument(h.base)
Interior.device_anabolic_begin!(h::AgentHabitat,w,p,dt)=Interior.device_anabolic_begin!(h.base,w,p,dt)
function Interior.device_anabolic_sense!(h::AgentHabitat,w,p,dt)
    Interior.device_anabolic_sense!(h.habitat,w,p,dt)
    L.policy_launch!(h,P.sense!,p.state,h.habitat.feedstock.arrays.fields,p.environment.fields,h.habitat.arrays.enzyme_bits,Interior.device_support_descriptor(p),p.environment.status,p.geometry,h.agent,h.steps,1;ndrange=Int(p.geometry.capacity))
end
function Interior.device_controller_override!(h::AgentHabitat,w,p,dt)
    L.policy_launch!(h,override!,p.state,p.environment.status,p.geometry,h.agent,h.channels;ndrange=Int(p.geometry.capacity))
end
function Interior.device_anabolic_prepare!(h::AgentHabitat,w,p,dt)
    Interior.device_anabolic_prepare!(h.habitat,w,p,dt)
    L.policy_launch!(h,requests!,p.state,p.environment.status,h.habitat.arrays,h.agent,h.habitat.parameters,dt;ndrange=Int(p.geometry.capacity))
end
function Interior.device_anabolic_inherit!(h::AgentHabitat,w,p,step)
    Interior.device_anabolic_inherit!(h.base,w,p,step)
    L.policy_launch!(h,birth_gift!,p.state,p.environment.status,p.parameters,h.agent;ndrange=Int(p.geometry.capacity))
end
function Interior.device_anabolic_finish!(h::AgentHabitat,w,p,dt)
    Interior.device_anabolic_finish!(h.base,w,p,dt)
    L.policy_launch!(h,track_funded!,p.state,p.environment.status,h.habitat.arrays,h.agent;ndrange=Int(p.geometry.capacity))
end
struct AgentEnv{W,P,H,Q,O}
    w::W
    p::P
    h::H
    options::Q
    origin::O
end
Base.getproperty(e::AgentEnv,k::Symbol)=k in (:w,:p,:h,:options,:origin) ? getfield(e,k) : getproperty(getfield(e,:h).agent,k)
function make_env(profile;backend=KernelAbstractions.CPU(),decision_time=.2,phase_time=8.,phase_offset=0,seed=211,channels=ACTIVE_CHANNELS,lower=LOWER,upper=UPPER)
    length(channels)==length(lower)==length(upper)||error("decoder layout mismatch")
    allunique(channels)&&all(k->k in 1:13,channels)||error("invalid active channels")
    all(isfinite,lower)&&all(isfinite,upper)&&all(lower .< upper)||error("invalid bounds")
    all(j->lower[j]>=(channels[j]>=11 ? 0 : -1) && upper[j]<=1,eachindex(channels))||error("action bounds outside native domain")
    isfinite(phase_time)&&phase_time>0&&isfinite(decision_time)&&decision_time>0||error("invalid environment clocks")
    f=N.fixture(profile;backend,seed);w,p=f.w,f.p;dt=.02
    isapprox(round(Int,phase_time/dt)*dt,phase_time)||error("phase clock must divide native dt")
    b=N.sidecar(p;dt,decision_time,phase_offset,learning=false,critic=false,sigma=0.,rate=0.,inheritance_fraction=0.,horizon=2.,act_time=32.)
    b.options=merge(b.options,(phase_stride=round(Int,phase_time/dt),))
    n=Int(p.cells.geometry.capacity);alloc(T,d...)=fill!(similar(p.cells.state.inputs,T,d...),zero(T))
    a=(ids=alloc(Int32,n),observed_step=fill!(alloc(Int32,n),-1),gift_recorded=alloc(Int32,length(b.events.child_id)),obs=alloc(Float32,n,6),x=alloc(Float32,25,n),incoming=alloc(Float32,n),valid=alloc(Float32,n),
        commands=alloc(Float32,n,13),command_ids=alloc(Int32,n),start_ids=alloc(Int32,n),start_value=alloc(Float64,n),gifts=alloc(Float64,n),funded=alloc(Float64,n),funded_ids=alloc(Int32,n),
        reward=alloc(Float32,n),outgoing=alloc(Float32,n),target=alloc(Float32,12,n),predvalid=alloc(Float32,n))
    h=AgentHabitat(b,a,channels);origin=N.C.A.audit(p)
    env=AgentEnv(w,p,h,(;dt,decision_time,channels,lower,upper,phase_time,phase_offset,seed),origin)
    # One declared native initialization step supplies actual local receptors.
    L.advance!(w,p,h,dt;dt);observe!(env);env
end
function observe!(env)
    p=env.p.cells;h=env.h
    L.policy_launch!(h,pack!,p.state,p.environment.status,p.parameters,h.state,h.agent,h.steps;ndrange=Int(p.geometry.capacity))
    env
end
function act!(env,raw)
    size(raw)==(length(env.options.channels),Int(env.p.cells.geometry.capacity))||error("raw action shape")
    eltype(raw)===Float32||error("raw action dtype must be Float32")
    typeof(KernelAbstractions.get_backend(raw))===typeof(KernelAbstractions.get_backend(env.x))||error("raw action backend mismatch")
    p=env.p.cells;h=env.h
    L.policy_launch!(h,install!,p.state,p.environment.status,p.parameters,h.agent,raw,h.channels,env.options.lower,env.options.upper;ndrange=Int(p.geometry.capacity));env
end
function advance!(env)
    L.advance!(env.w,env.p,env.h,env.options.decision_time;dt=env.options.dt)
    p=env.p.cells;h=env.h
    L.policy_launch!(h,outcome!,p.state,p.environment.status,p.parameters,h.agent;ndrange=Int(p.geometry.capacity))
    observe!(env)
    L.policy_launch!(h,targets!,h.agent;ndrange=Int(p.geometry.capacity));env
end
function summary(env)
    f=N.frame(env.w,env.p,env.h,env.origin)
    (;observation_timing="species receptors: last physics-step start; body: completed step; previous11 native regulated actions plus two material request commands",time=env.w.time,population=length(f.cells),births=f.births,deaths=f.deaths,
      ledger=f.ledger,funded=f.funded_secretion,active_channels=env.options.channels,
      scope="shared external training task; per-ID state; developmental accounting reward, not fitness")
end

end
