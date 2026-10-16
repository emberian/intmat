#!/usr/bin/env julia
"""Stoichiometric developmental value through shared paid decision cohorts."""
module NativeDevelopmentLiving
include("native_temporal_living.jl")
include("local_development_value.jl")
using .NativeTemporalLiving: Interior,KernelAbstractions,TOML,SHA,JSON3
using KernelAbstractions
const V=NativeTemporalLiving
const Value=LocalDevelopmentValue
@kernel function capital_start!(s,status,z,slot,coefficient)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1 && z.queue_ids[i,slot]==s.ids[i]
        z.queue_value[i,slot]+=coefficient*Float64(s.mass[i])
    end
end
@kernel function capital_birth!(s,status,z,events,capital,coefficient,adjust)
    child=@index(Global,Linear)
    if status[1]==0 && s.active[child]==1
        parent=s.new_parent[child];id=s.ids[child]
        if parent>0 && id<=length(capital.recorded) && events.child_id[id]==id && capital.recorded[id]==0
            gift=Float64(s.mass[child]);capital.mass[id]=gift;capital.recorded[id]=1
            if adjust && s.active[parent]==1 && z.ids[parent]==s.ids[parent]
                for k in axes(z.queue_ids,2)
                    z.queue_ids[parent,k]==s.ids[parent] && (z.queue_value[parent,k]-=coefficient*gift)
                end
            end
        end
    end
end
struct DevelopmentHabitat{H,E}
    base::H
    capital::E
    coefficient::Float64
end
Base.getproperty(h::DevelopmentHabitat,k::Symbol)=k in (:base,:capital,:coefficient) ? getfield(h,k) : getproperty(getfield(h,:base),k)
Base.setproperty!(h::DevelopmentHabitat,k::Symbol,v)=setproperty!(getfield(h,:base),k,v)
function sidecar(p;coefficient=nothing,kwargs...)
    base=V.sidecar(p;kwargs...)
    lambda=coefficient===nothing ? Float64(Value.replacement_coefficient(p.cells.parameters)) : Float64(coefficient)
    isfinite(lambda)&&lambda>=0||error("invalid developmental value coefficient")
    options=merge(base.options,(development_coefficient=lambda,))
    c=V.C.ConstructionHabitat(base.habitat,base.state,base.backend,options,base.original_masks,base.steps,base.launcher)
    living=V.N.LivingHabitat(V.D.CohortHabitat(V.J.JointHabitat(c,false)),base.events,base.inheritance_fraction,base.birth_transfer_adjusted)
    base=V.TemporalHabitat(living,base.temporal,base.reactive,base.swap,base.forced)
    capital=(mass=fill!(similar(base.events.gift),0.),recorded=fill!(similar(base.events.child_id),0))
    DevelopmentHabitat(base,capital,lambda)
end
Interior.device_anabolic_argument(h::DevelopmentHabitat)=Interior.device_anabolic_argument(h.base)
Interior.device_anabolic_begin!(h::DevelopmentHabitat,w,p,dt)=Interior.device_anabolic_begin!(h.base,w,p,dt)
Interior.device_anabolic_prepare!(h::DevelopmentHabitat,w,p,dt)=Interior.device_anabolic_prepare!(h.base,w,p,dt)
function Interior.device_anabolic_sense!(h::DevelopmentHabitat,w,p,dt)
    Interior.device_anabolic_sense!(h.base,w,p,dt)
    q=h.options
    if h.steps%q.decision_stride==0 && h.steps<q.act_stride && (q.learning||q.critic)
        slot=mod(h.steps÷q.decision_stride,q.queue_length)+1
        V.L.policy_launch!(h,capital_start!,p.state,p.environment.status,h.state,slot,h.coefficient;ndrange=Int(p.geometry.capacity))
    end
end
function Interior.device_anabolic_inherit!(h::DevelopmentHabitat,w,p,step)
    Interior.device_anabolic_inherit!(h.base,w,p,step)
    V.L.policy_launch!(h,capital_birth!,p.state,p.environment.status,h.state,h.events,h.capital,h.coefficient,h.birth_transfer_adjusted;ndrange=Int(p.geometry.capacity))
end
Interior.device_anabolic_finish!(h::DevelopmentHabitat,w,p,dt)=Interior.device_anabolic_finish!(h.base,w,p,dt)
function observation(w,p,h,origin)
    o=V.observation(w,p,h,origin)
    merge(o,(cells=[merge(c,(development_value=c.usable+h.coefficient*c.mass,)) for c in o.cells],))
end
function returns(initial,final,h,start)
    u=V.returns(initial,final,V.N.birth_events(h),start);byid=Dict(c.id=>c for c in final.cells)
    e=map(Array,h.events);capital=map(Array,h.capital)
    [let old=initial.cells[k],r=u[k]
        gifts=sum(capital.mass[j] for j in eachindex(e.child_id) if capital.recorded[j]==1 && e.parent_id[j]==r.id && e.time[j]>start;init=0.)
        dm=(r.alive ? byid[r.id].mass : 0.)-old.mass+gifts
        merge(r,(development_value=r.all_initial_value+h.coefficient*dm,partition_mass=gifts))
    end for k in eachindex(u)]
end
function frozen_probe(w,p,h,origin,duration,mode)
    h.backend isa KernelAbstractions.CPU||error("development frozen probes require CPU")
    b=deepcopy((w=w,p=p,h=h));a=b.h
    a.options=merge(a.options,(sigma=0f0,learning=false,critic=false,critic_learning=false));fill!(a.state.queue_ids,0)
    mode==:erased&&fill!(a.state.learned,0f0)
    base=V.TemporalHabitat(a.base.base,a.temporal,a.reactive,mode==:swapped,a.forced)
    probe=DevelopmentHabitat(base,a.capital,a.coefficient)
    initial=observation(b.w,b.p,probe,origin);start=b.w.time
    V.L.advance!(b.w,b.p,probe,duration;dt=.02)
    final=observation(b.w,b.p,probe,origin)
    (;mode,initial,final,returns=returns(initial,final,probe,start))
end
function run(profile,out;duration=24.,probe_time=4.,rate=.1,noise_seed=6201)
    ispath(out)&&error("fresh output required");mkpath(out)
    write(joinpath(out,"executed-driver.jl"),read(@__FILE__))
    arms=Dict{String,Any}()
    for(name,reactive,credit) in (("history",false,1.),("instantaneous",true,1.),("no_credit",false,0.),("reversed_credit",false,-1.))
        f=V.N.fixture(profile);w,p=f.w,f.p
        h=sidecar(p;reactive,rate=rate*credit,learning=credit!=0,seed=noise_seed,horizon=4.,act_time=duration,inheritance_fraction=0.,
            (Symbol(k)=>v for(k,v) in f.cfg["temporal_policy"])...)
        origin=V.C.A.audit(p);samples=Any[observation(w,p,h,origin)]
        for _ in 1:round(Int,duration/2)
            V.L.advance!(w,p,h,2.;dt=.02);push!(samples,observation(w,p,h,origin))
        end
        probes=[frozen_probe(w,p,h,origin,probe_time,mode) for mode in (:learned,:erased,:swapped)]
        z=h.state
        arms[name]=(;samples,probes,events=V.N.birth_events(h),coefficient=h.coefficient,options=h.options,
            learning_work=sum(z.paid),actor_updates=sum(z.actor_updates),critic_updates=sum(z.critic_updates),clipping=sum(z.clip_total),
            critic_max=maximum(abs,z.critic),learned_max=maximum(abs,z.learned),capital_gifts=Array(h.capital.mass))
        open(io->JSON3.write(io,arms[name]),joinpath(out,name*".json"),"w")
        println((;arm=name,population=length(last(samples).cells),work=sum(z.paid)));flush(stdout)
    end
    report=(;schema="developmental-living-learning-v1",duration,probe_time,rate,noise_seed,arms,recipe=TOML.parsefile(profile),
        source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
        scope="new stoichiometric developmental accounting objective; physical energy ledger unchanged; same actor/horizon/physics as inventory-only comparison")
    open(io->JSON3.write(io,report),joinpath(out,"report.json"),"w");report
end
end
