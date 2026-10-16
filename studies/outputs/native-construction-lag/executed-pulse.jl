#!/usr/bin/env julia
"""Marginal three-actuator lag map around actual learned joint-policy states."""
module NativeConstructionLag
using Interior, KernelAbstractions, SHA
using Interior: JSON3
include("native_joint_construction_credit.jl")
const J=NativeJointConstructionCredit
const C=J.C
const L=J.L
@kernel function perturb_drive!(s,status,z,saved,axis,amplitude,producer,restore)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1 && s.ids[i]==producer
        if restore
            z.drive[i,axis]=saved[i,axis]
        else
            saved[i,axis]=z.drive[i,axis]
            z.drive[i,axis]+=amplitude
        end
    end
end
struct PulseHabitat{H,A,D,P,S}
    base::H
    axis::A
    amplitude::D
    producer::P
    saved_drive::S
    origin_step::Int
    duration_steps::Int
end
const PULSE_FIELDS=(:base,:axis,:amplitude,:producer,:saved_drive,:origin_step,:duration_steps)
Base.getproperty(h::PulseHabitat,key::Symbol)=key in PULSE_FIELDS ? getfield(h,key) : getproperty(getfield(h,:base),key)
Base.setproperty!(h::PulseHabitat,key::Symbol,value)=setproperty!(getfield(h,:base),key,value)
function pulse(h;axis=1,amplitude=0.,producer=1,duration=.2)
    1<=axis<=3 || error("pulse axis must be investment, allocation or reserve")
    stride=round(Int,duration/Float64(h.options.dt))
    isapprox(stride*Float64(h.options.dt),duration;atol=1e-7) && stride>0 || error("pulse clock must divide dt")
    PulseHabitat(h,Int32(axis),Float32(amplitude),Int32(producer),zero(h.state.drive),h.steps,stride)
end
Interior.device_anabolic_argument(h::PulseHabitat)=Interior.device_anabolic_argument(h.base)
Interior.device_anabolic_begin!(h::PulseHabitat,w,p,dt)=Interior.device_anabolic_begin!(h.base,w,p,dt)
Interior.device_anabolic_sense!(h::PulseHabitat,w,p,dt)=Interior.device_anabolic_sense!(h.base,w,p,dt)
Interior.device_anabolic_inherit!(h::PulseHabitat,w,p,step)=Interior.device_anabolic_inherit!(h.base,w,p,step)
Interior.device_anabolic_finish!(h::PulseHabitat,w,p,dt)=Interior.device_anabolic_finish!(h.base,w,p,dt)
function Interior.device_anabolic_prepare!(h::PulseHabitat,w,p,dt)
    active=h.steps-h.origin_step<h.duration_steps
    if active
        L.policy_launch!(h,perturb_drive!,(active=p.state.active,ids=p.state.ids),p.environment.status,h.state,h.saved_drive,
            h.axis,h.amplitude,h.producer,false;ndrange=Int(p.geometry.capacity))
    end
    Interior.device_anabolic_prepare!(h.base,w,p,dt)
    if active
        # Restore the exact saved sample, avoiding Float32 add/subtract drift.
        # Physical requests and reserve actions retain the perturbed commands.
        L.policy_launch!(h,perturb_drive!,(active=p.state.active,ids=p.state.ids),p.environment.status,h.state,h.saved_drive,
            h.axis,h.amplitude,h.producer,true;ndrange=Int(p.geometry.capacity))
    end
end
end
