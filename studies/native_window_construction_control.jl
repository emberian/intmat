#!/usr/bin/env julia
"""Original two-time window credit, with the cohort experiment's act/drain schedule."""
module NativeWindowConstructionControl
using Interior
include("native_joint_construction_credit.jl")
const J=NativeJointConstructionCredit
const C=J.C
struct WindowHabitat{H}
    base::H
    act_stride::Int
end
Base.getproperty(h::WindowHabitat,k::Symbol)=k===:base || k===:act_stride ? getfield(h,k) : getproperty(getfield(h,:base),k)
Base.setproperty!(h::WindowHabitat,k::Symbol,v)=setproperty!(getfield(h,:base),k,v)
sidecar(p;act_time=16.,kwargs...)=begin
    b=J.sidecar(p;kwargs...);WindowHabitat(b,round(Int,act_time/b.options.dt))
end
Interior.device_anabolic_argument(h::WindowHabitat)=Interior.device_anabolic_argument(h.base)
function Interior.device_anabolic_begin!(h::WindowHabitat,w,p,dt)
    if h.steps>=h.act_stride
        h.options=merge(h.options,(sigma=0f0,learning=false,critic=false,critic_learning=false))
    end
    Interior.device_anabolic_begin!(h.base,w,p,dt)
end
Interior.device_anabolic_sense!(h::WindowHabitat,w,p,dt)=Interior.device_anabolic_sense!(h.base,w,p,dt)
Interior.device_anabolic_prepare!(h::WindowHabitat,w,p,dt)=Interior.device_anabolic_prepare!(h.base,w,p,dt)
Interior.device_anabolic_inherit!(h::WindowHabitat,w,p,step)=Interior.device_anabolic_inherit!(h.base,w,p,step)
Interior.device_anabolic_finish!(h::WindowHabitat,w,p,dt)=Interior.device_anabolic_finish!(h.base,w,p,dt)
function run(profile,out;kwargs...)
    C.run(profile,out;reserve_fraction=.2,make_sidecar=sidecar,driver_file=@__FILE__,
        policy_metadata=(law="window_credit_matched_act_drain",act_time=16.,drain_end=24.,memory=true),kwargs...)
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    NativeWindowConstructionControl.run(ARGS[1],ARGS[2])
end
