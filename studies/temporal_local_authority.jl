#!/usr/bin/env julia
"""Single-resident allocation intervention from a replayed completed checkpoint."""
module TemporalLocalAuthority
include("native_temporal_living.jl")
using .NativeTemporalLiving: Interior, KernelAbstractions, TOML, SHA, JSON3
using KernelAbstractions, Statistics, Serialization
const V=NativeTemporalLiving
@kernel function pulse!(s,status,z,saved,producer,amplitude,restore)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1 && (producer==0 || s.ids[i]==producer)
        if restore
            z.drive[i,2]=saved[i]
        else
            saved[i]=z.drive[i,2];z.drive[i,2]+=amplitude
        end
    end
end
struct Pulse{H,S}
    base::H
    saved::S
    producer::Int32
    amplitude::Float32
    stop_step::Int
end
Base.getproperty(h::Pulse,k::Symbol)=k in (:base,:saved,:producer,:amplitude,:stop_step) ? getfield(h,k) : getproperty(getfield(h,:base),k)
Base.setproperty!(h::Pulse,k::Symbol,v)=setproperty!(getfield(h,:base),k,v)
Interior.device_anabolic_argument(h::Pulse)=Interior.device_anabolic_argument(h.base)
Interior.device_anabolic_begin!(h::Pulse,w,p,dt)=Interior.device_anabolic_begin!(h.base,w,p,dt)
Interior.device_anabolic_sense!(h::Pulse,w,p,dt)=Interior.device_anabolic_sense!(h.base,w,p,dt)
Interior.device_anabolic_finish!(h::Pulse,w,p,dt)=Interior.device_anabolic_finish!(h.base,w,p,dt)
Interior.device_anabolic_inherit!(h::Pulse,w,p,step)=Interior.device_anabolic_inherit!(h.base,w,p,step)
function Interior.device_anabolic_prepare!(h::Pulse,w,p,dt)
    active=h.steps<h.stop_step
    s=(active=p.state.active,ids=p.state.ids)
    active && V.L.policy_launch!(h,pulse!,s,p.environment.status,h.state,h.saved,h.producer,h.amplitude,false;ndrange=Int(p.geometry.capacity))
    Interior.device_anabolic_prepare!(h.base,w,p,dt)
    active && V.L.policy_launch!(h,pulse!,s,p.environment.status,h.state,h.saved,h.producer,h.amplitude,true;ndrange=Int(p.geometry.capacity))
end
function observe(w,p,h,origin)
    o=V.observation(w,p,h,origin);s=p.cells.state;a=V.C.A.audit(p)
    cl=Float64.(s.ledger).+Float64.(s.ledger_low)
    cells=[(;id=c.id,usable=c.usable,balance=c.balance,total=c.total,trend=c.trend,mass=c.mass,
        energy=c.energy,stores=c.stores,position=vec(s.position[findfirst(==(c.id),s.ids),:]),
        uptake_a=cl[findfirst(==(c.id),s.ids),3],uptake_b=cl[findfirst(==(c.id),s.ids),4]) for c in o.cells]
    (;time=w.time,cells,ledger=o.ledger,funded=o.funded,requested=o.requested,births=o.births,deaths=o.deaths,
        catalyst=a.enzyme,bound=a.bound,converted=a.solubilized)
end
function run(profile,reference,out;checkpoint=12.,lags=(.2,.5,1.,2.,4.))
    ispath(out)&&error("fresh output required");mkpath(out)
    f=V.N.fixture(profile);w,p=f.w,f.p
    h=V.sidecar(p;rate=.1,seed=6201,horizon=4.,act_time=24.,inheritance_fraction=0.,
        (Symbol(k)=>v for(k,v) in f.cfg["temporal_policy"])...)
    origin=V.C.A.audit(p);V.L.advance!(w,p,h,checkpoint;dt=.02)
    initial=observe(w,p,h,origin)
    old=JSON3.read(read(reference,String));sample=only(filter(x->isapprox(x.time,checkpoint;atol=1e-8),old.samples))
    [c.id for c in initial.cells]==[Int(c.id) for c in sample.cells]||error("checkpoint identities differ")
    replay_error=maximum(maximum(abs.([c.mass-c0.mass,c.energy-c0.energy,c.usable-c0.usable,c.balance-Float32(c0.balance),c.trend-Float32(c0.trend)])) for(c,c0) in zip(initial.cells,sample.cells))
    replay_error<1e-12||error("checkpoint replay differs: $replay_error")
    order=sort(initial.cells;by=c->(c.balance,c.id));selected=order[[1,cld(length(order),2),length(order)]]
    h.options=merge(h.options,(learning=false,critic=false,critic_learning=false));fill!(h.state.queue_ids,0)
    serialize(joinpath(out,"private-checkpoint.bin"),(w=w,p=p,h=h))
    results=[]
    specifications=[(Int32(0),0f0)]
    append!(specifications,[(Int32(c.id),a) for c in selected for a in (-.3f0,.3f0)])
    append!(specifications,[(Int32(0),-.3f0),(Int32(0),.3f0)])
    for(producer,amplitude) in specifications
        b=deepcopy((w=w,p=p,h=h))
        pulse=Pulse(b.h,zeros(Float32,Int(b.p.cells.geometry.capacity)),producer,amplitude,b.h.steps+50)
        samples=[];previous=0.
        for lag in lags
            V.L.advance!(b.w,b.p,pulse,lag-previous;dt=.02);previous=lag
            o=observe(b.w,b.p,b.h,origin)
            returns=V.returns(initial,o,V.N.birth_events(b.h),checkpoint)
            push!(samples,(;lag,observation=o,returns))
        end
        push!(results,(;producer,amplitude,collective=producer==0 && amplitude!=0,samples))
        println((;producer,amplitude));flush(stdout)
    end
    report=(;schema="temporal-single-resident-authority-v1",checkpoint,replay_error,selection="minimum, median, maximum native raw balance, sorted before interventions; stable ID tie-break",selected,initial,results,
        pulse_time=1.,amplitude=.3,noise_seed=6201,recipe=f.cfg,source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
        scope="exact CPU replay of completed history checkpoint; frozen coefficients and critic, unchanged peer policy with common future innovations; source clocks/masks/remaining budgets preserved; collective interventions separately labeled")
    open(io->JSON3.write(io,report),joinpath(out,"report.json"),"w");report
end
end
