"""Read-only native sensing-hook recorder; run from a frozen source epoch."""
module RecordTemporalMaterial
include("native_temporal_living.jl")
using .NativeTemporalLiving: Interior, JSON3
const V=NativeTemporalLiving
mutable struct Recorder{H}
    inner::H
    rows::Vector{Any}
    seen::Set{Int32}
end
Interior.device_anabolic_argument(h::Recorder)=Interior.device_anabolic_argument(h.inner)
Interior.device_anabolic_begin!(h::Recorder,w,p,dt)=Interior.device_anabolic_begin!(h.inner,w,p,dt)
Interior.device_anabolic_prepare!(h::Recorder,w,p,dt)=Interior.device_anabolic_prepare!(h.inner,w,p,dt)
Interior.device_anabolic_finish!(h::Recorder,w,p,dt)=Interior.device_anabolic_finish!(h.inner,w,p,dt)
Interior.device_anabolic_inherit!(h::Recorder,w,p,step)=Interior.device_anabolic_inherit!(h.inner,w,p,step)
function Interior.device_anabolic_sense!(h::Recorder,w,p,dt)
    Interior.device_anabolic_sense!(h.inner,w,p,dt)
    a=h.inner;s=p.state;step=a.steps
    for i in eachindex(s.active)
        if s.active[i]==1
            id=s.ids[i];newborn=!(id in h.seen)
            if step%a.options.observation_stride==0 || newborn
                push!(h.rows,(;step,time=step*Float64(a.options.dt),id,slot=i,parent=s.parents[i],generation=s.generations[i],first=newborn,
                    observation=[V.C.feature(s,i,a.options.hidden,k) for k in 2:4],
                    native_fast=vec(copy(a.temporal.fast[i,:])),native_slow=vec(copy(a.temporal.slow[i,:])),
                    mass=s.mass[i],energy=s.energy[i],usable=V.L.local_value(s,p.parameters,i),position=vec(copy(s.position[i,:]))))
                push!(h.seen,id)
            end
        end
    end
end
function run(profile,out;duration=32.)
    ispath(out)&&error("fresh output required");mkpath(out)
    f=V.N.fixture(profile);w,p=f.w,f.p
    h=V.sidecar(p;rate=0.,learning=false,critic=false,sigma=0.,horizon=4.,act_time=duration,inheritance_fraction=0.,
        (Symbol(k)=>v for (k,v) in f.cfg["temporal_policy"])...)
    r=Recorder(h,Any[],Set{Int32}());origin=V.C.A.audit(p)
    for tick in 1:round(Int,duration/.2)
        V.L.advance!(w,p,r,.2;dt=.02)
        tick%40==0 && (println((;time=w.time,rows=length(r.rows),residents=count(==(1),p.cells.state.active)));flush(stdout))
    end
    final=V.observation(w,p,h,origin)
    report=(;schema="native-material-observations-v1",duration,dt=.02,observation_time=.2,
        timing="after native sensing and temporal-policy observation, before controller/uptake; step clock, not advance! endpoint; off-cadence first-ID samples included",
        scope="fixed zero learned allocation, no exploration/critic/actor updates; moving and native division remain active",
        channels=("raw_balance","raw_total","catalyst_total"),rows=r.rows,final,events=V.N.birth_events(h),recipe=f.cfg,
        source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN)
    open(io->JSON3.write(io,report),joinpath(out,"observations.json"),"w");report
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    RecordTemporalMaterial.run(ARGS[1],ARGS[2])
end
