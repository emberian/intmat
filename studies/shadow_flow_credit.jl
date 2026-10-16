"""Record causal flow baselines beside an unchanged native learner."""
module ShadowFlowCredit
include(get(ENV,"INTMAT_TEMPORAL_SOURCE",joinpath(@__DIR__,"native_temporal_living.jl")))
using .NativeTemporalLiving: Interior, JSON3
using Statistics
const V=NativeTemporalLiving
mutable struct Shadow{H}
    inner::H
    values::Dict{Int,Float64}
    clocks::Dict{Int,Int}
    flow::Dict{Int,Float64}
    gifts::Dict{Int,Float64}
    gift_ids::Set{Int}
    pending::Dict{Tuple{Int,Int},Any}
    rows::Vector{Any}
end
Shadow(h)=Shadow(h,Dict{Int,Float64}(),Dict{Int,Int}(),Dict{Int,Float64}(),Dict{Int,Float64}(),Set{Int}(),Dict{Tuple{Int,Int},Any}(),Any[])
Interior.device_anabolic_argument(h::Shadow)=Interior.device_anabolic_argument(h.inner)
Interior.device_anabolic_begin!(h::Shadow,w,p,dt)=Interior.device_anabolic_begin!(h.inner,w,p,dt)
Interior.device_anabolic_prepare!(h::Shadow,w,p,dt)=Interior.device_anabolic_prepare!(h.inner,w,p,dt)
function Interior.device_anabolic_inherit!(h::Shadow,w,p,step)
    Interior.device_anabolic_inherit!(h.inner,w,p,step)
    e=h.inner.events
    for k in eachindex(e.child_id)
        id=Int(e.child_id[k])
        if id!=0 && !(id in h.gift_ids)
            parent=Int(e.parent_id[k]);h.gifts[parent]=get(h.gifts,parent,0.)+e.gift[k];push!(h.gift_ids,id)
        end
    end
end
function Interior.device_anabolic_sense!(h::Shadow,w,p,dt)
    Interior.device_anabolic_sense!(h.inner,w,p,dt)
    a=h.inner;s=p.state;q=a.options;step=a.steps
    for i in eachindex(s.active)
        s.active[i]==1 || continue
        id=Int(s.ids[i]);new=!haskey(h.values,id)
        if new || step%q.observation_stride==0
            value=V.L.local_value(s,p.parameters,i)+get(h.gifts,id,0.)
            if new
                h.flow[id]=0.
            else
                interval=(step-h.clocks[id])*Float64(q.dt)
                alpha=-expm1(-interval/1.)
                h.flow[id]=(1-alpha)*h.flow[id]+alpha*(value-h.values[id])/interval
            end
            h.values[id]=value;h.clocks[id]=step
        end
        if step%q.decision_stride==0 && step<q.act_stride
            slot=mod(step÷q.decision_stride,q.queue_length)+1
            a.state.queue_ids[i,slot]==id || error("decision identity not recorded")
            h.pending[(id,step)]=(;baseline=q.horizon_stride*Float64(q.dt)*h.flow[id],
                original=a.state.queue_prediction[i,slot],score=a.state.queue_score[i,2,slot],
                features=vec(copy(a.state.queue_features[i,:,slot])),ready=a.state.queue_ready[i,slot])
        end
    end
end
function Interior.device_anabolic_finish!(h::Shadow,w,p,dt)
    Interior.device_anabolic_finish!(h.inner,w,p,dt)
    a=h.inner;q=a.options;step=a.steps
    if step>=q.horizon_stride && step%q.decision_stride==0
        for i in eachindex(p.state.active)
            if a.state.transaction[i]==1
                id=Int(p.state.ids[i]);key=(id,step-q.horizon_stride)
                haskey(h.pending,key)||error("matured decision missing shadow")
                saved=pop!(h.pending,key);target=a.state.last_reward[i]
                saved.original==a.state.prediction[i] || error("original baseline changed")
                push!(h.rows,(;id,origin_step=key[2],maturity_step=step,target,flow_baseline=saved.baseline,
                    original_baseline=saved.original,score=saved.score,features=saved.features,ready=saved.ready,
                    original_residual=a.state.prediction_error[i],flow_residual=target-saved.baseline))
            end
        end
    end
end
function metrics(rows;rate=.1)
    rms(v)=sqrt(mean(abs2,v))
    base=[r.original_residual for r in rows];flow=[r.flow_residual for r in rows]
    moments(err)=mean(sum((Float64(x)*r.score*e)^2 for x in r.features) for (r,e) in zip(rows,err))
    clips(err)=sum(abs(rate*Float64(x)*r.score*e)>.01 for (r,e) in zip(rows,err) for x in r.features)/sum(length(r.features) for r in rows)
    (;count=length(rows),target_rms=rms([r.target for r in rows]),original_rms=rms(base),flow_rms=rms(flow),
        original_score_second_moment=moments(base),flow_score_second_moment=moments(flow),
        original_clip_fraction=clips(base),flow_clip_fraction=clips(flow))
end
function run(profile,out;duration=24.)
    ispath(out)&&error("fresh output required");mkpath(out)
    f=V.N.fixture(profile);w,p=f.w,f.p
    h=V.sidecar(p;reactive=false,rate=.1,learning=true,seed=6201,horizon=4.,act_time=duration,inheritance_fraction=0.,
        (Symbol(k)=>v for (k,v) in f.cfg["temporal_policy"])...)
    s=Shadow(h);origin=V.C.A.audit(p)
    for _ in 1:round(Int,duration/2)
        V.L.advance!(w,p,s,2.;dt=.02)
        println((;time=w.time,matured=length(s.rows)));flush(stdout)
    end
    eligible=[r for r in s.rows if r.ready==1]
    report=(;schema="shadow-flow-credit-v1",rows=s.rows,all_windows=metrics(s.rows),eligible_windows=metrics(eligible),
        late_eligible=metrics([r for r in eligible if r.maturity_step>=600]),final=V.observation(w,p,h,origin),
        learning_work=sum(h.state.paid),clipping=sum(h.state.clip_total),actor_updates=sum(h.state.actor_updates),
        scope="original learner replayed unchanged; shadow baseline never enters actor or critic; same valid matured ID/cohort targets and actual scores",
        source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN)
    open(io->JSON3.write(io,report),joinpath(out,"report.json"),"w");report
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    ShadowFlowCredit.run(ARGS[1],ARGS[2])
end
