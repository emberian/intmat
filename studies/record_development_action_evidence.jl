#!/usr/bin/env julia
# Read-only observation of the immutable developmental source epoch.
length(ARGS)==3||error("SOURCE_EPOCH PROFILE FRESH_OUTPUT")
epoch,profile,out=abspath.(ARGS)
ispath(out)&&error("fresh output required");mkpath(out)
include(joinpath(epoch,"studies/native_development_living.jl"))
const W=NativeDevelopmentLiving
const V=W.V
f=V.N.fixture(profile);w,p=f.w,f.p
h=W.sidecar(p;reactive=false,rate=.1,learning=true,seed=6201,horizon=4.,act_time=24.,inheritance_fraction=0.,
    (Symbol(k)=>v for(k,v) in f.cfg["temporal_policy"])...)
origin=V.C.A.audit(p);records=Any[];graphs=Any[]
for tick in 1:24
    V.L.advance!(w,p,h,1.;dt=.02)
    z=h.state;s=p.cells.state;q=h.options
    if h.steps>=q.horizon_stride
        slot=mod((h.steps-q.horizon_stride)÷q.decision_stride,q.queue_length)+1
        for i in eachindex(s.active)
            if z.transaction[i]==1
                push!(records,(;id=Int(s.ids[i]),time=w.time,decision_time=w.time-4.,eligible=z.eligible[i]==1,
                    context=Float64.(vec(z.queue_features[i,2:3,slot])),innovation=Float64(z.queue_score[i,2,slot]*q.sigma*q.sigma),
                    target=z.last_reward[i],baseline=z.prediction[i],residual=z.prediction_error[i]))
            end
        end
    end
    edges=Tuple{Int,Int}[]
    for i in eachindex(s.active),k in axes(s.bonds,2)
        j=Int(s.bonds[i,k])
        if s.active[i]==1 && 0<j<=length(s.active) && s.active[j]==1 && s.ids[j]==s.bond_ids[i,k] && s.ids[i]<s.ids[j]
            push!(edges,(Int(s.ids[i]),Int(s.ids[j])))
        end
    end
    push!(graphs,(;time=w.time,edges))
    tick%8==0&&(println((;time=w.time,records=length(records)));flush(stdout))
end
report=(;schema="development-action-evidence-shadow-v1",records,graphs,final=W.observation(w,p,h,origin),
    source=V.Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=V.Interior.DEVICE_HABITAT_SOURCE_PIN,
    actor_updates=sum(h.state.actor_updates),learning_work=sum(h.state.paid),clipping=sum(h.state.clip_total),
    scope="read-only original developmental history replay; recorded actual matured pre-action features, innovations and original targets; no shared policy installed")
open(io->V.JSON3.write(io,report),joinpath(out,"report.json"),"w")
