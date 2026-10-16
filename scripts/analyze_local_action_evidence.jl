#!/usr/bin/env julia
using JSON3, LinearAlgebra, Statistics
include("../studies/local_action_evidence.jl")
using .LocalActionEvidence
length(ARGS)==3||error("SHADOW_REPORT ORIGINAL_ARM OUTPUT")
r=JSON3.read(read(ARGS[1],String));original=JSON3.read(read(ARGS[2],String))
r.final==last(original.samples)||error("final observation replay differs")
r.actor_updates==original.actor_updates||error("actor update mismatch")
r.learning_work==original.learning_work||error("payment mismatch")
r.clipping==original.clipping||error("clip mismatch")
train=[x for x in r.records if x.eligible && x.time<=16.0001]
held=[x for x in r.records if x.eligible && x.time>=20. && x.id in [c.id for c in original.samples[9].cells]]
byid=Dict{Int,Evidence}()
for x in train
    observe!(get!(Evidence,byid,Int(x.id)),x.innovation,x.context,x.residual)
end
graph=only(g for g in r.graphs if abs(g.time-16)<.001)
neighbors=Dict(i=>Int[] for i in keys(byid))
for e in graph.edges
    a,b=Int.(e)
    haskey(neighbors,a)&&haskey(byid,b)&&push!(neighbors[a],b)
    haskey(neighbors,b)&&haskey(byid,a)&&push!(neighbors[b],a)
end
results=[]
for ridge in (.01,.1,1.),mode in (:own,:bond)
    fits=Dict{Int,Any}()
    for id in keys(byid)
        ids=mode==:own ? [id] : [id;neighbors[id]]
        fits[id]=fit([byid[j] for j in ids],1;intercept_ridge=ridge,slope_ridge=ridge)
    end
    rows=[x for x in held if haskey(fits,Int(x.id))]
    error=mean((x.residual-x.innovation*(fits[Int(x.id)].intercept+dot(fits[Int(x.id)].slope,x.context)))^2 for x in rows)
    baseline=mean(x.residual^2 for x in rows)
    push!(results,(;ridge,mode,n=length(rows),residual_mse=error,zero_response_mse=baseline,
        median_smallest_information=median(minimum(eigvals(f.information))-ridge for f in values(fits)),
        median_slope_norm=median(norm(f.slope) for f in values(fits))))
end
report=(;schema="local-action-evidence-identification-v1",replay_exact=true,training_count=length(train),heldout_count=length(held),
    train_end=16,heldout_start=20,primary_ridge=.1,results,
    scope="fixed time16 bond neighborhood; train original eligible maturities through16 and predict20..24; descriptive response fit, not learned communication or causal benefit")
open(io->JSON3.pretty(io,report),ARGS[3],"w");println(report)
