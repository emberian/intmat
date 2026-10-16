"""Causal offline representation gate on real, identity-indexed native sensing."""
module MaterialMemoryTransfer
include("homeostatic_temporal_memory.jl")
using .HomeostaticTemporalMemory: JSON3
using LinearAlgebra, Statistics, SHA
const H=HomeostaticTemporalMemory

function encode(rows)
    n=length(rows);raw=zeros(n,3);fixed=zeros(n,6);adaptive=zeros(n,6);centers=zeros(n,3);gains=zeros(n,3)
    states=Dict{Int,Any}();last=Dict{Int,Int}();births=0
    for (j,row) in enumerate(rows)
        id=Int(row.id);step=Int(row.step);u=Float64.(row.observation)
        if !haskey(states,id)
            states[id]=[H.State(u[k],.25^2/2,0.,0.,Inf,0.) for k in 1:3]
            last[id]=step;births+=row.generation>0
        else
            step>last[id] || error("resident observations must advance strictly")
        end
        interval=(step-last[id])*.02
        for k in 1:3
            s=states[id][k]
            interval>0 && H.observe!(s,u[k],H.options(dt=interval))
            adaptive[j,2k-1]=s.fast;adaptive[j,2k]=s.slow
            centers[j,k]=s.center;gains[j,k]=H.gain(s,H.options())
            fixed[j,2k-1]=row.native_fast[k];fixed[j,2k]=row.native_slow[k]
        end
        raw[j,:].=u;last[id]=step
    end
    (;raw,fixed,adaptive,centers,gains,births,states)
end
function indices(rows;lag_steps=100)
    lookup=Dict((Int(r.id),Int(r.step))=>i for (i,r) in enumerate(rows))
    pairs=[(i,lookup[(Int(r.id),Int(r.step)+lag_steps)]) for (i,r) in enumerate(rows)
        if haskey(lookup,(Int(r.id),Int(r.step)+lag_steps)) && Int(r.step)%10==0]
    # Training labels end before16; held-out labels lie22..32.
    train=[p for p in pairs if 100<=Int(rows[p[1]].step)<700]
    test=[p for p in pairs if 1000<=Int(rows[p[1]].step)<1500]
    (;train,test)
end
function fit_predict(X,Y,train,test;lambda=.01)
    it=first.(train);jt=last.(train);ie=first.(test);je=last.(test)
    center=vec(mean(X[it,:];dims=1));scale=vec(std(X[it,:];dims=1,corrected=false));scale=max.(scale,1e-6)
    A=hcat(ones(length(it)),(X[it,:].-center')./scale');B=hcat(ones(length(ie)),(X[ie,:].-center')./scale')
    penalty=Matrix{Float64}(I,size(A,2),size(A,2));penalty[1,1]=0.
    weights=(A'A/length(it)+lambda*penalty)\(A'Y[jt,:]/length(it))
    prediction=B*weights
    (;prediction,weights,center,scale,mse=vec(mean((prediction.-Y[je,:]).^2;dims=1)))
end
function summary(x)
    (;minimum=minimum(x),q25=quantile(x,.25),median=median(x),q75=quantile(x,.75),maximum=maximum(x),std=std(x;corrected=false))
end
function run(input,out;manifest=nothing)
    isfile(out)&&error("fresh output required")
    d=JSON3.read(read(input,String));rows=d.rows;e=encode(rows);p=indices(rows)
    isempty(p.train) && error("no training pairs");isempty(p.test) && error("no held-out pairs")
    trainids=Set(Int(rows[i].id) for (i,_) in p.train)
    matrices=(raw_only=e.raw,raw_fixed=hcat(e.raw,e.fixed),raw_homeostatic=hcat(e.raw,e.adaptive),homeostatic_only=e.adaptive)
    models=Dict{String,Any}();predictions=Dict{String,Any}()
    for (name,X) in pairs(matrices)
        f=fit_predict(X,e.raw,p.train,p.test)
        groups=[(;name=label,count=count(mask),mse=any(mask) ? vec(mean((f.prediction[mask,:].-e.raw[last.(p.test)[mask],:]).^2;dims=1)) : nothing)
            for (label,mask) in (("seen_training_ids",[Int(rows[i].id) in trainids for (i,_) in p.test]),
                ("new_training_ids",[!(Int(rows[i].id) in trainids) for (i,_) in p.test]),
                ("crosses_source_switch",[Int(rows[i].step)÷400 != Int(rows[j].step)÷400 for (i,j) in p.test]),
                ("within_source_phase",[Int(rows[i].step)÷400 == Int(rows[j].step)÷400 for (i,j) in p.test]))]
        models[string(name)]=(;mse=f.mse,groups,weights=f.weights,training_center=f.center,training_scale=f.scale)
        predictions[string(name)]=f.prediction
    end
    ranges=[(;phase,channels=[summary(e.raw[[phase*400<=Int(r.step)<(phase+1)*400 for r in rows],k]) for k in 1:3],
        homeostatic_gains=[summary(e.gains[[phase*400<=Int(r.step)<(phase+1)*400 for r in rows],k]) for k in 1:3]) for phase in 0:3]
    gains=[H.gain(s,H.options()) for ss in values(e.states) for s in ss]
    report=(;schema="native-material-memory-transfer-v1",channels=d.channels,models,ranges,
        persistence_mse=vec(mean((e.raw[first.(p.test),:].-e.raw[last.(p.test),:]).^2;dims=1)),
        training_pairs=length(p.train),test_pairs=length(p.test),ids=length(e.states),newborn_ids=e.births,
        train_time=(2.,14.),train_target_time=(4.,16.),test_time=(20.,30.),test_target_time=(22.,32.),horizon=2.,ridge=.01,
        scope="offline same-ID local observation forecast; readout fitted once on early samples; no target/future input to filters; no reward or native policy changes",
        normalization="training-only affine column scaling for ridge; causal homeostat time20, variance prior .25²/2, floor .01²; newborn first-reading center and zero normalized memory",
        sample_histories=[(;id=r.id,parent=r.parent,generation=r.generation,time=r.time,first=r.first,observation=r.observation,center=vec(e.centers[i,:]),gain=vec(e.gains[i,:]))
            for (i,r) in enumerate(rows) if r.id in (1,13,25,35)],
        heldout_samples=[(;id=rows[i].id,time=rows[i].time,current=rows[i].observation,target=rows[j].observation,
            predictions=Dict(name=>vec(v[k,:]) for (name,v) in predictions)) for (k,(i,j)) in enumerate(p.test) if rows[i].id in (1,13,25,35)],
        final_homeostatic_gain=summary(gains),toy_work=sum(s.work for ss in values(e.states) for s in ss),
        observed_source=d.source,observed_habitat_source=d.habitat_source,observation_sha256=bytes2hex(sha256(read(input))),
        source_manifest=manifest===nothing ? nothing : JSON3.read(read(manifest,String)))
    open(io->JSON3.write(io,report),out,"w");report
end
end
