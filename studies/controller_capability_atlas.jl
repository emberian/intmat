#!/usr/bin/env julia
"""Driven capability assay using the actual native controller/regulation kernels.

Cells are clamped instruments here: no resource, mechanics or lifecycle advance.
This measures readable state, not ecological utility or single-arena throughput.
"""
module ControllerCapabilityAtlas
using Interior, KernelAbstractions, Random, Statistics, LinearAlgebra, JSON3, SHA

function bank(;n=8,hidden=8,seed=211,gain=1.,adaptation=.3,plasticity=.002)
    w=World(Grid(2,2,2;extent=(1,1,1));seed,fields=(),scratch_fields=())
    p=Interior.device_cellular_initialize(w,zeros(3,2n);capacity=2n,hidden,seed,
        support_strategy=:materialized,exchange_radius=.1,radius_bound=1.,
        response_time=.8,adaptation_time=8.,adaptation_strength=adaptation,
        plasticity_rate=plasticity,energy_type=Float64,inventory_type=Float64)
    s=p.state
    # Perturbed twins share every weight and forcing sample, differing only in
    # initial hidden state. Scale recurrent-to-hidden weights, not sensory gain.
    s.weights[:,1:hidden,1:hidden].*=Float32(gain)
    s.weights[n+1:2n,:,:].=s.weights[1:n,:,:]
    fill!(s.hidden,0);fill!(s.adaptation,0);fill!(s.marks,0)
    s.hidden[n+1:2n,:].=1f-3
    fill!(s.support_count,0) # no uptake/scatter; retain production neural update
    w,p
end

function drive!(p,u,dt)
    s=p.state;h=Int(p.geometry.hidden);n=Int(p.geometry.capacity)
    fill!(s.inputs,0)
    s.inputs[:,1:h].=s.hidden
    # Prescribed complementary bounded free-food observations and tonic bias.
    s.inputs[:,2h+1].=Float32(.5+.4u)
    s.inputs[:,2h+2].=Float32(.5-.4u)
    s.inputs[:,2h+20].=1
    Interior.dc_controller!(CPU(),64)(s,p.environment,p.geometry;
        ndrange=(n,h+Interior.CELLULAR_ACTIONS))
    Interior.dc_regulate!(CPU(),64)(s,p.environment,p.parameters,p.geometry,Float32(dt);ndrange=n)
    KernelAbstractions.synchronize(CPU())
    nothing
end

function sequence(;samples=1200,burn=200,seed=41)
    rand(MersenneTwister(seed),samples+burn).*2 .-1
end

function record(cfg,u;burn=200,dt=.05,hold=.4,n=8)
    isapprox(hold/dt,round(hold/dt);atol=1e-8) || error("hold must be a timestep multiple")
    _,p=bank(;n,cfg...);h=Int(p.geometry.hidden);steps=round(Int,hold/dt)
    # Features combine hidden state and slow adaptation, excluding direct input.
    features=zeros(length(u)-burn,2h,n)
    differences=zeros(length(u),n)
    saturation=zeros(length(u)-burn,n)
    for (t,value) in enumerate(u)
        for _ in 1:steps;drive!(p,value,dt);end
        s=p.state
        for i in 1:n
            differences[t,i]=sqrt(sum(abs2,Float64.(s.hidden[i,:]).-s.hidden[n+i,:])/h)
            if t>burn
                features[t-burn,1:h,i].=s.hidden[i,:]
                features[t-burn,h+1:2h,i].=s.adaptation[i,:]
                saturation[t-burn,i]=count(x->abs(x)>.95,s.hidden[i,:])/h
            end
        end
    end
    all(isfinite,features) || error("nonfinite controller features")
    (;features,differences,saturation,inputs=u[burn+1:end])
end

"""Training-only centering/scaling, fixed ridge; no test-set selection."""
function readout(train,test,ytrain,ytest;ridge=1e-3)
    μ=mean(train;dims=1);σ=max.(std(train;dims=1,corrected=false),1e-8)
    a=(train.-μ)./σ;b=(test.-μ)./σ;ym=mean(ytrain)
    weights=(a'a+ridge*size(a,1)*I)\(a'*(ytrain.-ym))
    predicted=b*weights .+ym
    denom=sum(abs2,ytest.-mean(ytest))
    denom>0 || error("constant held-out target")
    1-sum(abs2,predicted.-ytest)/denom
end

function score(train,test;lags=1:12)
    offset=max(maximum(lags),2);idx=(offset+1):length(train.inputs)
    idxt=(offset+1):length(test.inputs)
    target(u,ix,k)=k==0 ? u[ix.-1].*u[ix.-2] : u[ix.-k]
    baseline(u,ix)=hcat(u[ix],u[ix].^2)
    direct=Dict{String,Float64}();cells=Any[]
    for k in (collect(lags)...,0)
        label=k==0 ? "lag1_times_lag2" : "lag$(k)"
        direct[label]=readout(baseline(train.inputs,idx),baseline(test.inputs,idxt),
            target(train.inputs,idx,k),target(test.inputs,idxt,k))
    end
    for i in axes(train.features,3)
        scores=Dict{String,Float64}()
        for k in (collect(lags)...,0)
            label=k==0 ? "lag1_times_lag2" : "lag$(k)"
            scores[label]=readout(train.features[idx,:,i],test.features[idxt,:,i],
                target(train.inputs,idx,k),target(test.inputs,idxt,k))
        end
        # Circularly misalign held-out state with target. Preserve the feature
        # marginal and correlations, deliberately break time correspondence.
        shuffled=circshift(test.features[idxt,:,i],(length(idxt)÷2,0))
        control=readout(train.features[idx,:,i],shuffled,train.inputs[idx.-1],test.inputs[idxt.-1])
        push!(cells,(;founder=i,r2=scores,shifted_lag1_r2=control,
            hidden_saturation=mean(test.saturation[:,i]),
            twin_rms_at_burn=train.differences[size(train.differences,1)-size(train.features,1),i],
            twin_rms_final=train.differences[end,i]))
    end
    (;cells,direct_input_r2=direct)
end

function run(out;quick=false,dt=.05)
    ispath(out) && error("fresh output required");mkpath(out)
    driver=read(@__FILE__);write(joinpath(out,"executed-driver.jl"),driver)
    samples=quick ? 180 : 1200;n=quick ? 2 : 8;burn=200
    train=sequence(;samples,burn,seed=41);test=sequence(;samples,burn,seed=907)
    configs=[(;gain,adaptation,plasticity=.002) for gain in (.6,1.,1.4),adaptation in (0.,.3,.9)]
    results=Any[];started=time()
    for cfg in configs
        a=record(cfg,train;burn,dt,n);b=record(cfg,test;burn,dt,n)
        measured=score(a,b)
        item=(;config=cfg,measured...)
        push!(results,item)
        tag="gain$(cfg.gain)-adapt$(cfg.adaptation)"
        open(io->JSON3.write(io,(;config=cfg,train=a,test=b)),joinpath(out,tag*"-traces.json"),"w")
        println(tag," lag1 median R²=",round(median(c.r2["lag1"] for c in item.cells);digits=3));flush(stdout)
    end
    report=(;schema="interior-controller-capabilities-v1",dt,hold=.4,burn,samples,founders=n,
        train_seed=41,test_seed=907,founder_seed=211,ridge=1e-3,
        source_pin=Interior.DEVICE_CELLULAR_SOURCE_PIN,driver_sha256=bytes2hex(sha256(driver)),
        scope="prescribed-input native controller kernels on CPU; clamped physiology, no ecology; held-out linear readability, not information use or learning benefit",
        elapsed_seconds=time()-started,results)
    open(io->JSON3.pretty(io,report),joinpath(out,"report.json"),"w")
    report
end
end

if abspath(PROGRAM_FILE)==(@__FILE__) || (isdefined(Main,:INTERIOR_RUNNING_STUDY) && Main.INTERIOR_RUNNING_STUDY=="controller_capability_atlas")
    isempty(ARGS) && error("FRESH_OUTPUT [--quick] [--half-step]")
    ControllerCapabilityAtlas.run(ARGS[1];quick="--quick" in ARGS,dt="--half-step" in ARGS ? .025 : .05)
end
