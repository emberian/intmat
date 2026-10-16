#!/usr/bin/env julia
"""Uninterrupted paid action learning through the production anabolic hooks.

Decision and reward clocks are physical time, independent of the physics dt.
No simulator schedule or resource law is copied here.
"""
module NativeLifetimeCredit
using Interior, KernelAbstractions, Random, Statistics, SHA
using Interior: JSON3
include("native_action_credit.jl")

@inline local_value(s,q,i)=Float64(s.energy[i])+Float64(q.efficiency)*(Float64(s.metabolites[i,1])*q.potentials[1]+Float64(s.metabolites[i,2])*q.potentials[2])
@kernel function choose!(s,status,q,z,rows,h,decision,sigma,seed,clock,window_start,critic_mode)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1
        changed=z.ids[i]!=s.ids[i]
        if changed
            z.ids[i]=s.ids[i];z.ready[i]=0;z.baseline[i]=0.;z.previous[i]=local_value(s,q,i)
            z.paid[i]=0.;z.last_reward[i]=0.;z.older_reward[i]=0.;z.prediction[i]=0.
            for c in 1:8;z.critic[i,c]=0f0;z.features[i,c]=0f0;end
            for a in 1:length(rows),c in 1:size(z.learned,2)
                z.learned[i,c,a]=0f0;z.trace[i,c,a]=0f0
            end
        end
        # Prediction uses only observations and completed returns, before any
        # current-window noise is sampled. Features have fixed physical scales.
        if critic_mode && (window_start || changed)
            f1=s.inputs[i,2h+1];f2=s.inputs[i,2h+2]
            z.features[i,7]=f1-z.features[i,2];z.features[i,8]=f2-z.features[i,3]
            z.features[i,1]=1f0;z.features[i,2]=f1;z.features[i,3]=f2
            z.features[i,4]=Float32(local_value(s,q,i)/2)
            z.features[i,5]=Float32(z.last_reward[i]/.1)
            z.features[i,6]=Float32((z.last_reward[i]-z.older_reward[i])/.1)
            prediction=0f0
            for c in 1:8;prediction+=z.critic[i,c]*z.features[i,c];end
            z.prediction[i]=Float64(prediction)
        end
        for a in 1:length(rows)
            row=h+rows[a]
            if decision || changed
                mu=0f0
                for c in 1:size(z.learned,2)
                    mu+=(s.weights[i,c,row]+z.learned[i,c,a])*s.inputs[i,c]
                end
                noise=sigma*Interior.dc_normal(seed,s.ids[i],clock,100+a)
                z.drive[i,a]=mu+noise
                if sigma>0
                    for c in 1:size(z.learned,2)
                        z.trace[i,c,a]+=s.inputs[i,c]*noise/(sigma*sigma)
                    end
                end
            end
            # Hold the sampled pre-tanh drive between decisions. Learned marks
            # live separately; native hidden-row marks continue their own law.
            for c in 1:size(z.learned,2)
                s.marks[i,c,row]=-s.weights[i,c,row]+(c==size(z.learned,2) ? z.drive[i,a] : 0f0)
            end
        end
    end
end
@inline actor_proposal(old,trace,advantage,rate,bound)=clamp(old+clamp(rate*Float32(advantage)*trace,-.01f0,.01f0),-bound,bound)
@inline critic_proposal(old,x,step)=clamp(old+clamp(step*x,-.05f0,.05f0),-.3f0,.3f0)
@kernel function reward!(s,status,q,z,rate,bound,kappa,learning,critic_mode,critic_learning)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1 && z.ids[i]==s.ids[i]
        value=local_value(s,q,i);r=value-z.previous[i]
        prediction=critic_mode ? z.prediction[i] : z.baseline[i]
        advantage=r-prediction;cost=0.;canlearn=learning && z.ready[i]==1
        z.prediction_error[i]=advantage;z.ema_error[i]=r-z.baseline[i]
        moment=0.;clipped=Int32(0)
        for a in 1:size(z.learned,3),c in 1:size(z.learned,2)
            moment+=Float64(z.trace[i,c,a])^2
            clipped+=abs(rate*Float32(advantage)*z.trace[i,c,a])>.01f0
            if canlearn
                old=z.learned[i,c,a];proposed=actor_proposal(old,z.trace[i,c,a],advantage,rate,bound)
                cost+=Float64(kappa)*abs(Float64(proposed)-Float64(old))
            end
        end
        z.credit_moment[i]=advantage^2*moment;z.clipped[i]=clipped
        norm=1f-4
        for c in 1:8;norm+=z.features[i,c]^2;end
        criticstep=.3f0*Float32(r-z.prediction[i])/norm
        if critic_mode && critic_learning
            for c in 1:8
                old=z.critic[i,c];proposed=critic_proposal(old,z.features[i,c],criticstep)
                cost+=Float64(kappa)*abs(Float64(proposed)-Float64(old))
            end
        end
        # Actor and critic changes are one affordable material transaction.
        request=cost+2eps(s.energy[i])
        if cost>0 && s.energy[i]>=request
            nextenergy,paid=Interior.dc_energy_take(s.energy[i],request)
            if paid>=cost
                s.energy[i]=nextenergy;Interior.dc_ledger_add!(s.ledger,s.ledger_low,i,1,paid)
                z.paid[i]+=paid
                if canlearn
                    for a in 1:size(z.learned,3),c in 1:size(z.learned,2)
                        old=z.learned[i,c,a]
                        z.learned[i,c,a]=actor_proposal(old,z.trace[i,c,a],advantage,rate,bound)
                    end
                end
                if critic_mode && critic_learning
                    for c in 1:8
                        old=z.critic[i,c]
                        z.critic[i,c]=critic_proposal(old,z.features[i,c],criticstep)
                    end
                end
            end
        end
        # Previous value is before learning work; the next return includes cost.
        z.previous[i]=value;z.baseline[i]=z.ready[i]==1 ? .8z.baseline[i]+.2r : r
        z.ready[i]=1;z.older_reward[i]=z.last_reward[i];z.last_reward[i]=r
        for a in 1:size(z.trace,3),c in 1:size(z.trace,2);z.trace[i,c,a]=0f0;end
    end
end
mutable struct CreditHabitat{H,Z,B,Q}
    habitat::H
    state::Z
    backend::B
    options::Q
    steps::Int
end
function sidecar(p;dt=.02,decision_time=.2,reward_time=2.,rate=.003,sigma=.15,kappa=.001,seed=6201,learning=true,critic=false,critic_learning=critic)
    all(x->isfinite(x)&&x>0,(dt,decision_time,reward_time))||error("positive clocks required")
    for t in (decision_time,reward_time);isapprox(t/dt,round(t/dt);atol=1e-8)||error("clocks must divide physics dt");end
    isapprox(reward_time/decision_time,round(reward_time/decision_time);atol=1e-8)||error("reward clock must divide decision clock")
    s=p.cells.state;C=Int(p.cells.geometry.capacity);h=Int(p.cells.geometry.hidden);F=size(s.inputs,2);A=3
    alloc(T,dims...)=fill!(similar(s.inputs,T,dims...),zero(T))
    z=(learned=alloc(Float32,C,F,A),trace=alloc(Float32,C,F,A),drive=alloc(Float32,C,A),
       ids=alloc(Int32,C),ready=alloc(Int32,C),baseline=alloc(Float64,C),previous=alloc(Float64,C),paid=alloc(Float64,C),last_reward=alloc(Float64,C),older_reward=alloc(Float64,C),
       critic=alloc(Float32,C,8),features=alloc(Float32,C,8),prediction=alloc(Float64,C),
       prediction_error=alloc(Float64,C),ema_error=alloc(Float64,C),credit_moment=alloc(Float64,C),clipped=alloc(Int32,C))
    opts=(dt=Float32(dt),decision_stride=round(Int,decision_time/dt),reward_stride=round(Int,reward_time/dt),
          rate=Float32(rate),sigma=Float32(sigma),kappa=Float64(kappa),seed=UInt32(seed),rows=(4,5,7),hidden=h,learning,critic,critic_learning)
    CreditHabitat(p.habitat,z,p.cells.backend,opts,0)
end
Interior.device_anabolic_argument(h::CreditHabitat)=Interior.device_anabolic_argument(h.habitat)
Interior.device_anabolic_begin!(h::CreditHabitat,w,p,dt)=Interior.device_anabolic_begin!(h.habitat,w,p,dt)
Interior.device_anabolic_prepare!(h::CreditHabitat,w,p,dt)=Interior.device_anabolic_prepare!(h.habitat,w,p,dt)
Interior.device_anabolic_inherit!(h::CreditHabitat,w,p,step)=Interior.device_anabolic_inherit!(h.habitat,w,p,step)
function Interior.device_anabolic_sense!(h::CreditHabitat,w,p,dt)
    dt==h.options.dt||error("sidecar physics clock changed")
    Interior.device_anabolic_sense!(h.habitat,w,p,dt)
    q=h.options
    s=p.state;v=(active=s.active,ids=s.ids,energy=s.energy,metabolites=s.metabolites,weights=s.weights,marks=s.marks,inputs=s.inputs)
    chemistry=(efficiency=p.parameters.efficiency,potentials=p.parameters.potentials)
    choose!(h.backend,64)(v,p.environment.status,chemistry,h.state,q.rows,q.hidden,
        h.steps%q.decision_stride==0,q.sigma,q.seed,UInt32(h.steps÷q.decision_stride+1),h.steps%q.reward_stride==0,q.critic;ndrange=Int(p.geometry.capacity))
end
function Interior.device_anabolic_finish!(h::CreditHabitat,w,p,dt)
    Interior.device_anabolic_finish!(h.habitat,w,p,dt);h.steps+=1
    q=h.options
    if h.steps%q.reward_stride==0
        s=p.state;v=(active=s.active,ids=s.ids,energy=s.energy,metabolites=s.metabolites,ledger=s.ledger,ledger_low=s.ledger_low)
        chemistry=(efficiency=p.parameters.efficiency,potentials=p.parameters.potentials)
        reward!(h.backend,64)(v,p.environment.status,chemistry,h.state,q.rate,.3f0,q.kappa,q.learning,q.critic,q.critic_learning;ndrange=Int(p.geometry.capacity))
    end
end
function advance!(w,p,h,time;dt=.02)
    n=round(Int,time/dt);isapprox(n*dt,time;atol=1e-8)||error("duration must divide dt")
    for _ in 1:n;Interior.device_cellular_step!(w,p.cells,dt;anabolism=h,check_status=false);end
    Interior.device_cellular_status!(p.cells);w.time+=n*dt
end
function observation(w,p,h,origin)
    s=p.cells.state;a=NativeActionCredit.audit(p)
    material=a.material-origin.material-a.supplied+origin.supplied+a.removed-origin.removed
    potential=a.potential-origin.potential+a.work-origin.work+a.heat-origin.heat-a.supplied_potential+origin.supplied_potential+a.removed_potential-origin.removed_potential
    a.finite && abs(material)<1e-8 && abs(potential)<1e-8||error("paid learner accounting failure")
    (;time=w.time,ids=collect(s.ids[s.active.==1]),usable=NativeActionCredit.usable(p),
        learning_work=collect(h.state.paid[1:4]),work=a.work-origin.work,heat=a.heat-origin.heat,
        uptake_a=[Float64(s.ledger[i,3])+Float64(s.ledger_low[i,3]) for i in 1:4],
        uptake_b=[Float64(s.ledger[i,4])+Float64(s.ledger_low[i,4]) for i in 1:4],
        decisions=h.steps÷h.options.decision_stride,reward_windows=h.steps÷h.options.reward_stride,baseline=collect(h.state.baseline[1:4]),
        last_return=collect(h.state.last_reward[1:4]),
        prediction_error=collect(h.state.prediction_error[1:4]),ema_error=collect(h.state.ema_error[1:4]),
        credit_moment=collect(h.state.credit_moment[1:4]),clipped=collect(h.state.clipped[1:4]),critic_max=maximum(abs,h.state.critic),mark_max=maximum(abs,h.state.learned),material,potential)
end
function patch!(p)
    # Half-period displacement of free and raw food. Every species' total and
    # potential are preserved; catalyst/bound material/cell memories stay put.
    before=NativeActionCredit.audit(p);nx,ny,nz=Int.(p.cells.geometry.dims)
    for field in (p.cells.environment.fields,p.habitat.feedstock.arrays.fields)
        old=copy(field)
        for k in 1:2,z in 1:nz,y in 1:ny,x in 1:nx
            i=x+nx*((y-1)+ny*(z-1));j=mod1(x+nx÷2,nx)+nx*((y-1)+ny*(z-1))
            field[i,k]=old[j,k]
        end
    end
    after=NativeActionCredit.audit(p)
    abs(after.material-before.material)<1e-10 && abs(after.potential-before.potential)<1e-10||error("patch adds inventory")
end
function run(profile,out;seed=211,training=24.,test_time=4.,dt=.02,noise_seed=6201,critic=false)
    isapprox(training/2,round(training/2);atol=1e-8) && training>=2 || error("training must contain whole reward windows")
    ispath(out)&&error("fresh output required");mkpath(out)
    write(joinpath(out,"executed-driver.jl"),read(@__FILE__))
    f=NativeActionCredit.fixture(profile,(.1,.5);seed);results=Dict{String,Any}()
    for (arm,rate,learning) in (("credit",.003,true),("no_credit",0.,false),("reversed_credit",-.003,true))
        w=World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed);p=Interior.device_habitat_restore(w,f.checkpoint)
        h=sidecar(p;dt,rate,learning,seed=noise_seed,critic,critic_learning=critic);origin=NativeActionCredit.audit(p);trace=Any[]
        for _ in 1:round(Int,training/2)
            advance!(w,p,h,2.;dt);push!(trace,observation(w,p,h,origin))
            println(arm," time=",w.time," paid=",sum(h.state.paid));flush(stdout)
        end
        checkpoint=Interior.device_habitat_checkpoint(w,p);learned=copy(h.state.learned);branches=Dict{String,Any}()
        for erase in (false,true),patch in (false,true)
            wr=World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed);pr=Interior.device_habitat_restore(wr,checkpoint)
            hr=sidecar(pr;dt,rate=0.,sigma=0.,learning=false)
            # Restore identity before installing learned state, so the first
            # sensing call does not mistake the branch for a newborn controller.
            hr.state.ids.=pr.cells.state.ids
            for i in eachindex(hr.state.previous);hr.state.previous[i]=local_value(pr.cells.state,pr.cells.parameters,i);end
            !erase && (hr.state.learned.=learned)
            patch && patch!(pr)
            base=NativeActionCredit.audit(pr);u0=NativeActionCredit.usable(pr)
            advance!(wr,pr,hr,test_time;dt)
            name=(erase ? "erased" : "learned")*(patch ? "_patch" : "_sham")
            branches[name]=(;return_value=NativeActionCredit.usable(pr).-u0,observation(wr,pr,hr,base)...)
        end
        results[arm]=(;trace,branches,learned=learned[1:4,:,:])
        open(io->JSON3.write(io,results[arm]),joinpath(out,arm*".json"),"w")
        println(arm," complete paid=",sum(h.state.paid)," marks=",maximum(abs,learned));flush(stdout)
    end
    report=(;schema="native-lifetime-credit-v2",seed,noise_seed,critic,training,test_time,dt,decision_time=.2,reward_time=2.,kappa=.001,
        source=f.checkpoint.source,cellular_source=f.checkpoint.cells.source,driver_sha256=bytes2hex(sha256(read(joinpath(out,"executed-driver.jl")))),
        scope="uninterrupted nondividing closed 8^3 material; paid local return credit; same-state frozen learned/erased branches with conserved food displacement; no cooperation claim",results)
    open(io->JSON3.pretty(io,report),joinpath(out,"report.json"),"w")
    report
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    length(ARGS)>=2||error("PROFILE FRESH_OUTPUT [TRAINING=24] [TEST=4]")
    NativeLifetimeCredit.run(ARGS[1],ARGS[2];training=length(ARGS)>2 ? parse(Float64,ARGS[3]) : 24.,test_time=length(ARGS)>3 ? parse(Float64,ARGS[4]) : 4.)
end
