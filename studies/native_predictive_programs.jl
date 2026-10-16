#!/usr/bin/env julia
"""Paid local prediction driving authored, stateful material programs."""
module NativePredictivePrograms
include("native_living_construction.jl")
using .NativeLivingConstruction: Interior,KernelAbstractions,TOML,SHA,JSON3
using KernelAbstractions
const N=NativeLivingConstruction
const L=N.L
@inline receptor(x)=Float32(x)/(1f0+Float32(x))
@inline concentration(x)=clamp(x,0f0,.98f0)/max(.02f0,1f0-clamp(x,0f0,.98f0))
@inline program_action(program)=program==0 ? (0f0,.5f0) : program==1 ? (.3f0,.9f0) : program==2 ? (.3f0,.1f0) : (.3f0,.5f0)
@inline predictor_feature(m,i,c)=c==1 ? 1f0 : c<=3 ? m.obs[i,c-1] : m.fast[i,c-3]-m.slow[i,c-3]
@inline proposal(old,x,error,norm,rate)=clamp(old+clamp(rate*error*x/norm,-.02f0,.02f0),-1f0,1f0)
@inline function program_score(program,a,b,energy,mass,rawa,rawb,ea,eb,q,hq,opt)
    investment,allocation=program_action(program)
    amount=min(opt.dwell*hq.secretion_rate*mass*investment,2min(a,b),energy/max(hq.secretion_work,1f-6))
    # Explicit approximate conversion/capture model; no grid counterfactual is observed.
    lifetime=hq.decay_rate==0 ? opt.plan_horizon : -expm1(-hq.decay_rate*opt.plan_horizon)/hq.decay_rate
    benefit_a=opt.capture*amount*allocation*hq.catalytic_rate*lifetime*rawa/(hq.half_saturation+rawa)
    benefit_b=opt.capture*amount*(1-allocation)*hq.catalytic_rate*lifetime*rawb/(hq.half_saturation+rawb)
    # Native uptake throughput and storage headroom bound the optimistic capture.
    gaina=min(benefit_a,opt.plan_horizon*q.uptake_rate*mass*ea^q.enzyme_exponent,max(0f0,mass*q.store_capacity-a))
    gainb=min(benefit_b,opt.plan_horizon*q.uptake_rate*mass*eb^q.enzyme_exponent,max(0f0,mass*q.store_capacity-b))
    # Balanced precursor capacity proxy: no reward for merely depositing catalyst.
    2min(a-amount/2+gaina,b-amount/2+gainb)-2min(a,b)-amount*hq.secretion_work/(q.efficiency*(q.potentials[1]+q.potentials[2])/2)
end
@kernel function sense!(s,feed,food,enzyme,descriptor,status,g,m,step,stride)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1 && (step%stride==0 || m.ids[i]!=s.ids[i])
        for k in 1:2
            raw=0f0;free=0f0;cat=0f0
            for (site,weight) in Interior.device_support_edges(descriptor,i)
                w=weight/g.volume
                raw+=w*Float32(feed[site,k]);free+=w*Float32(food[site,k]);cat+=w*reinterpret(Float32,enzyme[site,k])
            end
            m.obs[i,k]=receptor(raw);m.obs[i,k+2]=receptor(free);m.obs[i,k+4]=receptor(cat)
        end
    end
end
@kernel function predict_choose!(s,status,q,hq,z,m,opt,step)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1 && (step%opt.stride==0 || m.ids[i]!=s.ids[i])
        newborn=m.ids[i]!=s.ids[i]
        if newborn
            m.ids[i]=s.ids[i];z.ids[i]=s.ids[i];z.birth_pending[i]=0
            m.program[i]=0;m.phase[i]=0;m.since[i]=step;m.decisions[i]=0;m.updates[i]=0;m.work[i]=0.;m.valid[i]=0
            m.from[i]=-1;m.target[i]=-1;m.error_from[i]=-1;m.error_target[i]=-1
            for k in 1:2
                m.fast[i,k]=m.obs[i,k];m.slow[i,k]=m.obs[i,k];m.previous[i,k]=m.obs[i,k]
                m.funded[i,k]=0.;m.error[i,k]=0f0
                for c in 1:5;m.weights[i,c,k]=0f0;end
            end
        end
        norm=1f-4
        for c in 1:5;norm+=m.features[i,c]^2;end
        cost=0.;errors=(m.obs[i,1]-m.previous[i,1],m.obs[i,2]-m.previous[i,2])
        if m.valid[i]==1 && step==m.target[i]
            for k in 1:2
                pred=0f0
                for c in 1:5;pred+=m.weights[i,c,k]*m.features[i,c];end
                error=errors[k]-pred;m.error[i,k]=m.obs[i,k]-m.forecast[i,k]
                if opt.learn
                    for c in 1:5
                        old=m.weights[i,c,k];next=proposal(old,m.features[i,c],error,norm,opt.rate)
                        cost+=opt.kappa*abs(Float64(next)-Float64(old))
                    end
                end
            end
            m.error_from[i]=m.from[i];m.error_target[i]=step
            request=cost+2eps(s.energy[i])
            if cost>0 && s.energy[i]>=request
                nextenergy,paid=Interior.dc_energy_take(s.energy[i],request)
                if paid>=cost
                    s.energy[i]=nextenergy;Interior.dc_ledger_add!(s.ledger,s.ledger_low,i,1,paid)
                    m.work[i]+=paid;m.updates[i]+=1
                    for k in 1:2
                        pred=0f0
                        for c in 1:5;pred+=m.weights[i,c,k]*m.features[i,c];end
                        error=errors[k]-pred
                        for c in 1:5;m.weights[i,c,k]=proposal(m.weights[i,c,k],m.features[i,c],error,norm,opt.rate);end
                    end
                end
            end
        end
        # Receptor/filter state follows native controller convention; learned marks pay turnover.
        interval=Float32(step-m.last_step[i])*opt.dt
        for k in 1:2
            m.fast[i,k]+=-expm1(-interval/.4f0)*(m.obs[i,k]-m.fast[i,k])
            m.slow[i,k]+=-expm1(-interval/2f0)*(m.obs[i,k]-m.slow[i,k])
        end
        for c in 1:5;m.features[i,c]=predictor_feature(m,i,c);end
        for k in 1:2
            delta=0f0
            for c in 1:5;delta+=m.weights[i,c,k]*m.features[i,c];end
            m.previous[i,k]=m.obs[i,k];m.forecast[i,k]=clamp(m.obs[i,k]+delta,0f0,.98f0)
        end
        m.last_step[i]=step;m.from[i]=step;m.target[i]=(step÷opt.stride+1)*opt.stride;m.valid[i]=1;m.decisions[i]+=1
        # A program may inspect every observation but must finish its minimum dwell.
        if m.phase[i]==1 && step-m.since[i]>=opt.dwell_stride
            m.phase[i]=2;m.since[i]=step
        elseif m.phase[i]==2 && step-m.since[i]>=opt.tail_stride
            m.phase[i]=0;m.program[i]=0;m.since[i]=step
        end
        if newborn || m.phase[i]==0
            ra=concentration(opt.predictive ? m.forecast[i,1] : m.obs[i,1])
            rb=concentration(opt.predictive ? m.forecast[i,2] : m.obs[i,2])
            best=0;score=0f0
            for program in 1:3
                candidate=program_score(program,Float32(s.metabolites[i,1]),Float32(s.metabolites[i,2]),Float32(s.energy[i]),s.mass[i],ra,rb,s.enzymes[i,1],s.enzymes[i,2],q,hq,opt)
                m.scores[i,program+1]=candidate
                if candidate>score+opt.margin;best=program;score=candidate;end
            end
            opt.forced>=0 && (best=opt.forced)
            if best!=m.program[i];m.since[i]=step;end
            m.program[i]=best;m.phase[i]=best==0 ? 0 : 1
        end
    end
end
@kernel function actuate!(s,status,a,z,m,hq,opt,dt)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1
        investment,allocation=program_action(m.phase[i]==1 ? m.program[i] : 0)
        Interior.dh_set_requests!(a,s,i,hq,dt,investment,allocation)
        s.actions[i,11]=2opt.reserve-1
        for k in 1:2;z.requested[i,k]+=Float64(a.requests[i,k]);end
    end
end
@kernel function funded!(s,status,a,m)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1
        for k in 1:2;m.funded[i,k]+=Float64(a.last_secretion[i,k]);end
    end
end
struct ProgramHabitat{H,M,O}
    base::H
    model::M
    program_options::O
end
Base.getproperty(h::ProgramHabitat,k::Symbol)=k in (:base,:model,:program_options) ? getfield(h,k) : getproperty(getfield(h,:base),k)
Base.setproperty!(h::ProgramHabitat,k::Symbol,v)=setproperty!(getfield(h,:base),k,v)
function sidecar(p;predictive=true,learn=true,rate=.2,kappa=.001,observation_time=.2,dwell=.4,tail=.8,plan_horizon=2.,capture=.25,reserve=.02,phase_time=8.,forced=-1,dt=.02)
    all(isfinite,(rate,kappa,observation_time,dwell,tail,plan_horizon,capture,reserve,phase_time,dt)) || error("nonfinite program option")
    rate>=0 && kappa>=0 && 0<=capture<=1 && 0<=reserve<=1 && dt>0 && tail>0 && plan_horizon>0 && phase_time>0 && forced in -1:3 || error("invalid program option")
    b=N.sidecar(p;dt,decision_time=observation_time,learning=false,critic=false,sigma=0.,rate=0.,inheritance_fraction=0.,horizon=2.,act_time=32.)
    b.options=merge(b.options,(phase_stride=round(Int,phase_time/dt),))
    n=Int(p.cells.geometry.capacity)
    alloc(T,d...)=fill!(similar(p.cells.state.inputs,T,d...),zero(T))
    m=(ids=alloc(Int32,n),obs=alloc(Float32,n,6),fast=alloc(Float32,n,2),slow=alloc(Float32,n,2),previous=alloc(Float32,n,2),forecast=alloc(Float32,n,2),error=alloc(Float32,n,2),
        weights=alloc(Float32,n,5,2),features=alloc(Float32,n,5),scores=alloc(Float32,n,4),funded=alloc(Float64,n,2),work=alloc(Float64,n),
        program=alloc(Int32,n),phase=alloc(Int32,n),since=alloc(Int32,n),decisions=alloc(Int32,n),updates=alloc(Int32,n),valid=alloc(Int32,n),last_step=alloc(Int32,n),
        from=alloc(Int32,n),target=alloc(Int32,n),error_from=alloc(Int32,n),error_target=alloc(Int32,n))
    stride=round(Int,observation_time/dt);dwell_stride=round(Int,dwell/dt);tail_stride=round(Int,tail/dt)
    stride>0&&isapprox(stride*dt,observation_time)&&dwell_stride>=stride && isapprox(dwell_stride*dt,dwell) && isapprox(tail_stride*dt,tail) && dwell_stride%stride==0 && tail_stride%stride==0 || error("invalid program clocks")
    opts=(;predictive,learn,rate=Float32(rate),kappa=Float64(kappa),dt=Float32(dt),stride,dwell_stride,tail_stride,dwell=Float32(dwell),plan_horizon=Float32(plan_horizon),capture=Float32(capture),reserve=Float32(reserve),margin=1f-6,forced)
    ProgramHabitat(b,m,opts)
end
Interior.device_anabolic_argument(h::ProgramHabitat)=Interior.device_anabolic_argument(h.base)
Interior.device_anabolic_begin!(h::ProgramHabitat,w,p,dt)=Interior.device_anabolic_begin!(h.base,w,p,dt)
Interior.device_anabolic_inherit!(h::ProgramHabitat,w,p,step)=Interior.device_anabolic_inherit!(h.base,w,p,step)
function Interior.device_anabolic_sense!(h::ProgramHabitat,w,p,dt)
    Interior.device_anabolic_sense!(h.habitat,w,p,dt)
    a=h.habitat.arrays;s=p.state;o=h.program_options
    L.policy_launch!(h,sense!,s,h.habitat.feedstock.arrays.fields,p.environment.fields,a.enzyme_bits,Interior.device_support_descriptor(p),p.environment.status,p.geometry,h.model,h.steps,o.stride;ndrange=Int(p.geometry.capacity))
    L.policy_launch!(h,predict_choose!,s,p.environment.status,p.parameters,h.habitat.parameters,h.state,h.model,o,h.steps;ndrange=Int(p.geometry.capacity))
end
function Interior.device_anabolic_prepare!(h::ProgramHabitat,w,p,dt)
    Interior.device_anabolic_prepare!(h.habitat,w,p,dt)
    L.policy_launch!(h,actuate!,p.state,p.environment.status,h.habitat.arrays,h.state,h.model,h.habitat.parameters,h.program_options,dt;ndrange=Int(p.geometry.capacity))
end
function Interior.device_anabolic_finish!(h::ProgramHabitat,w,p,dt)
    Interior.device_anabolic_finish!(h.base,w,p,dt)
    L.policy_launch!(h,funded!,p.state,p.environment.status,h.habitat.arrays,h.model;ndrange=Int(p.geometry.capacity))
end
@inline observed(m,s,i)=m.ids[i]==s.ids[i] && m.valid[i]==1
function observation(w,p,h)
    s=map(Array,p.cells.state);m=map(Array,h.model);z=map(Array,h.state);o=h.program_options
    cells=[(;id=s.ids[i],parent_id=s.parents[i],position=vec(s.position[i,:]),mass=s.mass[i],energy=s.energy[i],stores=vec(s.metabolites[i,:]),
        program=m.program[i],phase=m.phase[i],sensed_step=observed(m,s,i) ? m.last_step[i] : nothing,phase_since=m.since[i]*o.dt,decision_count=m.decisions[i],raw=observed(m,s,i) ? vec(m.obs[i,1:2]) : nothing,food=observed(m,s,i) ? vec(m.obs[i,3:4]) : nothing,catalyst=observed(m,s,i) ? vec(m.obs[i,5:6]) : nothing,
        fast_raw=vec(m.fast[i,:]),slow_raw=vec(m.slow[i,:]),raw_next=observed(m,s,i) ? vec(m.forecast[i,:]) : nothing,prediction_from_step=observed(m,s,i) ? m.from[i] : nothing,prediction_for_step=observed(m,s,i) ? m.target[i] : nothing,
        prediction_error=m.error_target[i]>=0 ? vec(m.error[i,:]) : nothing,error_from_step=m.error_from[i]>=0 ? m.error_from[i] : nothing,error_target_step=m.error_target[i]>=0 ? m.error_target[i] : nothing,prediction_updates=m.updates[i],prediction_work=m.work[i],
        investment=program_action(m.phase[i]==1 ? m.program[i] : 0)[1],allocation=program_action(m.phase[i]==1 ? m.program[i] : 0)[2],reserve=o.reserve,funded=vec(m.funded[i,:]),requested=vec(z.requested[i,:]),scores=vec(m.scores[i,:])) for i in eachindex(s.active) if s.active[i]==1]
    (;time=w.time,step=h.steps,units="raw/food/catalyst are concentration/(1+concentration); phase_since is model time; prediction stamps are physics steps",cells)
end
end
