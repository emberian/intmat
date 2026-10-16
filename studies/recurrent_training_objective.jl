"""External recurrent training objective; no differentiation through physics."""
module RecurrentTrainingObjective
using KernelAbstractions
export returns!, loss_gradient, bounded_action, log_sech2

@inline log_sech2(z)=oftype(z,2)*(log(oftype(z,2))-abs(z)-log1p(exp(-oftype(z,2)*abs(z))))
@inline bounded_action(z,lo,hi)=lo+(hi-lo)*(tanh(z)+1)/2

"""n-step bootstrap. continuation is zero at death/identity change, not time limits."""
function returns!(out,reward,continuation,valid,bootstrap;discount=.995)
    isfinite(discount)&&0<=discount<=1||throw(ArgumentError("discount must lie in [0,1]"))
    size(out)==size(reward)==size(continuation)==size(valid)||throw(DimensionMismatch())
    N,T=size(reward);length(bootstrap)==N||throw(DimensionMismatch())
    for n in 1:N
        r=bootstrap[n]
        for t in T:-1:1
            r=valid[n,t] ? reward[n,t]+discount*continuation[n,t]*r : zero(r)
            out[n,t]=r
        end
    end
    out
end

"""Output-level loss and derivative for a shared recurrent trunk.

y layout: A means, one value, P predicted observations. z is the fixed sampled
pre-transform action; noise is its fixed standard-normal innovation. Targets
and advantages are supplied frozen by caller. Prediction targets are next-step
observations. Entropy is a reparameterized Monte Carlo estimate of the actual
tanh/affine action distribution, not the latent Gaussian's entropy. Keep noise
fixed when finite-differencing; policy logprob holds z fixed instead.
"""
function loss_gradient(y,z,noise,targets,advantages,prediction,valid,predvalid;
        sigma=.3,lower=fill(-1.,size(z,1)),upper=fill(1.,size(z,1)),slope=ones(size(z,1)),
        value_weight=.5,prediction_weight=.1,entropy_weight=.01,
        prediction_scale=ones(size(prediction,1)))
    O,N,T=size(y);A=size(z,1);P=size(prediction,1)
    O==A+1+P||throw(DimensionMismatch("output schema"))
    size(z)==size(noise)==(A,N,T)||throw(DimensionMismatch("actions"))
    size(targets)==size(advantages)==size(valid)==size(predvalid)==(N,T)||throw(DimensionMismatch("masks"))
    size(prediction)==(P,N,T)||throw(DimensionMismatch("prediction"))
    length(lower)==length(upper)==length(slope)==A&&length(prediction_scale)==P||throw(DimensionMismatch("scale schema"))
    all(isfinite,(sigma,value_weight,prediction_weight,entropy_weight))&&sigma>0&&min(value_weight,prediction_weight,entropy_weight)>=0||throw(ArgumentError("finite nonnegative loss weights and positive sigma"))
    all(isfinite,lower)&&all(isfinite,upper)&&all(isfinite,slope)&&all(isfinite,prediction_scale)&&all(upper.>lower)&&all(slope.>0)&&all(prediction_scale.>0)||throw(ArgumentError("finite positive scales"))
    count=sum(valid);count>0||throw(ArgumentError("empty training batch"))
    g=zeros(eltype(y),size(y));actor=0.;value=0.;pred=0.;entropy=0.
    for t in 1:T,n in 1:N
        valid[n,t]||continue
        for a in 1:A
            mu=y[a,n,t];sample=z[a,n,t];scale=(upper[a]-lower[a])/2
            logp=-.5*((sample-mu)/sigma)^2-log(sigma)-.5log(2pi)-log(scale*slope[a])-log_sech2(slope[a]*sample)
            actor-=advantages[n,t]*logp
            u=mu+sigma*noise[a,n,t]
            entropy+=.5log(2pi*exp(1)*sigma^2)+log(scale*slope[a])+log_sech2(slope[a]*u)
            g[a,n,t]=(-advantages[n,t]*(sample-mu)/sigma^2+2entropy_weight*slope[a]*tanh(slope[a]*u))/count
        end
        err=y[A+1,n,t]-targets[n,t]
        value+=value_weight*err^2;g[A+1,n,t]=2value_weight*err/count
        if predvalid[n,t] && P>0
            for p in 1:P
                err=(y[A+1+p,n,t]-prediction[p,n,t])/prediction_scale[p]
                pred+=prediction_weight*err^2/P
                g[A+1+p,n,t]=2prediction_weight*err/(P*prediction_scale[p]*count)
            end
        end
    end
    (;loss=(actor+value+pred-entropy_weight*entropy)/count,gradient=g,
      actor=actor/count,value=value/count,prediction=pred/count,entropy=entropy/count,count)
end

@kernel function return_kernel!(out,reward,continuation,valid,bootstrap,discount,T)
    n=@index(Global,Linear)
    r=bootstrap[n]
    for t in T:-1:1
        r=valid[n,t] ? reward[n,t]+discount*continuation[n,t]*r : zero(r)
        out[n,t]=r
    end
end

@kernel function objective_kernel!(g,parts,y,z,noise,target,prediction,valid,predvalid,q,lower,upper,slope,scales)
    n,t=@index(Global,NTuple)
    for o in axes(y,1);g[o,n,t]=0;end
    for k in 1:4;parts[k,n,t]=0;end
    if valid[n,t]
        advantage=target[n,t]-y[q.A+1,n,t]
        actor=0f0;entropy=0f0;pred=0f0
        for a in 1:q.A
            mu=y[a,n,t];sample=z[a,n,t];scale=(upper[a]-lower[a])/2
            logp=-.5f0*((sample-mu)/q.sigma)^2-log(q.sigma)-.5f0*log(2f0*Float32(pi))-log(scale*slope[a])-log_sech2(slope[a]*sample)
            actor-=advantage*logp
            u=mu+q.sigma*noise[a,n,t]
            entropy+=.5f0*log(2f0*Float32(pi)*exp(1f0)*q.sigma^2)+log(scale*slope[a])+log_sech2(slope[a]*u)
            g[a,n,t]=(-advantage*(sample-mu)/q.sigma^2+2*q.entropy_weight*slope[a]*tanh(slope[a]*u))/q.count
        end
        err=y[q.A+1,n,t]-target[n,t]
        value=q.value_weight*err^2;g[q.A+1,n,t]=2*q.value_weight*err/q.count
        if predvalid[n,t] && q.P>0
            for p in 1:q.P
                err=(y[q.A+1+p,n,t]-prediction[p,n,t])/scales[p]
                pred+=q.prediction_weight*err^2/q.P
                g[q.A+1+p,n,t]=2*q.prediction_weight*err/(q.P*scales[p]*q.count)
            end
        end
        parts[1,n,t]=actor/q.count;parts[2,n,t]=value/q.count
        parts[3,n,t]=pred/q.count;parts[4,n,t]=entropy/q.count
    end
end

"""Preallocated native objective buffers; caller uses the same backend as y."""
function workspace(y)
    O,N,T=size(y)
    (;gradient=similar(y),parts=similar(y,4,N,T),targets=similar(y,N,T))
end

"""Native derivative. Advantages are frozen target minus current behavior value.

continuation here is OUTGOING reward continuation, unlike the RNN's incoming
memory-reset mask. count must count valid decision samples, not capacity slots.
"""
function evaluate!(w,y,z,noise,reward,continuation,valid,bootstrap,prediction,predvalid;
        count,A=size(z,1),sigma=.3,discount=.995,lower,upper,slope,
        prediction_scale,value_weight=.5,prediction_weight=.1,entropy_weight=.01)
    count isa Integer&&0<count<=size(reward,1)*size(reward,2)||throw(ArgumentError("valid decision count out of range"))
    backend=get_backend(y);N,T=size(reward)
    P=size(prediction,1)
    size(y)==(A+1+P,N,T)&&size(z)==size(noise)==(A,N,T)||throw(DimensionMismatch("head/action layout"))
    size(continuation)==size(valid)==size(predvalid)==(N,T)&&length(bootstrap)==N&&size(prediction)==(P,N,T)||throw(DimensionMismatch("rollout layout"))
    size(w.gradient)==size(y)&&size(w.targets)==(N,T)&&size(w.parts)==(4,N,T)||throw(DimensionMismatch("objective workspace"))
    length(lower)==length(upper)==length(slope)==A&&length(prediction_scale)==P||throw(DimensionMismatch("scale layout"))
    all(isfinite,Float32.((sigma,discount,value_weight,prediction_weight,entropy_weight)))&&Float32(sigma)>0&&0<=discount<=1&&min(value_weight,prediction_weight,entropy_weight)>=0||throw(ArgumentError("finite representable objective options"))
    all(a->get_backend(a)==backend,(z,noise,reward,continuation,valid,bootstrap,prediction,predvalid,lower,upper,slope,prediction_scale))||throw(ArgumentError("objective backend mismatch"))
    q=(;A,P=size(prediction,1),count=Float32(count),sigma=Float32(sigma),
       value_weight=Float32(value_weight),prediction_weight=Float32(prediction_weight),entropy_weight=Float32(entropy_weight))
    return_kernel!(backend)(w.targets,reward,continuation,valid,bootstrap,Float32(discount),T;ndrange=N)
    objective_kernel!(backend)(w.gradient,w.parts,y,z,noise,w.targets,prediction,valid,predvalid,q,lower,upper,slope,prediction_scale;ndrange=(N,T))
    w
end
end
