"""Shared external training. Physical simulation is outside this derivative path."""
module RecurrentTraining
include("shared_recurrent_runtime.jl")
include("recurrent_training_objective.jl")
using .SharedRecurrentRuntime: KernelAbstractions
using Random, Statistics
const R=SharedRecurrentRuntime
const O=RecurrentTrainingObjective
export initialize, trainer, train_step!, epoch_specs

function initialize(;seed=921,D=25,H=64,A=3,P=12)
    theta=R.initialize(MersenneTwister(seed),D,H,A+1+P)
    theta.Wy .*= .05f0
    # Native initial bounded means: investment .15, allocation .5, reserve .035.
    A==3 && (theta.by[3]=Float32(log((.035-.02)/(.25-.035))))
    theta
end

function trainer(theta,B,T;rate=.001,sigma=.3,discount=.995,prediction_weight=.1,
        lower=Float32[0,0,.02],upper=Float32[.3,1,.25],slope=fill(.5f0,3),
        prediction_scale=ones(Float32,12))
    all(isfinite,Float32.((rate,sigma,discount,prediction_weight)))&&rate>0&&Float32(sigma)>0&&0<=discount<=1&&prediction_weight>=0||throw(ArgumentError("finite valid training options"))
    length(lower)==length(upper)==length(slope)==3&&length(prediction_scale)==size(theta.Wy,1)-4||throw(DimensionMismatch("trainer output schema"))
    all(isfinite,lower)&&all(isfinite,upper)&&all(isfinite,slope)&&all(isfinite,prediction_scale)&&all(upper.>lower)&&all(slope.>0)&&all(prediction_scale.>0)||throw(ArgumentError("finite positive transform/prediction scales"))
    lo,hi,sl,ps=Float32.(lower),Float32.(upper),Float32.(slope),Float32.(prediction_scale)
    all(isfinite,lo)&&all(isfinite,hi)&&all(isfinite,sl)&&all(isfinite,ps)&&all(hi.>lo)&&all(sl.>0)&&all(ps.>0)||throw(ArgumentError("scales must remain finite positive in Float32"))
    cache=R.allocate(theta,B,T)
    backend=KernelAbstractions.get_backend(theta.Wx)
    adapt(v)=begin a=KernelAbstractions.zeros(backend,Float32,length(v));copyto!(a,v);a end
    (;cache,objective=O.workspace(cache.y),optimizer=R.optimizer(theta;rate,clip=1.),
      sigma,discount,prediction_weight,lower=adapt(lower),upper=adapt(upper),slope=adapt(slope),
      prediction_scale=adapt(prediction_scale))
end

"""A fresh rollout supplies incoming memory masks and outgoing return masks.

Required arrays: x,h0,memory_continuation,z,noise,reward,return_continuation,
valid,bootstrap,prediction,predvalid; count is the number of valid decisions.
All arrays share theta's backend. Predictions are passive on-policy forecasts.
"""
function train_step!(theta,t,batch;update=true)
    R.forward!(t.cache,theta,batch.x,batch.h0,batch.memory_continuation)
    O.evaluate!(t.objective,t.cache.y,batch.z,batch.noise,batch.reward,
        batch.return_continuation,batch.valid,batch.bootstrap,batch.prediction,batch.predvalid;
        count=batch.count,sigma=t.sigma,discount=t.discount,lower=t.lower,upper=t.upper,
        slope=t.slope,prediction_scale=t.prediction_scale,prediction_weight=t.prediction_weight)
    R.backward!(t.cache,theta,batch.x,batch.h0,batch.memory_continuation,t.objective.gradient)
    update && R.update!(theta,t.cache.g,t.optimizer)
    R.check_status!(t.optimizer)
    parts=vec(sum(Array(t.objective.parts);dims=(2,3)))
    (;actor_loss=parts[1],value_loss=parts[2],prediction_loss=parts[3],entropy=parts[4],
      loss=parts[1]+parts[2]+parts[3]-.01*parts[4],valid_decisions=batch.count,
      reward_sum=sum(Array(batch.reward)),updated=update)
end

"""Declared independent training/held-out environments, with no phase input."""
function epoch_specs()
    (training=[(;seed=211,phase_time=6.),(;seed=213,phase_time=8.),(;seed=215,phase_time=10.)],
     heldout=[(;seed=212,phase_time=7.),(;seed=214,phase_time=9.)])
end
end
