#!/usr/bin/env julia
"""External shared recurrent training in independent finite native environments."""
module TrainNativeRecurrent
include("recurrent_training.jl")
include("native_recurrent_agent.jl")
using .RecurrentTraining: KernelAbstractions
using KernelAbstractions, Random, Serialization, SHA
const RT=RecurrentTraining
const R=RT.R
const A=NativeRecurrentAgent

@kernel function sample!(u,noise,zsave,nsave,y,ids,valid,sigma,seed,clock,t)
    a,n=@index(Global,NTuple)
    e=valid[n]!=0 ? A.Interior.dc_normal(seed,ids[n],clock,400+a) : 0f0
    noise[a,n]=e;u[a,n]=y[a,n]+sigma*e
    zsave[a,n,t]=u[a,n];nsave[a,n,t]=e
end
@kernel function bootstrap!(out,y,valid,valueindex)
    n=@index(Global,Linear)
    out[n]=valid[n]!=0 ? y[valueindex,n] : 0f0
end

function buffers(theta,B,T)
    D,H,O=R.dimensions(theta);P=O-4
    alloc(d...)=fill!(similar(theta.Wx,Float32,d...),0f0)
    bools(d...)=fill!(similar(theta.Wx,Bool,d...),false)
    (;x=alloc(D,B,T),h0=alloc(H,B),memory_continuation=bools(B,T),incoming=bools(B),
      z=alloc(3,B,T),noise=alloc(3,B,T),reward=alloc(B,T),return_continuation=alloc(B,T),
      valid=bools(B,T),bootstrap=alloc(B),prediction=alloc(P,B,T),predvalid=bools(B,T),
      h=alloc(H,B),next_h=alloc(H,B),y=alloc(O,B),u=alloc(3,B),eps=alloc(3,B),
      probe_h=alloc(H,B),probe_y=alloc(O,B))
end

function rollout!(env,theta,b;T=size(b.x,3),sigma=.3,seed=921,clock=0)
    backend=get_backend(theta.Wx);copyto!(b.h0,b.h)
    for t in 1:T
        copyto!(view(b.x,:,:,t),env.x)
        copyto!(view(b.memory_continuation,:,t),env.incoming)
        copyto!(view(b.valid,:,t),env.valid)
        copyto!(b.incoming,env.incoming)
        R.step!(b.next_h,b.y,theta,env.x,b.h,b.incoming)
        sample!(backend)(b.u,b.eps,b.z,b.noise,b.y,env.ids,env.valid,Float32(sigma),UInt32(seed),UInt32(clock+t),t;ndrange=size(b.u))
        A.act!(env,b.u);A.advance!(env)
        copyto!(view(b.reward,:,t),env.reward)
        copyto!(view(b.return_continuation,:,t),env.outgoing)
        copyto!(view(b.prediction,:,:,t),env.target)
        copyto!(view(b.predvalid,:,t),env.predvalid)
        copyto!(b.h,b.next_h)
    end
    copyto!(b.incoming,env.incoming)
    R.step!(b.probe_h,b.probe_y,theta,env.x,b.h,b.incoming)
    bootstrap!(backend)(b.bootstrap,b.probe_y,env.valid,4;ndrange=length(b.bootstrap))
    synchronize(backend)
    merge(b,(;count=count(Array(b.valid))))
end

function recarry!(b,trainer,theta)
    # Recompute the just-collected prefix under current theta; the starting h0
    # still summarizes older history under earlier weights (finite burn-in).
    R.forward!(trainer.cache,theta,b.x,b.h0,b.memory_continuation)
    copyto!(b.h,view(trainer.cache.h,:,:,size(b.x,3)))
end

function evaluate(profile,theta,spec;blocks=8,T=8,seed=991,backend=CPU())
    env=A.make_env(profile;backend,seed=spec.seed,phase_time=spec.phase_time)
    b=buffers(theta,size(env.x,2),T);total=0.;decisions=0
    for block in 1:blocks
        data=rollout!(env,theta,b;T,sigma=0.,seed,clock=(block-1)*T)
        total+=sum(Array(data.reward));decisions+=data.count
    end
    (;environment=spec,blocks,T,total_reward=total,valid_decisions=decisions,
      reward_per_decision=total/max(1,decisions),physical=A.summary(env))
end

function run(profile,out;H=64,T=8,train_blocks=8,eval_blocks=8,rate=.001,seed=921,
        specs=RT.epoch_specs(),backend=CPU())
    ispath(out)&&error("fresh output required");mkpath(out)
    write(joinpath(out,"executed-driver.jl"),read(@__FILE__))
    theta=R.to_backend(backend,RT.initialize(;H,seed));initial=deepcopy(theta);training=Any[]
    opt=nothing
    for (epoch,spec) in enumerate(specs.training)
        env=A.make_env(profile;backend,seed=spec.seed,phase_time=spec.phase_time)
        b=buffers(theta,size(env.x,2),T);trainer=RT.trainer(theta,size(env.x,2),T;rate)
        # Optimizer moments persist across environment resets along with theta.
        if opt!==nothing;trainer=merge(trainer,(;optimizer=opt));end
        opt=trainer.optimizer
        for block in 1:train_blocks
            started=time_ns()
            batch=rollout!(env,theta,b;T,seed=seed+epoch,clock=(block-1)*T)
            stats=RT.train_step!(theta,trainer,batch)
            recarry!(b,trainer,theta);synchronize(get_backend(theta.Wx))
            row=(;epoch,environment=spec,block,stats,elapsed_seconds=(time_ns()-started)/1e9)
            push!(training,row);println(row);flush(stdout)
        end
    end
    heldout=[(;environment=spec,learned=evaluate(profile,theta,spec;blocks=eval_blocks,T,backend),
       initial=evaluate(profile,initial,spec;blocks=eval_blocks,T,backend)) for spec in specs.heldout]
    report=(;schema="native-shared-recurrent-training-v1",H,T,train_blocks,eval_blocks,rate,seed,
      training,heldout,provenance="external shared theta optimizer; private per-ID hidden state; native physical actuation costs; no paid biological parameter adaptation",
      auxiliary="passive on-policy next-observation forecast from shared h, not action-conditioned or counterfactual",
      driver_sha256=bytes2hex(SHA.sha256(read(joinpath(out,"executed-driver.jl")))))
    serialize(joinpath(out,"parameters.bin"),(;theta=map(Array,theta),initial=map(Array,initial)))
    open(io->A.Interior.JSON3.write(io,report),joinpath(out,"report.json"),"w")
    report
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    TrainNativeRecurrent.run(ARGS[1],ARGS[2];H=length(ARGS)>2 ? parse(Int,ARGS[3]) : 64,
      train_blocks=length(ARGS)>3 ? parse(Int,ARGS[4]) : 8,
      eval_blocks=length(ARGS)>4 ? parse(Int,ARGS[5]) : 8)
end
