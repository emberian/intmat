#!/usr/bin/env julia
# Bounded whole-pipeline attribution, using the real native adapter and trainer.
include(joinpath(@__DIR__,"..","scripts","device_cellular_setup.jl"))
using KernelAbstractions,Statistics
include("train_native_recurrent.jl")
const TNR=TrainNativeRecurrent
const R=TNR.R
const RT=TNR.RT
const A=TNR.A
function stage(f)
    synchronize(backend);started=time_ns();result=f();synchronize(backend)
    (;seconds=(time_ns()-started)/1e9,result)
end
function objective_parity(theta,trainer,batch)
    c=trainer.cache;R.forward!(c,theta,batch.x,batch.h0,batch.memory_continuation)
    O=RT.O
    O.evaluate!(trainer.objective,c.y,batch.z,batch.noise,batch.reward,batch.return_continuation,batch.valid,batch.bootstrap,batch.prediction,batch.predvalid;
        count=batch.count,sigma=trainer.sigma,discount=trainer.discount,lower=trainer.lower,upper=trainer.upper,slope=trainer.slope,
        prediction_scale=trainer.prediction_scale,prediction_weight=trainer.prediction_weight)
    y=Float64.(Array(c.y));valid=Array(batch.valid);predvalid=Array(batch.predvalid)
    returns=zeros(Float64,size(valid));O.returns!(returns,Float64.(Array(batch.reward)),Float64.(Array(batch.return_continuation)),valid,Float64.(Array(batch.bootstrap));discount=trainer.discount)
    advantages=returns.-y[4,:,:]
    ref=O.loss_gradient(y,Float64.(Array(batch.z)),Float64.(Array(batch.noise)),returns,advantages,Float64.(Array(batch.prediction)),valid,predvalid;
        sigma=trainer.sigma,lower=Array(trainer.lower),upper=Array(trainer.upper),slope=Array(trainer.slope),
        prediction_scale=Array(trainer.prediction_scale),prediction_weight=trainer.prediction_weight)
    actual=Float64.(Array(trainer.objective.gradient));error=maximum(abs.(actual.-ref.gradient);init=0.)/max(1.,maximum(abs,ref.gradient;init=0.))
    error<=5e-5||Base.error("native objective derivative mismatch $error")
    (;scaled_error=error,valid_decisions=batch.count)
end
function run(out)
    rows=Any[]
    profile=joinpath(@__DIR__,"..","experiments","native","living-patch-parallel-support-64.toml")
    for (H,Tsteps) in ((64,8),(128,32))
        println((;stage="fixture",H,Tsteps));flush(stdout)
        env=A.make_env(profile;backend,seed=211,phase_time=6.)
        theta=R.to_backend(backend,RT.initialize(;H,seed=921));b=TNR.buffers(theta,size(env.x,2),Tsteps)
        trainer=RT.trainer(theta,size(env.x,2),Tsteps)
        for block in 1:2
            batch=TNR.rollout!(env,theta,b;T=Tsteps,seed=921,clock=(block-1)*Tsteps)
            RT.train_step!(theta,trainer,batch);TNR.recarry!(b,trainer,theta)
        end
        synchronize(backend);samples=Any[];checks=Any[]
        for block in 1:3
            # Separate synchronized phases attribute costs; outer uninterrupted
            # wall includes these diagnostic barriers and report construction.
            started=time_ns()
            rollout=stage(()->TNR.rollout!(env,theta,b;T=Tsteps,seed=921,clock=(block+1)*Tsteps))
            batch=rollout.result
            train=stage(()->RT.train_step!(theta,trainer,batch))
            carry=stage(()->TNR.recarry!(b,trainer,theta))
            elapsed=(time_ns()-started)/1e9
            push!(samples,(;rollout_seconds=rollout.seconds,train_seconds=train.seconds,recarry_seconds=carry.seconds,total_seconds=elapsed,
                valid_decisions=batch.count,physical_steps=round(Int,env.options.decision_time/env.options.dt)*Tsteps,
                simulated_time=Tsteps*env.options.decision_time,stats=train.result))
            println((;H,Tsteps,block,rollout_ms=1000rollout.seconds,train_ms=1000train.seconds,carry_ms=1000carry.seconds));flush(stdout)
        end
        push!(checks,objective_parity(theta,trainer,merge(b,(;count=count(Array(b.valid))))))
        observer=stage(()->A.summary(env))
        total=sum(s.total_seconds for s in samples);phys=sum(s.physical_steps for s in samples);decisions=sum(s.valid_decisions for s in samples)
        row=(;H,Tsteps,capacity=size(env.x,2),samples,objective_checks=checks,
            physical_steps_per_second=phys/total,valid_decisions_per_second=decisions/total,
            simulated_time_per_wall_second=sum(s.simulated_time for s in samples)/total,
            rollout_fraction=sum(s.rollout_seconds for s in samples)/total,observer_seconds=observer.seconds,physical=observer.result,
            source=Interior.DEVICE_CELLULAR_SOURCE_PIN,runtime_module_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"shared_recurrent_runtime.jl")))))
        push!(rows,row)
        open(io->JSON3.write(io,(;schema="native-recurrent-pipeline-profile-v1",runtime=device_cellular_runtime(),
            resources="CPU10–11 SMT siblings/one P-core;Julia2;8GiB/no swap;nice15;exclusive GPU",
            caveat="Single-world full native rollout+training+recarry. Diagnostic phase barriers included; constructor/JIT/two warmblocks and finalobserver excluded from rates. Short correctness/profile run, not evidence of learning benefit.",rows)),out,"w")
        println((;stage="done",H,Tsteps,physical_steps_per_second=row.physical_steps_per_second,rollout_fraction=row.rollout_fraction));flush(stdout)
        env=nothing;theta=nothing;b=nothing;trainer=nothing;GC.gc()
    end
end
run(ARGS[1])
