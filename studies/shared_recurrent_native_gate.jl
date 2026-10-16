#!/usr/bin/env julia
# Runtime numerical/scaling gate only. Full simulation+learning throughput is a
# separate outer-driver measurement; these timings must not be called arena speed.
include(joinpath(@__DIR__,"..","scripts","device_cellular_setup.jl"))
using KernelAbstractions,Statistics
include("shared_recurrent_runtime.jl")
const R=SharedRecurrentRuntime
function relative_error(a,b)
    x=Array(a);y=Array(b)
    maximum(abs.(Float64.(x).-Float64.(y));init=0.)/max(1.,maximum(abs,Float64.(y);init=0.))
end
function iteration!(c,theta,x,h0,mask,dy,opt)
    R.forward!(c,theta,x,h0,mask);R.backward!(c,theta,x,h0,mask,dy)
    R.update!(theta,c.g,opt);R.check_status!(opt)
end
bytes(x::AbstractArray)=length(x)*sizeof(eltype(x))
bytes(x::NamedTuple)=sum(bytes,values(x);init=0)
bytes(x)=0
function run(out)
    rows=Any[];D=25;O=16
    for (H,B,Tsteps) in ((16,7,8),(64,512,8),(128,512,32))
        println((;stage="fixture",H,B,Tsteps));flush(stdout)
        rng=MersenneTwister(633+H);original=R.initialize(rng,D,H,O)
        x0=randn(rng,Float32,D,B,Tsteps);h00=.1f0*randn(rng,Float32,H,B)
        dy0=randn(rng,Float32,O,B,Tsteps)/(B*Tsteps);mask0=ones(Int32,B,Tsteps)
        mask0[1:3:end,1].=0;mask0[2:3:end,cld(Tsteps,2)].=0
        theta64=map(a->Float64.(a),original);c64=R.allocate(theta64,B,Tsteps)
        R.forward!(c64,theta64,Float64.(x0),Float64.(h00),mask0)
        R.backward!(c64,theta64,Float64.(x0),Float64.(h00),mask0,Float64.(dy0))
        theta=R.to_backend(backend,original);c=R.allocate(theta,B,Tsteps)
        transfer(a)=begin b=KernelAbstractions.zeros(backend,eltype(a),size(a)...);copyto!(b,a);b end
        x=transfer(x0);h0=transfer(h00);dy=transfer(dy0);mask=transfer(mask0)
        R.forward!(c,theta,x,h0,mask);R.backward!(c,theta,x,h0,mask,dy);KernelAbstractions.synchronize(backend)
        errors=Dict("h"=>relative_error(c.h,c64.h),"y"=>relative_error(c.y,c64.y),"dh0"=>relative_error(c.dh0,c64.dh0))
        for k in keys(theta);errors["g.$k"]=relative_error(c.g[k],c64.g[k]);end
        maximum(values(errors))<=5e-5||error("native BPTT mismatch $errors")
        initial=Array(c.dh0);all(iszero,initial[:,1:3:end])||error("identity reset gradient leaked")
        opt=R.optimizer(theta);opt64=R.optimizer(theta64)
        R.update!(theta,c.g,opt);R.check_status!(opt);R.update!(theta64,c64.g,opt64);R.check_status!(opt64)
        for k in keys(theta);errors["adam.$k"]=relative_error(theta[k],theta64[k]);end
        maximum(values(errors))<=5e-5||error("native Adam mismatch $errors")
        hn=similar(h0);yn=similar(c.y,Float32,O,B)
        R.step!(hn,yn,theta,view(x,:,:,1),h0,view(mask,:,1));KernelAbstractions.synchronize(backend)
        # Warm the exact component timing call path; no simulation is included.
        iteration!(c,theta,x,h0,mask,dy,opt)
        samples=Any[]
        for _ in 1:3
            GC.gc();KernelAbstractions.synchronize(backend)
            t=@timed iteration!(c,theta,x,h0,mask,dy,opt)
            t.compile_time==0||error("timed compilation $(t.compile_time)")
            push!(samples,(;seconds=t.time,gc_seconds=t.gctime,host_allocated_bytes=t.bytes))
        end
        allocated=bytes(theta)+bytes(c)+bytes(opt.m)+bytes(opt.v)+bytes(opt.norm)+bytes(opt.status)+bytes(x)+bytes(h0)+bytes(dy)+bytes(mask)+bytes(hn)+bytes(yn)
        row=(;H,B,Tsteps,D,O,errors,device_array_bytes=allocated,component_samples=samples,
            component_block_seconds=median(s.seconds for s in samples),valid_sample_steps=B*Tsteps,
            caveat="Components only: forward+BPTT+Adam+compactstatus. Excludes physics, rollout/objective, observers/transfers and resource contention. Not a full-pipeline throughput result.")
        push!(rows,row)
        open(io->JSON3.write(io,(;schema="shared-recurrent-native-gate-v1",runtime=device_cellular_runtime(),
            module_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"shared_recurrent_runtime.jl")))),
            resources="CPU10–11 SMT siblings, one physical P-core; Julia2;8GiB/no swap;nice15;exclusive GPU",rows)),out,"w")
        println((;H,B,Tsteps,max_error=maximum(values(errors)),component_ms=1000row.component_block_seconds,device_mib=allocated/2.0^20));flush(stdout)
    end
end
run(ARGS[1])
