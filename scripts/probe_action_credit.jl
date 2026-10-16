#!/usr/bin/env julia
include("device_cellular_setup.jl")
using KernelAbstractions
include(joinpath(@__DIR__,"../studies/native_action_credit.jl"))
length(ARGS)==2||error("PROFILE FRESH_OUTPUT_JSON")
ispath(ARGS[2])&&error("fresh output required")
function residual(a,b)
    (;material=b.material-a.material-b.supplied+a.supplied+b.removed-a.removed,
      potential=b.potential-a.potential+b.work-a.work+b.heat-a.heat-b.supplied_potential+a.supplied_potential+b.removed_potential-a.removed_potential)
end
results=Any[]
for precision in ("Float64","Float32"),context in ((.1,.5),(.5,.1),(.2,.4),(.4,.2))
    f=NativeActionCredit.fixture(ARGS[1],context);cfg=f.config
    cfg["native"]["inventory_type"]=precision;cfg["native"]["energy_type"]=precision
    cw,cp,cf=NativeActionCredit.native_constructed_garden(CPU(),cfg)
    gw,gp,gf=NativeActionCredit.native_constructed_garden(backend,cfg)
    start=time()
    for _ in 1:200
        Interior.device_cellular_step!(cw,cp.cells,.02;anabolism=cp.habitat,check_status=false)
        Interior.device_cellular_step!(gw,gp.cells,.02;anabolism=gp.habitat,check_status=false)
    end
    Interior.device_cellular_status!(cp.cells);Interior.device_cellular_status!(gp.cells)
    cs=Interior.device_habitat_snapshot(cp);gs=Interior.device_habitat_snapshot(gp)
    errors=Dict(string(k)=>maximum(abs,Float64.(cs.cells.state[k]).-Float64.(gs.cells.state[k])) for k in (:mass,:metabolites,:energy,:hidden,:actions))
    cr=residual(cf.initial,Interior.device_habitat_audit(cs));gr=residual(gf.initial,Interior.device_habitat_audit(gs))
    tolerance=precision=="Float64" ? 1e-8 : 1e-4
    println((;precision,context,cr,gr,errors));flush(stdout)
    open(io->JSON3.pretty(io,(;results,precision,context,cr,gr,errors)),ARGS[2]*".progress","w")
    maximum(values(errors))<1e-4 || error("CPU/GPU mismatch")
    closure_pass=maximum(abs,values(cr))<tolerance && maximum(abs,values(gr))<tolerance
    precision=="Float64" && !closure_pass && error("wide closure failure")
    push!(results,(;precision,context,closure_pass,closure_tolerance=tolerance,steps=200,dt=.02,errors,cpu_residual=cr,gpu_residual=gr,elapsed_seconds=time()-start))
    println(last(results));flush(stdout)
end
# The new sidecar kernels must compile on the native GPU with identical inputs.
rng=MersenneTwister(41);h=8;rows=(4,5,7);shape=(64,36,3)
learned=randn(rng,Float32,shape).*.01f0;noise=randn(rng,Float32,64,3).*.15f0
inputs=randn(rng,Float32,64,36);trace=zeros(Float32,shape);adv=randn(rng,Float32,64)
marks=zeros(Float32,64,36,h+Interior.CELLULAR_ACTIONS)
dlearned=AMDGPU.ROCArray(learned);dnoise=AMDGPU.ROCArray(noise);dinputs=AMDGPU.ROCArray(inputs)
dtrace=AMDGPU.ROCArray(trace);dadv=AMDGPU.ROCArray(adv);dmarks=AMDGPU.ROCArray(marks)
for (be,m,l,n,x,t,a) in ((CPU(),marks,learned,noise,inputs,trace,adv),(backend,dmarks,dlearned,dnoise,dinputs,dtrace,dadv))
    NativeActionCredit.install!(be,64)(m,l,n,h,rows;ndrange=shape)
    NativeActionCredit.eligible!(be,64)(t,x,n,inv(.15f0^2);ndrange=shape)
    NativeActionCredit.reinforce!(be,64)(l,t,a,.003f0,.3f0;ndrange=shape)
    KernelAbstractions.synchronize(be)
end
sidecar=(marks=maximum(abs,marks.-Array(dmarks)),trace=maximum(abs,trace.-Array(dtrace)),learned=maximum(abs,learned.-Array(dlearned)))
maximum(values(sidecar))<1e-5||error("sidecar parity failure")
report=(;schema="native-action-credit-gpu-gate-v1",results,sidecar,runtime=device_cellular_runtime(),
    source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
    driver_sha256=bytes2hex(Interior.SHA.sha256(read(@__FILE__))),
    scope="eight bounded CPU/AMD paired coupled windows plus sidecar kernel parity; Float32 closure is explicitly reported pass/fail, only wide closure gates acceptance; no full resident-learning throughput claim")
open(io->JSON3.pretty(io,report),ARGS[2],"w")
