#!/usr/bin/env julia
include("device_cellular_setup.jl")
length(ARGS) in (2,3)||error("PROFILE FRESH_OUTPUT_JSON [joint]")
joint=length(ARGS)==3 && ARGS[3]=="joint"
include(joinpath(@__DIR__,joint ? "../studies/native_joint_construction_credit.jl" : "../studies/native_construction_credit.jl"))
ispath(ARGS[2])&&error("fresh output required")
N=joint ? NativeJointConstructionCredit.C : NativeConstructionCredit
make_sidecar=joint ? NativeJointConstructionCredit.sidecar : N.sidecar
results=[]
for feed in (false,true),dt in (.02,.01)
    f=N.fixture(ARGS[1];feed,reserve_fraction=.2);runs=[]
    for be in (CPU(),backend)
        w=World(f.grid,be,NamedTuple(),NamedTuple(),0.,f.seed);p=Interior.device_habitat_restore(w,f.checkpoint)
        h=make_sidecar(p;dt);origin=N.A.audit(p);trace=[]
        for _ in 1:3
            N.L.advance!(w,p,h,2.;dt);push!(trace,N.observation(w,p,h,origin))
        end
        cp=Interior.device_habitat_checkpoint(w,p)
        wr=World(f.grid,be,NamedTuple(),NamedTuple(),0.,f.seed);pr=Interior.device_habitat_restore(wr,cp)
        N.A.audit(pr)==N.A.audit(p)||error("phase checkpoint accounting mismatch")
        push!(runs,(;trace,learned=Array(h.state.learned),critic=Array(h.state.critic)))
    end
    cpu,gpu=runs;errors=Dict{String,Float64}()
    for field in (:usable,:learning_work,:uptake_a,:uptake_b,:supplied,:supplied_potential,:work,:heat,:secretion_a,:secretion_b,:enzyme,:bound,:solubilized,:cumulative_requests)
        errors[string(field)]=maximum(maximum(abs,getproperty(c,field).-getproperty(g,field)) for (c,g) in zip(cpu.trace,gpu.trace))
    end
    errors["learned"]=maximum(abs,cpu.learned.-gpu.learned);errors["critic"]=maximum(abs,cpu.critic.-gpu.critic)
    maximum(values(errors))<1e-4||error("construction CPU/AMD mismatch")
    all(t.all_residents_valid for r in runs for t in r.trace)||error("resident identity changed")
    residual=maximum(max(abs(t.material),abs(t.potential)) for r in runs for t in r.trace)
    push!(results,(;feed,dt,errors,residual,cpu=cpu.trace,gpu=gpu.trace))
    println((;feed,dt,residual,maximum_error=maximum(values(errors))));flush(stdout)
end
report=(;schema="native-construction-credit-gpu-v1",joint,results,runtime=device_cellular_runtime(),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,
    habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,study_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,joint ? "../studies/native_joint_construction_credit.jl" : "../studies/native_construction_credit.jl")))),
    scope="paid habitat policy, reserve fraction .2 founder, optional joint reserve action with recurrent features, finite raw sources on/off, two timesteps, swapped-mask checkpoint restore, six-time native parity")
open(io->JSON3.pretty(io,report),ARGS[2],"w")
