#!/usr/bin/env julia
include("device_cellular_setup.jl")
include(joinpath(@__DIR__,"../studies/native_cohort_construction_credit.jl"))
length(ARGS)==2||error("PROFILE FRESH_OUTPUT_JSON")
ispath(ARGS[2])&&error("fresh output required")
D=NativeCohortConstructionCredit;N=D.C
results=[]
for (horizon,dt,feed,duration) in ((2.,.02,true,6.),(2.,.01,true,6.),(2.,.02,false,6.),(8.,.02,true,18.))
    f=N.fixture(ARGS[1];feed,reserve_fraction=.2);runs=[]
    for be in (CPU(),backend)
        w=World(f.grid,be,NamedTuple(),NamedTuple(),0.,f.seed);p=Interior.device_habitat_restore(w,f.checkpoint)
        h=D.sidecar(p;dt,horizon);origin=N.A.audit(p);trace=[]
        for _ in 1:round(Int,duration/2)
            N.L.advance!(w,p,h,2.;dt);push!(trace,N.observation(w,p,h,origin))
        end
        cp=Interior.device_habitat_checkpoint(w,p)
        wr=World(f.grid,be,NamedTuple(),NamedTuple(),w.time,f.seed);pr=Interior.device_habitat_restore(wr,cp)
        N.A.audit(pr)==N.A.audit(p)||error("cohort physical checkpoint mismatch")
        push!(runs,(;trace,learned=Array(h.state.learned),critic=Array(h.state.critic)))
    end
    cpu,gpu=runs;errors=Dict{String,Float64}()
    for field in (:usable,:learning_work,:uptake_a,:uptake_b,:supplied,:work,:heat,:secretion_a,:secretion_b,:enzyme,:bound,:solubilized,:cumulative_requests,:funded_increment)
        errors[string(field)]=maximum(maximum(abs,getproperty(c,field).-getproperty(g,field)) for (c,g) in zip(cpu.trace,gpu.trace))
    end
    for field in (:cohorts_enqueued,:cohorts_matured,:cohorts_pending,:cohorts_invalidated,:calibration_skips,:actor_updates,:critic_updates,:clip_total)
        all(getproperty(c,field)==getproperty(g,field) for (c,g) in zip(cpu.trace,gpu.trace))||error("cohort count mismatch: $field")
    end
    errors["learned"]=maximum(abs,cpu.learned.-gpu.learned);errors["critic"]=maximum(abs,cpu.critic.-gpu.critic)
    maximum(values(errors))<1e-4||error("cohort CPU/AMD mismatch")
    all(t.all_residents_valid for r in runs for t in r.trace)||error("resident identity changed")
    residual=maximum(max(abs(t.material),abs(t.potential)) for r in runs for t in r.trace)
    push!(results,(;horizon,dt,feed,duration,errors,residual,cpu=cpu.trace,gpu=gpu.trace))
    println((;horizon,dt,feed,residual,maximum_error=maximum(values(errors))));flush(stdout)
end
report=(;schema="native-cohort-construction-gpu-v1",results,runtime=device_cellular_runtime(),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,
    habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,study_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"../studies/native_cohort_construction_credit.jl")))),
    scope="decision cohorts, actual H2 and H8 queue lengths, native paid actor/predictor updates, two physics timesteps, finite feed/no-feed, and H8 deterministic-drain transition")
open(io->JSON3.pretty(io,report),ARGS[2],"w")
