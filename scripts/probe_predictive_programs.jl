#!/usr/bin/env julia
include("device_cellular_setup.jl")
include("../studies/native_predictive_programs.jl")
const P=NativePredictivePrograms
length(ARGS)==2||error("PROFILE OUTPUT")
ispath(ARGS[2])&&error("fresh output required")
cfg=P.TOML.parsefile(ARGS[1]);cfg["cellular"]["biomass"]=1.799
fixture=tempname()*".toml";open(io->P.TOML.print(io,cfg),fixture,"w")
runs=[]
for be in (CPU(),backend)
    f=P.N.fixture(fixture;backend=be);w,p=f.w,f.p;h=P.sidecar(p;forced=1)
    origin=P.N.C.A.audit(p)
    P.L.advance!(w,p,h,2.;dt=.02)
    frame=P.N.frame(w,p,h,origin)
    push!(runs,(;frame,model=map(Array,h.model),events=P.N.birth_events(h),observation=P.observation(w,p,h)))
end
cpu,gpu=runs;errors=Dict{String,Float64}()
for k in keys(cpu.model)
    x=getproperty(cpu.model,k);y=getproperty(gpu.model,k)
    if eltype(x)<:Integer
        x==y||error("integer model state differs: $k")
    else
        errors[string(k)]=maximum(abs,x.-y;init=0.)
    end
end
[c.id for c in cpu.frame.cells]==[c.id for c in gpu.frame.cells]||error("identity mismatch")
for k in (:mass,:energy,:stores,:position)
    errors[string(k)]=maximum(maximum(abs,getproperty(a,k).-getproperty(b,k)) for(a,b) in zip(cpu.frame.cells,gpu.frame.cells))
end
for k in keys(cpu.frame.material)
    errors["material_"*string(k)]=maximum(abs,getproperty(cpu.frame.material,k).-getproperty(gpu.frame.material,k))
end
sum(gpu.model.work)>0&&sum(gpu.model.updates)>0||error("no paid predictor updates")
sum(gpu.model.funded)>0||error("no funded construction")
length(gpu.events)>0||error("no actual births in native gate")
maximum(values(errors))<1e-4||error("native mismatch $errors")
residual=maximum(max(abs(r.frame.ledger.material_residual),abs(r.frame.ledger.potential_residual)) for r in runs)
residual<1e-8||error("ledger closure")
report=(;schema="native-predictive-programs-gate-v1",errors,residual,births=length(gpu.events),work=sum(gpu.model.work),updates=sum(gpu.model.updates),funded=sum(gpu.model.funded),
    runtime=device_cellular_runtime(),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
    scope="two-time CPU/AMD uninterrupted parity with actual births, paid prediction, build dwell and harvest tail; not a learning-benefit comparison")
open(io->JSON3.pretty(io,report),ARGS[2],"w");println(report)
