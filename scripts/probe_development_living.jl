#!/usr/bin/env julia
include("device_cellular_setup.jl")
include(joinpath(@__DIR__,"../studies/native_development_living.jl"))
length(ARGS)==2||error("PROFILE FRESH_OUTPUT")
ispath(ARGS[2])&&error("fresh output required")
W=NativeDevelopmentLiving;V=W.V
cfg=V.TOML.parsefile(ARGS[1]);cfg["cellular"]["biomass"]=1.79
fixture=tempname()*".toml";open(io->V.TOML.print(io,cfg),fixture,"w")
runs=[]
for be in (CPU(),backend)
    f=V.N.fixture(fixture;backend=be);w,p=f.w,f.p
    h=W.sidecar(p;rate=.1,horizon=4.,act_time=12.,inheritance_fraction=.5,
        (Symbol(k)=>v for(k,v) in cfg["temporal_policy"])...)
    copyto!(h.state.ids,p.cells.state.ids)
    weights=zeros(Float32,size(h.state.learned));weights[:,:,2].=.001f0;copyto!(h.state.learned,weights)
    origin=V.C.A.audit(p);frames=[]
    for duration in (.42,3.58,5.)
        V.L.advance!(w,p,h,duration;dt=.02)
        push!(frames,V.N.frame(w,p,h,origin))
    end
    push!(runs,(;frames,events=V.N.birth_events(h),temporal=map(Array,h.temporal),private=map(Array,h.state),capital=map(Array,h.capital),coefficient=h.coefficient))
end
cpu,gpu=runs;errors=Dict{String,Float64}()
for (a,b) in zip(cpu.frames,gpu.frames)
    [c.id for c in a.cells]==[c.id for c in b.cells]||error("resident identity differs")
    for field in (:position,:radius,:mass,:energy,:stores,:learning_work)
        errors[string(field)]=max(get(errors,string(field),0.),maximum(maximum(abs,getproperty(c,field).-getproperty(d,field)) for(c,d) in zip(a.cells,b.cells)))
    end
    for field in keys(a.material)
        errors[string(field)]=max(get(errors,string(field),0.),maximum(abs,getproperty(a.material,field).-getproperty(b.material,field)))
    end
end
for section in (:temporal,:private,:capital)
    a=getproperty(cpu,section);b=getproperty(gpu,section)
    for k in keys(a)
        x=getproperty(a,k);y=getproperty(b,k)
        if eltype(x)<:Integer
            x==y||error("integer state mismatch $section.$k")
        else
            errors["$section.$k"]=maximum(abs,x.-y;init=0.)
        end
    end
end
sum(gpu.private.actor_updates)>0||error("no paid actor updates")
sum(gpu.private.paid)>0||error("no work paid")
length(cpu.events)==length(gpu.events)>0||error("births absent/different")
any(e.copy_funded && e.copied_l1>0 && e.copy_work>0 for e in gpu.events)||error("nonzero paid inheritance absent")
maximum(values(errors))<1e-4||error("parity tolerance failed: $errors")
residual=maximum(max(abs(f.ledger.material_residual),abs(f.ledger.potential_residual)) for r in runs for f in r.frames)
residual<1e-8||error("physical ledger residual")
report=(;schema="developmental-living-native-gate-v1",coefficient=gpu.coefficient,errors,residual,actor_updates=sum(gpu.private.actor_updates),learning_work=sum(gpu.private.paid),
    cpu_births=cpu.events,gpu_births=gpu.events,recipe=cfg,runtime=device_cellular_runtime(),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
    scope="stoichiometric developmental-value objective; uninterrupted 9-time CPU/AMD parity; near-division fixture and nonzero initial allocation readout exercise paid birth; filter cadence .2, action cadence1, horizon4; no GPU fork claim")
open(io->JSON3.pretty(io,report),ARGS[2],"w");println((;maximum_error=maximum(values(errors)),residual,births=length(gpu.events)))
