#!/usr/bin/env julia
include("device_cellular_setup.jl")
include(joinpath(@__DIR__,"../studies/native_living_construction.jl"))
length(ARGS)==2||error("PROFILE FRESH_OUTPUT_JSON")
ispath(ARGS[2])&&error("fresh output required")
M=NativeLivingConstruction;runs=[]
for be in (CPU(),backend)
    f=M.fixture(ARGS[1];backend=be);w,p=f.w,f.p
    h=M.sidecar(p;inheritance_fraction=.5,act_time=8.)
    # Deliberately nonzero initial readout for a copy-payment parity fixture,
    # not a trained-policy or inherited-benefit claim.
    copyto!(h.state.ids,p.cells.state.ids);fill!(h.state.learned,.001f0)
    origin=M.C.A.audit(p);frames=[]
    for _ in 1:4
        M.L.advance!(w,p,h,1.);push!(frames,M.frame(w,p,h,origin))
    end
    push!(runs,(;frames,events=M.birth_events(h),learned=Array(h.state.learned),paid=Array(h.state.paid)))
end
cpu,gpu=runs;errors=Dict{String,Float64}()
for (a,b) in zip(cpu.frames,gpu.frames)
    [c.id for c in a.cells]==[c.id for c in b.cells]||error("native living IDs differ")
    for field in (:position,:radius,:mass,:energy,:stores,:learning_work)
        err=maximum(maximum(abs,getproperty(c,field).-getproperty(d,field)) for (c,d) in zip(a.cells,b.cells))
        errors[string(field)]=max(get(errors,string(field),0.),err)
    end
    for field in keys(a.material)
        errors[string(field)]=max(get(errors,string(field),0.),maximum(abs,getproperty(a.material,field).-getproperty(b.material,field)))
    end
end
length(cpu.events)==length(gpu.events)>0||error("native living birth count differs or absent")
all(e.copied_l1>0 && e.copy_work>0 && e.copy_funded for e in gpu.events)||error("copy payment fixture not exercised")
for (a,b) in zip(cpu.events,gpu.events)
    a.parent_id==b.parent_id && a.child_id==b.child_id && a.adjusted_pending_count==b.adjusted_pending_count||error("birth identity/credit differs")
end
errors["learned"]=maximum(abs,cpu.learned.-gpu.learned);errors["paid"]=maximum(abs,cpu.paid.-gpu.paid)
maximum(values(errors))<1e-4||error("native living CPU/AMD mismatch")
residual=maximum(max(abs(f.ledger.material_residual),abs(f.ledger.potential_residual)) for r in runs for f in r.frames)
report=(;schema="native-living-construction-gpu-v1",errors,residual,cpu_births=cpu.events,gpu_births=gpu.events,runtime=device_cellular_runtime(),
    source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
    scope="16³ moving/dividing native world, four-time coupled parity, deliberately nonzero initial readout to exercise paid child copy; no trained-policy claim")
open(io->JSON3.pretty(io,report),ARGS[2],"w")
println((;maximum_error=maximum(values(errors)),residual,births=length(gpu.events)))
