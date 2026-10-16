#!/usr/bin/env julia
# Small routing/restore gate; throughput evidence remains the frozen prototype run.
include("living_gpu_profile.jl")
function production_gate(out)
    rows=Any[]
    for founders in (24,46)
        println((;stage="warm",founders));flush(stdout)
        cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","native","living-patch-parallel-support-64.toml"))
        cfg["cellular"]["founders"]=founders;cfg["native"]["support_execution"]="serial"
        path=tempname()*".toml";open(io->TOML.print(io,cfg),path,"w")
        f=N.fixture(path;backend);h=N.sidecar(f.p;dt=.02,horizon=2.,act_time=100.,inheritance_fraction=.5,launcher=PolicyLaunch(nothing))
        original=(;w=f.w,p=f.p,h);advance!(original,220)
        physical=Interior.device_habitat_checkpoint(original.w,original.p)
        candidate=merge(physical,(cells=merge(physical.cells,(support_execution=:parallel,)),))
        a=clone(original,physical);b=clone(original,candidate)
        advance!(a,8);advance!(b,8);checked=parity(a,b)
        sum(Array(b.h.state.actor_updates))>0||error("no paid updates")
        Interior.device_support_execution(b.p.cells)===:parallel||error("parallel routing lost")
        restored=clone(b,Interior.device_habitat_checkpoint(b.w,b.p))
        advance!(b,3);advance!(restored,3);continued=parity(b,restored)
        push!(rows,(;founders,exact_elements=checked,restore_elements=continued,births=Array(b.p.cells.environment.metadata)[3],
            paid_actor_updates=sum(Array(b.h.state.actor_updates)),execution=Interior.device_support_execution(restored.p.cells)))
        open(io->JSON3.write(io,(;schema="living-support-production-gate-v1",runtime=device_cellular_runtime(),
            source=Interior.DEVICE_CELLULAR_SOURCE_PIN,resources="CPU10–11 SMT siblings/one P-core;8GiB;no swap;nice15;exclusive GPU",rows)),out,"w")
        println(last(rows));flush(stdout)
    end
end
production_gate(ARGS[1])
