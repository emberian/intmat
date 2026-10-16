#!/usr/bin/env julia
"""Same-physical-law initialization and short-step audit across spatial refinements."""
module LivingResolutionAudit
using Interior, KernelAbstractions, TOML
using Interior: JSON3
include("native_living_construction.jl")
const N=NativeLivingConstruction
function run(profile,out)
    ispath(out)&&error("fresh output required");mkpath(out);rows=Any[]
    cfg=TOML.parsefile(profile)
    for n in (16,32,64)
        c=deepcopy(cfg);c["grid"]["shape"]=[n,n,n]
        recipe=joinpath(out,"recipe-$n.toml");open(io->TOML.print(io,c),recipe,"w")
        f=N.fixture(recipe);w,p=f.w,f.p;h=N.sidecar(p;dt=c["native"]["dt"],horizon=2.,act_time=2.,inheritance_fraction=.5)
        before=Interior.device_habitat_snapshot(p);a=Interior.device_habitat_audit(before)
        Interior.device_anabolic_sense!(h,w,p.cells,h.options.dt)
        sensed=Array(p.cells.state.inputs)[1:24,32:34]
        start=time();N.L.advance!(w,p,h,.2;dt=c["native"]["dt"]);elapsed=time()-start
        after=Interior.device_habitat_snapshot(p);fr=N.frame(w,p,h,a)
        push!(rows,(;resolution=n,volume=p.cells.geometry.volume,stability_limit=Interior.stability_limit(w,p),dt=c["native"]["dt"],
            initial_material=a.material,initial_potential=a.potential,raw_totals=vec(sum(before.feedstock.fields;dims=1)),
            source_world_rates=collect(p.habitat.feedstock.parameters.source_rates).*collect(p.habitat.feedstock.parameters.source_totals),
            source_budget=vec(before.feedstock.resource_ledger[:,1]+before.feedstock.resource_ledger_low[:,1]),founder_positions=before.cells.state.position[1:24,:],
            material_receptors=sensed,support_storage_bytes=sizeof(p.cells.state.sites)+sizeof(p.cells.state.support),
            elapsed_seconds=elapsed,final_population=length(fr.cells),ledger=fr.ledger))
        println((;resolution=n,elapsed_seconds=elapsed));flush(stdout)
    end
    open(io->JSON3.write(io,(;schema="living-resolution-audit-v1",rows)),joinpath(out,"report.json"),"w")
end
end
if abspath(PROGRAM_FILE)==(@__FILE__);LivingResolutionAudit.run(ARGS[1],ARGS[2]);end
