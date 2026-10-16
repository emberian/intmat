#!/usr/bin/env julia
"""Replay the declared scene to inspect actual actuator funding and local access."""
module LivingConstructionAuthority
using Interior, KernelAbstractions, Statistics
using Interior: JSON3
include("native_living_construction.jl")
const M=NativeLivingConstruction
function run(profile,out;duration=32.)
    ispath(out)&&error("fresh report required")
    f=M.fixture(profile);w,p=f.w,f.p;h=M.sidecar(p;act_time=duration,inheritance_fraction=.5);origin=M.C.A.audit(p);trace=[]
    for _ in 1:round(Int,duration/2)
        M.L.advance!(w,p,h,2.)
        x=Interior.device_habitat_snapshot(p);s=x.cells.state;a=Interior.device_habitat_audit(x);active=findall(==(1),s.active);q=p.cells.parameters
        cl=Float64.(s.ledger).+Float64.(s.ledger_low);hl=Float64.(x.arrays.cell_ledger).+Float64.(x.arrays.cell_ledger_low)
        hidden=Int(p.cells.geometry.hidden)
        cells=[(;id=s.ids[i],generation=s.generations[i],mass=s.mass[i],energy=s.energy[i],stores=vec(s.metabolites[i,:]),
            reserve_target=Float64(q.energy_capacity)*s.mass[i]*(1+s.actions[i,11])/2,
            reserve_fraction=(1+s.actions[i,11])/2,growth_enzyme=s.enzymes[i,3],oxidation_enzyme=s.enzymes[i,4],
            raw_sensor=s.inputs[i,2hidden+17],catalyst_sensor=s.inputs[i,2hidden+18],
            growth_requested=hl[i,1],growth_allocated=hl[i,4],secretion_requested=hl[i,2]+hl[i,3],secretion_funded=hl[i,7]+hl[i,8],
            uptake=cl[i,3]+cl[i,4],oxidized=cl[i,7]) for i in active]
        push!(trace,(;time=w.time,cells,global_raw_concentration=sum(x.feedstock.fields[:,1:2])/prod(f.cfg["grid"]["extent"]),
            total_biomass=sum(s.mass[active]),solubilized=a.solubilized,enzyme=a.enzyme,bound=a.bound,
            material_residual=a.material-origin.material-a.supplied+origin.supplied+a.removed-origin.removed,
            potential_residual=a.potential-origin.potential+a.work-origin.work+a.heat-origin.heat-a.supplied_potential+origin.supplied_potential))
    end
    report=(;schema="living-construction-authority-v1",duration,source=Interior.DEVICE_CELLULAR_SOURCE_PIN,trace)
    open(io->JSON3.write(io,report),out,"w");report
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    LivingConstructionAuthority.run(ARGS[1],ARGS[2])
end
