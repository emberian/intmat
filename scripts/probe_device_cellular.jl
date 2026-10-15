#!/usr/bin/env julia
length(ARGS)>0 && (ENV["DEVICE_BACKEND"]=ARGS[1])
include("device_cellular_setup.jl")
using Interior: DeviceCellular,device_cellular_step!,device_cellular_status!,device_cellular_statistics,device_cellular_snapshot
out=length(ARGS)>1 ? ARGS[2] : "outputs/device-cellular-probe"
ispath(out) && error("Use a fresh output directory")
function fixture()
    three=get(ENV,"DEVICE_PROBE_3D","1")=="1"
    w=World(Grid(8,8,three ? 8 : 1;extent=(8,8,three ? 8 : 1)));fill!(w.fields.nutrient,0)
    r=ResourcePool(w;potentials=(4,4,0),initial=(.2,.2,0),diffusion=(.02,.02,.01))
    x=three ? Float64[2 3.2 6;2 2 6;1 1.2 5] : Float64[2 3.2 6;2 2 6;0 0 0]
    c=CellularEcology(w,x,r;capacity=8,biomass=2,energy=1,mutation=0,developmental_partition=0)
    c.biomass[3]=.1
    get(ENV,"DEVICE_PROBE_RANDOM","1")=="0" && fill!(c.weights,0)
    h=size(c.hidden,1)
    c.weights[h+8,end,:].=2;c.weights[h+9,end,:].=2;c.weights[h+10,end,:].=2
    add_cell_bond!(c.bonds,1,2,1.8)
    w,c
end
wc,c=fixture();pc=DeviceCellular(wc,c)
wg=World(wc.grid,backend,NamedTuple(),NamedTuple(),0.,wc.seed);pg=DeviceCellular(wg,c)
initial=device_cellular_statistics(pg)
initial_audit=Interior.device_cellular_audit(device_cellular_snapshot(pg))
for _ in 1:parse(Int,get(ENV,"DEVICE_PROBE_STEPS","50"))
    device_cellular_step!(wc,pc,.01)
    device_cellular_step!(wg,pg,.01)
end
xc=device_cellular_snapshot(pc);xg=device_cellular_snapshot(pg)
errors=(position=maximum(abs.(xc.state.position.-xg.state.position)),
    metabolites=maximum(abs.(xc.state.metabolites.-xg.state.metabolites)),
    hidden=maximum(abs.(xc.state.hidden.-xg.state.hidden)),resources=maximum(abs.(xc.environment.fields.-xg.environment.fields)))
final=device_cellular_statistics(pg)
final_audit=Interior.device_cellular_audit(xg)
mkpath(out)
passed=all(x->x<=1e-4,values(errors)) && final_audit.finite
result=(passed=passed,tolerance=1e-4,runtime=device_cellular_runtime(),hidden_maximum=maximum(abs,xg.state.hidden),weight_maximum=maximum(abs,xg.state.weights),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,backend=kind,initial=initial,final=final,initial_audit=initial_audit,final_audit=final_audit,errors=errors,
    material_error=final_audit.material-initial_audit.material-final_audit.supplied+initial_audit.supplied+final_audit.removed-initial_audit.removed,
    potential_error=final_audit.potential-initial_audit.potential-final_audit.supplied_potential+initial_audit.supplied_potential+final_audit.removed_potential-initial_audit.removed_potential+final_audit.work-initial_audit.work+final_audit.heat-initial_audit.heat,
    model_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"..","src","device_cellular.jl")))),
    kernels_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"..","src","device_cellular_kernels.jl")))))
open(joinpath(out,"fixture.json"),"w") do io;JSON3.pretty(io,result);end
println(result)
passed || error("CPU/device primary array tolerance exceeded; inspect fixture.json")
