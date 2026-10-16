include("illuminated_memory.jl")
using LinearAlgebra
root=ARGS[1];config=ARGS[2];baseline=ARGS[3]
mkpath(root);c=TOML.parsefile(config)
# Zero absorption and zero power have distinct controls, both exact in this model.
d,_,_=initial_lenses(16)
q,l,_=absorbed_power(d;grid=3,ds=.004,alpha=0.,power=100.)
@assert iszero(sum(q)) && iszero(l.absorbed) && abs(l.emergent+l.unresolved-100)<1e-10
q,l,_=absorbed_power(d;grid=3,ds=.004,alpha=2.,power=0.)
@assert iszero(sum(q)) && iszero(l.absorbed) && iszero(l.emergent)
paths=Dict("coupled"=>baseline)
for mode in ("dark","frozen","half_refresh")
    cc=copy(c);cc["mode"]=mode=="half_refresh" ? "coupled" : mode
    mode=="half_refresh"&&(cc["refresh_steps"]=c["refresh_steps"]÷2)
    cc["id"]="illuminated-memory-"*mode
    cfg=joinpath(root,mode*".toml");open(io->TOML.print(io,cc),cfg,"w")
    out=joinpath(root,mode);illuminated(cfg,out);paths[mode]=out
end
reports=Dict(key=>JSON3.read(read(joinpath(path,"run.json"),String),Dict{String,Any}) for (key,path) in paths)
function finalfield(key,name="solid")
    r=reports[key];loadfield(paths[key],last(r["frames"]),name,r["shape"][1])
end
solid=Dict(key=>finalfield(key) for key in keys(paths));difference(a,b)=norm(Float64.(a).-Float64.(b))/max(norm(Float64.(b)),eps())
metrics=Dict("coupled_vs_dark_relative_solid"=>difference(solid["coupled"],solid["dark"]),"coupled_vs_frozen_relative_solid"=>difference(solid["coupled"],solid["frozen"]),"half_refresh_relative_solid"=>difference(solid["half_refresh"],solid["coupled"]))
for (key,r) in reports
    @assert abs(r["energy"]["residual"])<1e-7
    @assert abs(last(r["history"])["mass_error"])/r["initial_material_mass"]<1e-6
end
@assert metrics["coupled_vs_dark_relative_solid"]>1e-4
@assert metrics["coupled_vs_frozen_relative_solid"]>1e-5
# Independently halve ray integration step on final coupled state.
kwargs=(grid=c["ray_grid"],alpha=c["alpha"],power=c["power"],contrast=c["contrast"])
qa,la,ea=absorbed_power(solid["coupled"];kwargs...,ds=c["ray_ds"])
qb,lb,eb=absorbed_power(solid["coupled"];kwargs...,ds=c["ray_ds"]/2)
metrics["ray_halfstep_relative_dose"]=difference(qb,qa)
metrics["ray_halfstep_absorbed_power_relative"]=abs(lb.absorbed-la.absorbed)/la.absorbed
metrics["ray_halfstep_max_endpoint"]=maximum(norm(collect(a).-collect(b)) for (a,b) in zip(ea,eb))
r=reports["coupled"];e0=first(r["ledger"])["endpoints"];ef=last(r["ledger"])["endpoints"]
metrics["coupled_first_last_mean_endpoint_change"]=mean(norm(a.-b) for (a,b) in zip(e0,ef))
open(io->JSON3.pretty(io,(metrics=metrics,paths=paths,zero_absorption_pass=true,zero_power_pass=true)),joinpath(root,"checks.json"),"w")
println("CHECKS ",JSON3.write(metrics))
