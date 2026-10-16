# Invariants, zero-interface selectivity, time refinement and independent edge removals.
include(joinpath(@__DIR__,"../studies/brine_memory.jl"))
function brine_checks(out)
 ispath(out)&&error("Fresh check output required");mkpath(out)
 w=World(Grid(3,3,3));a=fill(2f0,3,3,3);b=fill(1f0,3,3,3);m=fill(.25f0,3,3,3);g=fill(.2f0,3,3,3)
 before=(copy(a),copy(b),copy(m));Interior.launch!(w,brine_react!,a,b,m,g,10f0,3f0,0f0,.35f0,.8333333f0)
 @assert a==before[1]&&b==before[2]&&m==before[3]
 fill!(g,1f0);mass0=sum(Float64,a)+sum(Float64,b)+sum(Float64,m)
 Interior.launch!(w,brine_react!,a,b,m,g,10f0,3f0,0f0,.35f0,.8333333f0)
 @assert minimum(a)>=0&&minimum(b)>=0&&minimum(m)>=0&&all(==(0),b)
 masserr=sum(Float64,a)+sum(Float64,b)+sum(Float64,m)-mass0;@assert masserr==0
 c=TOML.parsefile(length(ARGS)>1 ? ARGS[2] : joinpath(@__DIR__,"../experiments/art/brine-lace-pilot.toml"));reports=Dict{String,Any}()
 for arm in ("closed","no-mineral-feedback","phase-only","transport-only")
  reports[arm]=brine_arm(c,joinpath(out,arm),arm)
 end
 half=copy(c);half["dt"]/=2;half["steps"]*=2;reports["half-dt"]=brine_arm(half,joinpath(out,"half-dt"),"closed")
 finalfield(arm,name)=copy(reinterpret(Float32,read(joinpath(out,arm,last(reports[arm].frames).directory,name*".f32"))))
 phi=finalfield("closed","phi");mineral=finalfield("closed","mineral")
 ablations=Dict(arm=>(phase_max_difference=maximum(abs.(phi.-finalfield(arm,"phi"))),mineral_max_difference=maximum(abs.(mineral.-finalfield(arm,"mineral"))),mineral_mass=last(reports[arm].history).mineral_mass) for arm in ("no-mineral-feedback","phase-only","transport-only"))
 @assert ablations["phase-only"].mineral_max_difference>1e-3
 @assert ablations["transport-only"].phase_max_difference>1e-3
 phase_rms=sqrt(mean(abs2,phi.-finalfield("half-dt","phi")));mineral_rms=sqrt(mean(abs2,mineral.-finalfield("half-dt","mineral")))
 m0=last(reports["closed"].history).mineral_mass;mhalf=last(reports["half-dt"].history).mineral_mass
 relative_mineral_change=abs(mhalf-m0)/mhalf;@assert relative_mineral_change<.05&&phase_rms<.05
 result=(schema="brine-causal-and-timestep-checks",zero_interface_reaction_exactly_disabled=true,large_reaction_step_mass_error=masserr,physical_nucleation_gradient=c["nucleation_gradient"],physical_saturation_gradient=c["saturation_gradient"],dt=c["dt"],half_dt=half["dt"],time=reports["closed"].time,phase_half_dt_rms=phase_rms,mineral_half_dt_rms=mineral_rms,mineral_mass_half_dt_relative_change=relative_mineral_change,ablations=ablations,history=Dict(arm=>last(r.history) for (arm,r) in reports),interpretation="Independent removals of mineral→phase and mineral→transport edges, same initial finite inventories. Global signed-phase and A+B+mineral conservation are checked independently of local update identities. Same-grid half-dt comparison is a temporal sensitivity check, not spatial convergence or a thermodynamic calibration.")
 open(io->JSON3.pretty(io,result),joinpath(out,"checks.json"),"w");println(JSON3.write(result))
end
brine_checks(ARGS[1])
