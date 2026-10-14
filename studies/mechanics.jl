using Interior, JSON3, Printf
isdefined(Interior,:BondNetwork) || Base.include(Interior,joinpath(@__DIR__,"..","src","mechanics.jl"))
out=length(ARGS)>0 ? abspath(ARGS[1]) : joinpath(@__DIR__,"..","outputs","mechanics")
mkpath(out);w=World(Grid(64,64,64);seed=77)
p=Interior.BondNetwork(w;size=(10,10,10),extent=(1.1,1.1,1.1),mass=.02,stiffness=8.,damping=12.,stretch_rate=.19,shear_rate=.035,failure_strain=.17,weak_seam=.34,jitter=.1)
frames=Any[];initial=Interior.mechanical_energy(p).total
function capture(index)
    name=@sprintf("frame-%04d",index)
    Interior.export_bond_graph(p,joinpath(out,name*".obj"))
    open(joinpath(out,name*"-graph.json"),"w") do io
        JSON3.pretty(io,Dict("array_order"=>"column-major", "position_shape"=>collect(size(p.position)),"bonds_shape"=>collect(size(p.bonds)),"position"=>p.position,"bonds"=>p.bonds,"node_ids"=>p.node_ids,"bond_ids"=>p.bond_ids,"damage"=>p.damage,"peak_strain"=>p.peak_strain,"threshold"=>p.threshold))
    end
    energy=Interior.mechanical_energy(p)
    push!(frames,Dict("name"=>name,"time"=>w.time,"broken"=>count(==(1),p.damage),"energy"=>energy,"energy_balance_error"=>energy.total-initial-energy.external_work))
    fields=Interior.mechanical_optical_fields(w,p;radius=.011)
    export_volume(fields,joinpath(out,name*"-volume"))
    render_volume(fields,joinpath(out,name*".ppm");settings=RenderSettings(width=256,height=256,steps=96,emission=3.,absorption=1.0,scattering=.15,exposure=1.3))
end
capture(0)
for frame in 1:10
    for i in 1:350;step!(w,p,.0005);end
    capture(frame)
    println("mechanical frame ",frame," broken ",count(==(1),p.damage)," / ",length(p.damage))
end
open(joinpath(out,"run.json"),"w") do io
    JSON3.pretty(io,Dict("title"=>"A yielding web remembers its tear","model"=>"disordered central-force graph, irreversible strain fracture; no continuum calibration or repair", "frames"=>frames,"seed"=>77,"dt"=>.0005,"nodes"=>length(p.node_ids),"bonds"=>length(p.bond_ids),"loading"=>Dict("stretch_speed"=>p.stretch_rate,"shear_speed"=>p.shear_rate),"fracture_reservoir"=>p.damage_energy,"damping_reservoir"=>p.damping_energy,"external_work"=>p.external_work))
end
println("Saved yielding specimen ",out)
