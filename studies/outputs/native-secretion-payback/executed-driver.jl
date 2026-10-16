#!/usr/bin/env julia
"""Paid secretion pulse versus no-secretion on the shared constructive schedule."""
module NativeSecretionPayback
using Interior, KernelAbstractions, Statistics, SHA
using Interior: JSON3
include("native_action_credit.jl")
@kernel function pulse!(a,s,status,rate,exponent,dt,investment,producer)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1
        value=s.ids[i]==producer ? investment : 0f0
        # This overrides only the investment actuator. Production substrate,
        # energy affordability, conservative deposition and ledgers remain active.
        a.expression[i,1]=value;a.expression[i,2]=.5f0;a.expression[i,3]=.5f0
        request=dt*rate*s.mass[i]*value*.5f0^exponent
        a.requests[i,1]=request;a.requests[i,2]=request
    end
end
mutable struct PulseHabitat{H,B}
    habitat::H
    backend::B
    dt::Float32
    duration_steps::Int
    investment::Float32
    producer::Int32
    steps::Int
end
Interior.device_anabolic_argument(h::PulseHabitat)=Interior.device_anabolic_argument(h.habitat)
Interior.device_anabolic_begin!(h::PulseHabitat,w,p,dt)=Interior.device_anabolic_begin!(h.habitat,w,p,dt)
Interior.device_anabolic_sense!(h::PulseHabitat,w,p,dt)=Interior.device_anabolic_sense!(h.habitat,w,p,dt)
Interior.device_anabolic_inherit!(h::PulseHabitat,w,p,step)=Interior.device_anabolic_inherit!(h.habitat,w,p,step)
function Interior.device_anabolic_prepare!(h::PulseHabitat,w,p,dt)
    dt==h.dt||error("pulse physics clock changed")
    Interior.device_anabolic_prepare!(h.habitat,w,p,dt)
    a=h.habitat.arrays;s=p.state;q=h.habitat.parameters
    pulse!(h.backend,64)((expression=a.expression,requests=a.requests),(active=s.active,ids=s.ids,mass=s.mass),
        p.environment.status,q.secretion_rate,q.allocation_exponent,dt,
        h.steps<h.duration_steps ? h.investment : 0f0,h.producer;ndrange=Int(p.geometry.capacity))
end
function Interior.device_anabolic_finish!(h::PulseHabitat,w,p,dt)
    Interior.device_anabolic_finish!(h.habitat,w,p,dt);h.steps+=1
end
function fixture(profile;seed=211,raw=2.,free=.03)
    cfg=NativeActionCredit.fixture(profile,(free,free);seed).config
    cfg["cellular"]["speed"]=0.;cfg["feedstock"]["initial"]=[raw,raw]
    cfg["feedstock"]["niche_contrast"]=0.
    w,p,_=NativeActionCredit.native_constructed_garden(CPU(),cfg)
    (;grid=w.grid,seed=w.seed,checkpoint=Interior.device_habitat_checkpoint(w,p),config=cfg)
end
function observe(p,origin,time)
    snap=Interior.device_habitat_snapshot(p);a=Interior.device_habitat_audit(snap)
    s=snap.cells.state;h=snap.arrays
    material=a.material-origin.material-a.supplied+origin.supplied+a.removed-origin.removed
    potential=a.potential-origin.potential+a.work-origin.work+a.heat-origin.heat-a.supplied_potential+origin.supplied_potential+a.removed_potential-origin.removed_potential
    a.finite && abs(material)<1e-8 && abs(potential)<1e-8 || error("secretion pulse accounting failed")
    all(s.active[1:4].==1) && s.ids[1:4]==Int32[1,2,3,4]||error("fixed-resident fixture changed identity")
    ledger=Float64.(s.ledger).+Float64.(s.ledger_low)
    hledger=Float64.(h.cell_ledger).+Float64.(h.cell_ledger_low)
    q=p.cells.parameters
    usable=[Float64(s.energy[i])+Float64(q.efficiency)*sum(Float64(s.metabolites[i,k])*q.potentials[k] for k in 1:2) for i in 1:4]
    (;time,usable,uptake_a=ledger[1:4,3],uptake_b=ledger[1:4,4],
      work=ledger[1:4,1].+hledger[1:4,9],heat=ledger[1:4,2].+hledger[1:4,10],
      secreted_a=hledger[1:4,7],secreted_b=hledger[1:4,8],enzyme=a.enzyme,bound=a.bound,solubilized=a.solubilized,
      raw=[sum(Float64,snap.feedstock.fields[:,k]) for k in 1:2],
      free=[sum(Float64,snap.cells.environment.fields[:,k]) for k in 1:2],
      permeability_min=minimum(h.permeability),material,potential)
end
function material_frame(p,time)
    x=Interior.device_habitat_snapshot(p);a=x.arrays;v=Float64(p.cells.geometry.volume)
    (;time,shape=Int.(p.cells.geometry.dims),extent=Float64.(p.cells.geometry.extent),
      positions=x.cells.state.position[1:4,:],radii=x.cells.state.radii[1:4],
      raw=vec(sum(Float64.(x.feedstock.fields[:,1:2]);dims=2))./v,
      catalyst=vec(sum(Float64.(reinterpret(Float32,a.enzyme_bits));dims=2))./v,
      bound=vec(sum(Float64.(a.bound);dims=2))./v,
      free=vec(sum(Float64.(x.cells.environment.fields[:,1:2]);dims=2))./v)
end
function branch(f;backend=CPU(),investment=.8,duration=.5,horizon=16.,dt=.02,capture=.1,frames=false)
    for t in (duration,horizon,capture);isapprox(t/dt,round(t/dt);atol=1e-8)||error("clocks must divide dt");end
    w=World(f.grid,backend,NamedTuple(),NamedTuple(),0.,f.seed);p=Interior.device_habitat_restore(w,f.checkpoint)
    h=PulseHabitat(p.habitat,backend,Float32(dt),round(Int,duration/dt),Float32(investment),Int32(1),0)
    origin=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p));trace=Any[observe(p,origin,0.)]
    fields=frames ? Any[material_frame(p,0.)] : Any[]
    frame_steps=[round(Int,t/dt) for t in (.2,.4,1.,4.,horizon)]
    for step in 1:round(Int,horizon/dt)
        Interior.device_cellular_step!(w,p.cells,dt;anabolism=h,check_status=false)
        if step%round(Int,capture/dt)==0;push!(trace,observe(p,origin,step*dt));end
        frames && step in frame_steps && push!(fields,material_frame(p,step*dt))
    end
    w.time=horizon
    (;trace,fields,final=Interior.device_habitat_snapshot(p))
end
function compare(control,pulse;duration=.5)
    differences=map(control.trace,pulse.trace) do c,p
        delta=p.usable.-c.usable
        (;time=c.time,private=delta[1],neighbors=sum(delta[2:4]),total=sum(delta),
          uptake_a=p.uptake_a.-c.uptake_a,uptake_b=p.uptake_b.-c.uptake_b,
          work=p.work.-c.work,enzyme=p.enzyme.-c.enzyme,bound=p.bound.-c.bound,solubilized=p.solubilized.-c.solubilized)
    end
    function payback(key)
        indices=findall(d->d.time>0 && getproperty(d,key)>1e-8,differences)
        first=isempty(indices) ? nothing : differences[indices[1]].time
        stable=findfirst(i->differences[i].time>0 && all(j->getproperty(differences[j],key)>1e-8,i:length(differences)),eachindex(differences))
        post=findfirst(d->d.time>=duration && getproperty(d,key)>1e-8,differences)
        (;first,after_command=post===nothing ? nothing : differences[post].time,sustained_until_end=stable===nothing ? nothing : differences[stable].time)
    end
    (;differences,private_payback=payback(:private),group_payback=payback(:total),neighbor_benefit=payback(:neighbors))
end
function run(profile,out;seed=211,horizon=16.,duration=.5)
    ispath(out)&&error("fresh output required");mkpath(out)
    write(joinpath(out,"executed-driver.jl"),read(@__FILE__));results=Any[]
    for raw in (0.,2.)
        f=fixture(profile;seed,raw)
        frames=seed==211 && raw==2.
        control=branch(f;investment=0.,duration,horizon,frames);pulse=branch(f;duration,horizon,frames)
        result=(;raw,config=f.config,control=control.trace,pulse=pulse.trace,compare(control,pulse;duration)...,
          source=f.checkpoint.cells.source,habitat_source=f.checkpoint.source)
        if frames
            open(io->JSON3.write(io,(;control=control.fields,pulse=pulse.fields)),joinpath(out,"material-frames.json"),"w")
        end
        push!(results,result)
        println((;raw,private_payback=result.private_payback,group_payback=result.group_payback,final=last(result.differences)));flush(stdout)
    end
    report=(;schema="native-secretion-payback-v2",seed,horizon,duration,investment=.8,allocation=.5,dt=.02,
      driver_sha256=bytes2hex(sha256(read(joinpath(out,"executed-driver.jl")))),results,
      scope="programmed single-resident paid secretion pulse versus matched no-secretion action; closed conserved chemistry; immobile nondividing four-founder material; no learning or cooperation claim")
    open(io->JSON3.pretty(io,report),joinpath(out,"report.json"),"w");report
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    length(ARGS)==2||error("PROFILE FRESH_OUTPUT")
    NativeSecretionPayback.run(ARGS...)
end
