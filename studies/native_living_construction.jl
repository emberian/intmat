#!/usr/bin/env julia
"""Moving/dividing native matter with paid private learning and explicit inheritance."""
module NativeLivingConstruction
using Interior, KernelAbstractions, TOML, SHA
using Interior: JSON3
include("native_cohort_construction_credit.jl")
const D=NativeCohortConstructionCredit
const J=D.J
const C=D.C
const L=D.L

@kernel function birth!(s,status,q,z,events,clock,fraction,adjust)
    child=@index(Global,Linear)
    if status[1]==0 && s.active[child]==1
        parent=s.new_parent[child]
        if parent>0
            id=s.ids[child]
            if id>length(events.child_id)
                KernelAbstractions.@atomic status[1]|=Int32(8)
            elseif events.child_id[id]==0
                parentvalid=s.active[parent]==1 && z.ids[parent]==s.ids[parent]
                gift=L.local_value(s,q,child) # immediately after native partition, before readout-copy work
                pending=Int32(0)
                if adjust && parentvalid
                    for k in 1:size(z.queue_ids,2)
                        if z.queue_ids[parent,k]==s.ids[parent]
                            z.queue_value[parent,k]-=gift;pending+=1
                        end
                    end
                end
                # No old occupant's return, predictor, action or private marks survive.
                z.ids[child]=id;z.ready[child]=0;z.baseline[child]=0.;z.last_reward[child]=0.;z.older_reward[child]=0.
                z.paid[child]=0.;z.birth_pending[child]=1
                for k in 1:size(z.queue_ids,2)
                    z.invalidated[child]+=z.queue_ids[child,k]!=0;z.queue_ids[child,k]=0
                end
                for c in 1:8;z.critic[child,c]=0f0;z.features[child,c]=0f0;z.current_features[child,c]=0f0;end
                amount=0.
                for a in 1:3,c in 1:size(z.learned,2)
                    value=parentvalid ? fraction*z.learned[parent,c,a] : 0f0
                    amount+=abs(Float64(value));z.learned[child,c,a]=0f0;z.trace[child,c,a]=0f0
                end
                cost=q.kappa*amount;paid=0.;funded=cost==0
                request=cost+2eps(s.energy[child])
                if cost>0 && s.energy[child]>=request
                    nextenergy,debit=Interior.dc_energy_take(s.energy[child],request)
                    if debit>=cost
                        s.energy[child]=nextenergy;Interior.dc_ledger_add!(s.ledger,s.ledger_low,child,1,debit)
                        z.paid[child]+=debit;paid=debit;funded=true
                    end
                end
                if funded && parentvalid
                    for a in 1:3,c in 1:size(z.learned,2);z.learned[child,c,a]=fraction*z.learned[parent,c,a];end
                end
                events.child_id[id]=id;events.parent_id[id]=s.ids[parent];events.time[id]=clock
                events.gift[id]=gift;events.copied_l1[id]=funded ? amount : 0.;events.copy_work[id]=paid
                events.funded[id]=funded;events.pending[id]=pending
            end
        end
    end
end
@kernel function newborn_observation!(s,status,q,z,h,memory)
    i=@index(Global,Linear)
    if status[1]==0 && s.active[i]==1 && z.birth_pending[i]==1 && z.ids[i]==s.ids[i]
        # Sensing has now run at the child's actual location. No invented trend.
        z.current_features[i,2]=C.feature(s,i,h,2);z.current_features[i,3]=C.feature(s,i,h,3)
        for a in 1:3
            mu=a==1 ? -.8f0 : a==2 ? 0f0 : -0.6931472f0
            for c in 1:size(z.learned,2);mu+=z.learned[i,c,a]*J.feature(s,q,i,h,c,memory);end
            z.drive[i,a]=mu # first off-cadence action deterministic, no unrecorded exploration
        end
        z.birth_pending[i]=0
    end
end
struct LivingHabitat{H,E}
    base::H
    events::E
    inheritance_fraction::Float32
    birth_transfer_adjusted::Bool
end
function Base.getproperty(h::LivingHabitat,k::Symbol)
    k in (:base,:events,:inheritance_fraction,:birth_transfer_adjusted) ? getfield(h,k) : getproperty(getfield(h,:base),k)
end
Base.setproperty!(h::LivingHabitat,k::Symbol,v)=setproperty!(getfield(h,:base),k,v)
function sidecar(p;inheritance_fraction=0.,birth_transfer_adjusted=true,event_capacity=4096,kwargs...)
    0<=inheritance_fraction<=1||error("inheritance fraction outside[0,1]")
    b=D.sidecar(p;kwargs...);n=Int(p.cells.geometry.capacity)
    alloc(T,d...)=fill!(similar(p.cells.state.inputs,T,d...),zero(T))
    z=merge(b.state,(birth_pending=alloc(Int32,n),))
    events=(child_id=alloc(Int32,event_capacity),parent_id=alloc(Int32,event_capacity),time=alloc(Float64,event_capacity),gift=alloc(Float64,event_capacity),
        copied_l1=alloc(Float64,event_capacity),copy_work=alloc(Float64,event_capacity),funded=alloc(Int32,event_capacity),pending=alloc(Int32,event_capacity))
    c=C.ConstructionHabitat(b.habitat,z,b.backend,b.options,b.original_masks,b.steps,b.launcher)
    LivingHabitat(D.CohortHabitat(J.JointHabitat(c,b.memory)),events,Float32(inheritance_fraction),birth_transfer_adjusted)
end
Interior.device_anabolic_argument(h::LivingHabitat)=Interior.device_anabolic_argument(h.base)
Interior.device_anabolic_begin!(h::LivingHabitat,w,p,dt)=Interior.device_anabolic_begin!(h.base,w,p,dt)
Interior.device_anabolic_prepare!(h::LivingHabitat,w,p,dt)=Interior.device_anabolic_prepare!(h.base,w,p,dt)
Interior.device_anabolic_finish!(h::LivingHabitat,w,p,dt)=Interior.device_anabolic_finish!(h.base,w,p,dt)
function Interior.device_anabolic_inherit!(h::LivingHabitat,w,p,step)
    Interior.device_anabolic_inherit!(h.base,w,p,step)
    s=p.state;q=p.parameters
    view=(active=s.active,ids=s.ids,new_parent=s.new_parent,energy=s.energy,metabolites=s.metabolites,ledger=s.ledger,ledger_low=s.ledger_low)
    chemistry=(efficiency=q.efficiency,potentials=q.potentials,kappa=h.options.kappa)
    L.policy_launch!(h,birth!,view,p.environment.status,chemistry,h.state,h.events,(h.steps+1)*Float64(h.options.dt),h.inheritance_fraction,h.birth_transfer_adjusted;ndrange=Int(p.geometry.capacity))
end
function Interior.device_anabolic_sense!(h::LivingHabitat,w,p,dt)
    dt==h.options.dt||error("living clock changed")
    Interior.device_anabolic_sense!(h.habitat,w,p,dt)
    s=p.state;q=p.parameters
    v=(active=s.active,ids=s.ids,energy=s.energy,metabolites=s.metabolites,mass=s.mass,inputs=s.inputs,hidden=s.hidden,adaptation=s.adaptation)
    chemistry=(efficiency=q.efficiency,potentials=q.potentials,energy_capacity=q.energy_capacity,store_capacity=q.store_capacity)
    L.policy_launch!(h,newborn_observation!,v,p.environment.status,chemistry,h.state,h.options.hidden,h.memory;ndrange=Int(p.geometry.capacity))
    D.choose_step!(h,p)
end
function fixture(profile;backend=CPU(),seed=nothing)
    cfg=TOML.parsefile(profile);seed!==nothing&&(cfg["seed"]=seed)
    w,p,capacity=C.A.native_constructed_garden(backend,cfg)
    (;w,p,cfg,capacity)
end
function frame(w,p,h,origin)
    x=Interior.device_habitat_snapshot(p);s=x.cells.state;a=Interior.device_habitat_audit(x);z=map(Array,h.state)
    material_residual=a.material-origin.material-a.supplied+origin.supplied+a.removed-origin.removed
    potential_residual=a.potential-origin.potential+a.work-origin.work+a.heat-origin.heat-a.supplied_potential+origin.supplied_potential+a.removed_potential-origin.removed_potential
    a.finite && abs(material_residual)<1e-8 && abs(potential_residual)<1e-8||error("living physical closure")
    cells=[(;id=s.ids[i],parent_id=s.parents[i],lineage=s.lineages[i],generation=s.generations[i],position=vec(s.position[i,:]),radius=s.radii[i],
        mass=s.mass[i],energy=s.energy[i],stores=vec(s.metabolites[i,:]),learned_l2=sqrt(sum(abs2,z.learned[i,:,:])),learning_work=z.paid[i]) for i in eachindex(s.active) if s.active[i]==1]
    material=(raw_a=x.feedstock.fields[:,1],raw_b=x.feedstock.fields[:,2],
        catalyst_a=reinterpret(Float32,x.arrays.enzyme_bits[:,1]),catalyst_b=reinterpret(Float32,x.arrays.enzyme_bits[:,2]),
        bound_a=x.arrays.bound[:,1],bound_b=x.arrays.bound[:,2],free_a=x.cells.environment.fields[:,1],free_b=x.cells.environment.fields[:,2])
    ledger=(material=a.material,potential=a.potential,work=a.work-origin.work,heat=a.heat-origin.heat,supplied=a.supplied-origin.supplied,
        removed=a.removed-origin.removed,supplied_potential=a.supplied_potential-origin.supplied_potential,removed_potential=a.removed_potential-origin.removed_potential,material_residual,potential_residual)
    (;time=w.time,cells,material,ledger,funded_secretion=a.secreted,cumulative_requests=vec(sum(z.requested;dims=1)),births=x.cells.environment.metadata[3],deaths=x.cells.environment.metadata[4])
end
function birth_events(h)
    e=map(Array,h.events)
    [(;time=e.time[i],parent_id=e.parent_id[i],child_id=e.child_id[i],child_usable_before_copy=e.gift[i],copied_l1=e.copied_l1[i],
        copy_work=e.copy_work[i],copy_funded=e.funded[i]==1,adjusted_pending_count=e.pending[i]) for i in eachindex(e.child_id) if e.child_id[i]!=0]
end
function run(profile,out;backend=CPU(),duration=8.,capture_time=1.,inheritance_fraction=0.,birth_transfer_adjusted=true,horizon=2.,rate=.003,noise_seed=6201)
    ispath(out)&&error("fresh output required");mkpath(out)
    write(joinpath(out,"executed-driver.jl"),read(@__FILE__))
    f=fixture(profile;backend);w,p=f.w,f.p;dt=f.cfg["native"]["dt"]
    h=sidecar(p;dt,horizon,act_time=duration,rate,learning=rate!=0,seed=noise_seed,inheritance_fraction,birth_transfer_adjusted)
    origin=C.A.audit(p);frames=Any[frame(w,p,h,origin)];started=time()
    for _ in 1:round(Int,duration/capture_time)
        L.advance!(w,p,h,capture_time;dt);push!(frames,frame(w,p,h,origin))
        println((;time=w.time,cells=length(last(frames).cells),births=last(frames).births,deaths=last(frames).deaths));flush(stdout)
    end
    metadata=(shape=f.cfg["grid"]["shape"],extent=f.cfg["grid"]["extent"],origin=[0.,0.,0.],dt,seed=f.cfg["seed"],source_rates=f.cfg["feedstock"]["source_rates"],
        source_budgets=f.cfg["feedstock"]["source_budgets"],horizon,actor_rate=rate,sigma=.3,inheritance_fraction,birth_transfer_adjusted,
        material_units="amount per voxel",field_order="x-fast",elapsed_seconds=time()-started,capacity=f.capacity,source=Interior.DEVICE_CELLULAR_SOURCE_PIN,
        habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,driver_sha256=bytes2hex(sha256(read(joinpath(out,"executed-driver.jl")))))
    report=(;schema="living-construction-v1",metadata,frames,birth_events=birth_events(h))
    open(io->JSON3.write(io,report),joinpath(out,"report.json"),"w");report
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    NativeLivingConstruction.run(ARGS[1],ARGS[2];duration=length(ARGS)>2 ? parse(Float64,ARGS[3]) : 8.)
end
