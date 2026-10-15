include("cellular_initialization.jl")

include("cellular_storage.jl")

"""Growing material cells with an inherited regulatory network and open chemical ports.

Positions are physical periodic coordinates in `[0, extent)`. Two imported
precursors must meet locally to build biomass; either can instead be oxidized.
Transporter, synthesis and oxidation capacities share one regulated enzyme
budget. Persistent contacts carry force, metabolites and developmental cues.
This is a deliberately synthetic physiology, not a calibrated biological model.

Chemical work is an explicit reservoir. `work_spent` includes maintenance,
synthesis and motor work; it is not all labelled heat. Mechanical potential has
its own coarse overdamped diagnostics, rather than a false thermodynamic closure
across changing cell radii. Division retains developmental state; only genotype
mutation introduces new hereditary information.
"""
mutable struct CellularEcology
    ids::Vector{Int}
    parents::Vector{Int}
    lineages::Vector{Int}
    generations::Vector{Int}
    positions::Matrix{Float64}
    biomass::Vector{Float64}
    metabolites::Matrix{Float64}
    energy::Vector{Float64}
    hidden::Matrix{Float32}
    adaptation::Matrix{Float32}
    weights::AbstractArray{Float32,3}
    plastic_marks::AbstractArray{Float32,3}
    last_actions::Matrix{Float32}
    age::Vector{Float64}
    resources::ResourcePool
    medium::Any
    bonds::CellBonds
    stencil::CellStencil
    seed::Int
    next_id::Int
    capacity::Int
    births::Int
    deaths::Int
    blocked_births::Int
    work_spent::Float64
    chemical_heat::Float64
    initial_energy::Float64
    consumed::Vector{Float64}
    junction_transferred::Float64
    synthesized::Float64
    oxidized::Float64
    shed::Float64
    inert_residual::Float64
    mechanical::Any
    parameters::NamedTuple
end

const CELLULAR_ACTIONS=11
const CELLULAR_SENSORS=20

function CellularEcology(w::World,positions::AbstractMatrix,resources::ResourcePool;
        medium=nothing,seed=w.seed,capacity=4096,hidden=8,biomass=1.,energy=1.,genome_storage=:dense,
        metabolites=.12,exchange_radius=1.,density=1.,store_capacity=.8,energy_capacity=3.,
        uptake_rate=1.2,half_saturation=.2,enzyme_exponent=1.4,
        growth_rate=.8,catabolic_rate=1.4,growth_work=.4,efficiency=.65,
        basal_rate=.035,neural_cost=.0002,starvation_rate=.16,
        minimum_mass=.25,division_mass=1.8,speed=.7,motor_cost=.5,
        junction_rate=.8,repulsion=8.,stiffness=3.,break_strain=.8,max_degree=6,
        response_time=.8,adaptation_time=6.,adaptation_strength=.3,
        plasticity_rate=.002,plasticity_decay=.02,mutation=.08,
        mutation_probability=.03,developmental_partition=.03,packing_barrier=.6)
    n=size(positions,2); h=Int(hidden)
    size(positions,1)==3 && all(isfinite,positions) || throw(ArgumentError("physical positions must be finite 3×N"))
    0<=n<=capacity && h>=2 || throw(ArgumentError("invalid cellular capacity/hidden size"))
    length(resources.fields)==3 && resources.potentials[3]==0 || throw(ArgumentError("cells require two energetic precursors and one zero-potential waste species"))
    pars=cellular_parameters(;exchange_radius,density,store_capacity,energy_capacity,uptake_rate,half_saturation,enzyme_exponent,growth_rate,catabolic_rate,growth_work,efficiency,basal_rate,neural_cost,starvation_rate,minimum_mass,division_mass,speed,motor_cost,junction_rate,repulsion,stiffness,break_strain,max_degree,response_time,adaptation_time,adaptation_strength,plasticity_rate,plasticity_decay,mutation,mutation_probability,developmental_partition,packing_barrier)
    cellular_validate_founder_inventory(pars,biomass,energy,metabolites;minimum_mass,store_capacity,energy_capacity)
    weights=cellular_founder_weights(seed,h,n)
    xyz=Float64.(positions)
    for d in 1:3; xyz[d,:].=d==3 && w.grid.size[3]==1 ? 0. : mod.(xyz[d,:],Float64(w.grid.extent[d])); end
    p=CellularEcology(collect(1:n),zeros(Int,n),collect(1:n),zeros(Int,n),xyz,
        fill(Float64(biomass),n),fill(Float64(metabolites),2,n),fill(Float64(energy),n),
        zeros(Float32,h,n),zeros(Float32,h,n),weights,zeros(Float32,size(weights)),
        zeros(Float32,CELLULAR_ACTIONS,n),zeros(n),resources,medium,CellBonds(),
        CellStencil(w.grid,xyz,exchange_radius),Int(seed),n+1,Int(capacity),0,0,0,
        0.,0.,n*Float64(energy),zeros(2),0.,0.,0.,0.,0.,nothing,pars)
    storage=Symbol(genome_storage)
    storage in (:dense,:indexed) || throw(ArgumentError("genome_storage must be dense or indexed"))
    storage===:indexed && _cellular_enable_indexed_storage!(p)
    validate(w,p)
    p
end

function validate(w::World,p::CellularEcology)
    n=length(p.ids); h=size(p.hidden,1)
    _cellular_validate_storage(p)
    w.backend isa CPU || throw(ArgumentError("CellularEcology requires CPU fields"))
    q=p.parameters
    all(x->isfinite(x)&&x>=0,values(q)) && all(>(0),(q.exchange_radius,q.density,q.store_capacity,q.energy_capacity,
        q.half_saturation,q.response_time,q.adaptation_time,q.motor_cost,q.minimum_mass,q.division_mass)) &&
        0<=q.efficiency<=1 && 0<=q.mutation_probability<=1 && q.minimum_mass<q.division_mass/2 &&
        q.max_degree isa Integer || throw(ArgumentError("invalid cellular parameters"))
    length(p.resources.fields)==length(p.resources.potentials)==3 && p.resources.potentials[3]==0 &&
        all(x->isfinite(x)&&x>=0,p.resources.potentials) && all(a->size(a)==w.grid.size,p.resources.fields) || throw(ArgumentError("invalid cellular resource ports"))
    n<=p.capacity && length(unique(p.ids))==n && all(<(p.next_id),p.ids) || throw(ArgumentError("invalid cellular IDs/capacity"))
    size(p.positions)==(3,n) && size(p.metabolites)==(2,n) && size(p.hidden)==size(p.adaptation)==(h,n) &&
        size(p.weights)==size(p.plastic_marks)==(h+CELLULAR_ACTIONS,2h+CELLULAR_SENSORS,n) && size(p.last_actions)==(CELLULAR_ACTIONS,n) || throw(DimensionMismatch("cellular arrays"))
    all(v->length(v)==n,(p.parents,p.lineages,p.generations,p.biomass,p.energy,p.age)) || throw(DimensionMismatch("cellular vectors"))
    all(a->all(isfinite,a),(p.positions,p.hidden,p.adaptation,p.weights,p.plastic_marks,p.last_actions)) &&
        all(a->all(x->isfinite(x)&&x>=0,a),(p.biomass,p.energy,p.metabolites,p.age)) || throw(ArgumentError("invalid cellular state"))
    all(>(0),p.biomass) || throw(ArgumentError("living cells require positive biomass"))
    p.resources.medium===nothing || throw(ArgumentError("cellular porosity requires ResourcePool.medium=nothing; pass the medium to CellularEcology instead"))
    size(p.resources.permeability)==w.grid.size && p.stencil.grid.size==w.grid.size && p.stencil.grid.extent==w.grid.extent || throw(DimensionMismatch("cellular grid"))
    nothing
end
stability_limit(w,p::CellularEcology)=min(.25,p.parameters.response_time/2)

cellular_radii(w,p) = cellular_radii(w,p,p.parameters.density)

function cellular_radii(w,p,density)
    volume=(p.biomass.+vec(sum(p.metabolites;dims=1)))./density
    w.grid.size[3]>1 ? cbrt.(volume.*(3/(4pi))) : sqrt.(volume./pi)
end

function _cellular_select!(p,indices)
    for name in (:ids,:parents,:lineages,:generations,:biomass,:energy,:age)
        setfield!(p,name,getfield(p,name)[indices])
    end
    for name in (:positions,:metabolites,:hidden,:adaptation,:last_actions)
        setfield!(p,name,getfield(p,name)[:,indices])
    end
    _cellular_select_storage!(p,indices)
    p
end
function reorder_ecology!(p::CellularEcology,indices)
    sort(collect(indices))==collect(eachindex(p.ids)) || throw(ArgumentError("invalid cellular permutation"))
    _cellular_select!(p,indices)
    rebuild_cell_stencil!(p.stencil,p.positions)
    p
end

_cellular_inputs(w,p,concentrations,gradients,radii) =
    _cellular_inputs(w,p,concentrations,gradients,radii,p.parameters,p.medium)

function _cellular_inputs(w,p,concentrations,gradients,radii,q,medium)
    n=length(p.ids); h=size(p.hidden,1); input=zeros(Float32,2h+CELLULAR_SENSORS,n)
    input[1:h,:].=p.hidden
    index=Dict(id=>i for (i,id) in enumerate(p.ids)); degree=zeros(Int,n); tension=zeros(n); pressure=zeros(n)
    for (i,j) in cellular_pairs(w.grid,p.positions,radii)
        distance=_cellular_length(_cellular_delta(w.grid,p.positions,i,j))
        overlap=max(0.,radii[i]+radii[j]-distance)
        # Dimensionless contact compression is observable even without adhesion.
        pressure[i]+=overlap/max(radii[i],eps()); pressure[j]+=overlap/max(radii[j],eps())
    end
    for (edge,rest) in zip(p.bonds.endpoints,p.bonds.rest_lengths)
        i=get(index,edge[1],0); j=get(index,edge[2],0); (i==0||j==0) && continue
        delta=_cellular_delta(w.grid,p.positions,i,j)
        strain=max(0.,sqrt(sum(abs2,delta))/max(rest,eps())-1)
        gate=((p.last_actions[9,i]+1)/2)*((p.last_actions[9,j]+1)/2)
        for row in 1:h
            input[h+row,i]+=gate*p.hidden[row,j]
            input[h+row,j]+=gate*p.hidden[row,i]
        end
        degree[i]+=1; degree[j]+=1; tension[i]+=strain; tension[j]+=strain
    end
    for i in 1:n
        degree[i]>0 && (input[h+1:2h,i]./=degree[i])
        base=2h
        for s in 1:2
            input[base+s,i]=_neural_squash(concentrations[s,i])
            for d in 1:3; input[base+2+3(s-1)+d,i]=_neural_squash(gradients[s,d,i]*q.exchange_radius); end
        end
        input[base+9,i]=Float32(p.metabolites[1,i]/max(p.biomass[i]*q.store_capacity,eps()))
        input[base+10,i]=Float32(p.metabolites[2,i]/max(p.biomass[i]*q.store_capacity,eps()))
        input[base+11,i]=_neural_squash(p.energy[i]/max(p.biomass[i],eps()))
        input[base+12,i]=_neural_squash(p.biomass[i]/q.division_mass)
        input[base+13,i]=Float32(degree[i]/max(q.max_degree,1))
        input[base+14,i]=_neural_squash(tension[i])
        input[base+15,i]=_neural_squash(pressure[i])
        if medium!==nothing
            cell=CartesianIndex(ntuple(d->dimsafe(p.positions[d,i]/w.grid.spacing[d],w.grid.size[d]),3))
            input[base+16,i]=_neural_squash(medium.phase.phi[cell])
            input[base+17,i]=_neural_squash(medium.phase.interface[cell])
            input[base+18,i]=_neural_squash(medium.signal[cell])
        end
        input[base+19,i]=_neural_squash(p.age[i]/20)
        input[base+20,i]=1
    end
    input
end

# Local additions return their unrepresentable remainder to the cell. Matter is
# never charged for a Float32 field increment that did not actually occur.
function _cellular_waste!(p,i,amount)
    deposit_cell_resource!(p.resources,3,p.stencil,i,amount)
end

# Degradation is paid by actual intracellular endpoint changes. Computing
# amount*donor/stock can overdraw a tiny donor (including a subnormal) by an ulp.
# Form a bounded fraction first, then account what the stored endpoint lost.
# Waste that its Float32 field cannot receive remains explicit inert inventory.
function _cellular_degrade_precursors!(p,i,amount)
    stock=sum(view(p.metabolites,:,i))
    isfinite(stock) && isfinite(amount) && 0<=amount<=stock || throw(ArgumentError("invalid precursor degradation amount"))
    all(x->isfinite(x)&&x>=0,view(p.metabolites,:,i)) || throw(ArgumentError("invalid precursor donor"))
    amount==0 && return (;mass=0.,potential=0.)
    fraction=amount/stock;mass=0.;potential=0.
    for s in 1:2
        before=p.metabolites[s,i]
        after=before-min(before,fraction*before)
        removed=before-after;p.metabolites[s,i]=after
        mass+=removed;potential+=removed*p.resources.potentials[s]
    end
    deposited=_cellular_waste!(p,i,mass)
    p.inert_residual+=mass-deposited
    (;mass,potential)
end

cellular_porosity!(w,p) = cellular_porosity!(w,p,p.parameters,p.medium)

function cellular_porosity!(w,p,q,medium)
    r=p.resources; fill!(r.permeability,0)
    s=p.stencil
    for i in eachindex(p.ids), edge in s.offsets[i]:s.offsets[i+1]-1
        r.permeability[s.sites[edge]]+=Float32((p.biomass[i]+sum(view(p.metabolites,:,i)))*s.weights[edge]/s.cell_volume)
    end
    for site in eachindex(r.permeability)
        background=medium===nothing ? 1f0 : medium.permeability[site]
        r.permeability[site]=Float32(background/(1+q.packing_barrier*r.permeability[site]))
    end
    nothing
end

# Keep the checkpoint-compatible parameter field abstract, but specialize the
# numerical loop on its concrete NamedTuple at this dispatch boundary.
ecology_step!(w::World,p::CellularEcology,dt::Real;anabolism=nothing) =
    _ecology_step!(w,p,dt,p.parameters,anabolism)
_ecology_step!(w::World,p::CellularEcology,dt::Real,q) = _ecology_step!(w,p,dt,q,nothing)

function _ecology_step!(w::World,p::CellularEcology,dt::Real,q,anabolism)
    check_dt(w,p,dt)
    # Resolve storage only after a possible logical-cell permutation.
    issorted(p.ids) || reorder_ecology!(p,sortperm(p.ids))
    _ecology_step_storage!(w,p,dt,q,p.weights,p.plastic_marks,anabolism)
end

function _ecology_step_storage!(w::World,p::CellularEcology,dt::Real,q,weights,marks,anabolism)
    delta=Float64(Float32(dt))
    n=length(p.ids); h=size(p.hidden,1)
    rebuild_cell_stencil!(p.stencil,p.positions)
    concentrations=sample_concentrations(p.resources,p.stencil)
    gradients=sample_concentration_gradients(p.resources,p.stencil)
    radii=cellular_radii(w,p,q.density)
    inputs=_cellular_inputs(w,p,concentrations,gradients,radii)
    outputs=zeros(Float32,h+CELLULAR_ACTIONS,n)
    neural_controller!(outputs,weights,marks,inputs)
    enzymes=zeros(Float64,4,n); requests=zeros(Float64,3,n)
    for i in 1:n
        # Smooth regulatory allocation, with a finite total catalytic capacity.
        for e in 1:4; enzymes[e,i]=exp(2Float64(outputs[h+3+e,i])); end
        enzymes[:,i]./=sum(view(enzymes,:,i))
        room=max(0.,q.store_capacity*p.biomass[i]-sum(view(p.metabolites,:,i)))
        for s in 1:2
            concentration=concentrations[s,i]
            requests[s,i]=delta*q.uptake_rate*p.biomass[i]*enzymes[s,i]^q.enzyme_exponent*concentration/(q.half_saturation+concentration)
        end
        total=requests[1,i]+requests[2,i]
        total>room && (requests[:,i].*=room/total)
        alpha=-expm1(-delta/q.response_time); beta=-expm1(-delta/q.adaptation_time)
        for j in 1:h
            # Apply adaptation inside the activation, not as unbounded state drift.
            target=tanh(atanh(clamp(outputs[j,i],-.999999f0,.999999f0))-q.adaptation_strength*p.adaptation[j,i])
            p.hidden[j,i]+=Float32(alpha*(target-p.hidden[j,i]))
            p.adaptation[j,i]+=Float32(beta*(p.hidden[j,i]-p.adaptation[j,i]))
        end
        p.last_actions[:,i].=view(outputs,h+1:h+CELLULAR_ACTIONS,i)
        # Bounded, local developmental plasticity; no immediate-fitness reward is
        # smuggled in. Its persistence and rates are part of the explicit model.
        @inbounds for col in axes(inputs,1),row in 1:h
            mark=marks[row,col,i]; activity=p.hidden[row,i]
            marks[row,col,i]=clamp(mark+Float32(delta*(q.plasticity_rate*(activity*inputs[col,i]-activity^2*mark)-q.plasticity_decay*mark)),-.3f0,.3f0)
        end
    end
    grants=allocate_cell_resources!(p.resources,p.stencil,requests)
    for i in 1:n,s in 1:2
        p.metabolites[s,i]+=grants[s,i]; p.consumed[s]+=grants[s,i]
    end
    openness=Float64.((p.last_actions[9,:].+1)./2)
    exchange=exchange_cell_metabolites!(p.metabolites,p.bonds,p.ids,p.biomass,delta;conductance=q.junction_rate,openness)
    p.junction_transferred+=exchange.transferred
    for i in 1:n
        # Catabolism frees work from precursor potential; all spent precursor
        # mass becomes the explicitly zero-potential environmental waste species.
        stock=sum(view(p.metabolites,:,i))
        potential=sum(p.metabolites[s,i]*p.resources.potentials[s] for s in 1:2)
        # Respiratory capacity is not a command to burn every incoming molecule.
        # Homeostasis draws toward an evolvable reserve target plus upcoming work.
        reserve=q.energy_capacity*p.biomass[i]*(p.last_actions[11,i]+1)/2
        planned=delta*((q.basal_rate+q.growth_work*q.growth_rate*enzymes[3,i])*p.biomass[i]+q.neural_cost*h+
            q.motor_cost*q.speed^2*sum(abs2,view(p.last_actions,1:3,i)))
        planned=_cellular_planned_work(planned,w,p,i,delta,q,anabolism)
        demand=max(0.,reserve+planned-p.energy[i])
        fuel_needed=stock>0 && q.efficiency*potential>0 ? demand*stock/(q.efficiency*potential) : 0.
        oxidize=min(stock,delta*q.catabolic_rate*p.biomass[i]*enzymes[4,i],fuel_needed)
        if oxidize>0
            degraded=_cellular_degrade_precursors!(p,i,oxidize)
            p.energy[i]+=q.efficiency*degraded.potential
            p.chemical_heat+=(1-q.efficiency)*degraded.potential;p.oxidized+=degraded.mass
        end
        maintenance=delta*(q.basal_rate*p.biomass[i]+q.neural_cost*h)
        paid=min(p.energy[i],maintenance); p.energy[i]-=paid; p.work_spent+=paid
        if paid<maintenance
            damaged=min(p.biomass[i],delta*q.starvation_rate*p.biomass[i]*(1-paid/maintenance))
            actual=_cellular_waste!(p,i,damaged); p.biomass[i]-=actual; p.shed+=actual
        end
        # Loss of body capacity degrades the excess stored precursor to inert
        # residue; its potential is dissipated, never credited as useful work.
        stock=sum(view(p.metabolites,:,i))
        excess=max(0.,stock-q.store_capacity*p.biomass[i])
        if excess>0
            degraded=_cellular_degrade_precursors!(p,i,excess)
            p.chemical_heat+=degraded.potential
        end
        # Equal amounts of complementary precursors become material structure.
        # Their remaining potential is dissipated; synthesis additionally uses work.
        _cellular_anabolic_build!(w,p,i,delta,
            delta*q.growth_rate*p.biomass[i]*enzymes[3,i],q,anabolism)
        overflow=max(0.,p.energy[i]-q.energy_capacity*p.biomass[i])
        p.energy[i]-=overflow; p.chemical_heat+=overflow
        p.age[i]+=delta
    end
    # A terminal body leaves before mechanics: complete starvation may produce
    # exactly zero radius, and dead cells cannot propel or transmit junctions.
    dead=findall(<(q.minimum_mass),p.biomass)
    for i in dead
        mass=p.biomass[i]+sum(view(p.metabolites,:,i))
        actual=_cellular_waste!(p,i,mass)
        p.inert_residual+=mass-actual
        p.chemical_heat+=sum(p.metabolites[s,i]*p.resources.potentials[s] for s in 1:2)+p.energy[i]
        p.energy[i]=0; p.deaths+=1
    end
    if !isempty(dead)
        _cellular_select!(p,setdiff(collect(eachindex(p.ids)),dead))
        living=Set(p.ids); keep=[all(id->id in living,edge) for edge in p.bonds.endpoints]
        p.bonds.endpoints=p.bonds.endpoints[keep]; p.bonds.rest_lengths=p.bonds.rest_lengths[keep]
    end
    n=length(p.ids)
    radii=cellular_radii(w,p,q.density)
    motility=Float64.(p.last_actions[1:3,:]).*q.speed
    w.grid.size[3]==1 && (motility[3,:].=0)
    before=sum(p.energy)
    p.mechanical=cellular_mechanics!(w.grid,p.ids,p.positions,radii,p.bonds,delta;
        adhesion=Float64.((p.last_actions[8,:].+1)./2),motility,
        repulsion=q.repulsion,stiffness=q.stiffness,break_strain=q.break_strain,max_degree=q.max_degree,
        energy=p.energy,motor_cost=q.motor_cost)
    p.mechanical=merge(p.mechanical,(ids=copy(p.ids),))
    p.work_spent+=before-sum(p.energy)
    # All lifecycle deposits use the *new* occupied support.
    rebuild_cell_stencil!(p.stencil,p.positions)
    parents=Int[]
    for i in 1:n
        if p.biomass[i]>=q.division_mass && p.last_actions[10,i]>0
            push!(parents,i)
        end
    end
    # The cap is a computational guard. When reached, unbiased keyed priorities
    # choose among simultaneously eligible parents; report every blocked event.
    slots=max(0,p.capacity-n)
    if length(parents)>slots
        sort!(parents;by=i->electrical_uniform(UInt64(p.seed)⊻UInt64(p.ids[i]),floor(Int,w.time/delta)+1))
        p.blocked_births+=length(parents)-slots; resize!(parents,slots)
    end
    sort!(parents;by=i->p.ids[i])
    _cellular_divide!(w,p,parents)
    rebuild_cell_stencil!(p.stencil,p.positions)
    cellular_porosity!(w,p)
    p
end

# The ordinary law remains an exact, concrete dispatch path. Optional processes
# may provide a separately explicit anabolic allocator without a global registry.
_cellular_planned_work(planned,w,p,i,delta,q,::Nothing)=planned

function _cellular_anabolic_build!(w,p,i,delta,requested_growth,q,::Nothing)
    build=min(2minimum(view(p.metabolites,:,i)),requested_growth)
    q.growth_work>0 && (build=min(build,p.energy[i]/q.growth_work))
    p.metabolites[:,i].-=build/2; p.biomass[i]+=build
    work=build*q.growth_work; p.energy[i]-=work; p.work_spent+=work
    p.chemical_heat+=build/2*sum(p.resources.potentials[1:2]); p.synthesized+=build
    nothing
end

_cellular_divide!(w,p,parents) = _cellular_divide!(w,p,parents,p.parameters)

function _cellular_divide!(w,p,parents,q)
    isempty(parents) && return p
    n=length(p.ids); births=length(parents); h=size(p.hidden,1)
    # One allocation per state array for the complete birth batch.
    for name in (:positions,:metabolites,:hidden,:adaptation,:last_actions)
        old=getfield(p,name); setfield!(p,name,_cellular_append_slots(old,births))
    end
    _cellular_grow_genomes!(p,births)
    # Re-seeding restores the complete stream, including Gaussian caches. Keep
    # one local generator per birth batch instead of allocating its large state
    # for every daughter; the daughter-specific seed and draw order are unchanged.
    rng=MersenneTwister(0)
    for (offset,i) in enumerate(parents)
        j=n+offset; id=p.next_id; p.next_id+=1
        Random.seed!(rng,string(p.seed,":division:",id,":",p.ids[i]))
        push!(p.ids,id); push!(p.parents,p.ids[i]); push!(p.lineages,p.lineages[i]); push!(p.generations,p.generations[i]+1)
        p.biomass[i]/=2; p.energy[i]/=2; p.metabolites[:,i]./=2
        push!(p.biomass,p.biomass[i]); push!(p.energy,p.energy[i]); _cellular_copy_slot!(p.metabolites,j,i)
        p.age[i]=0; push!(p.age,0.)
        _cellular_copy_slot!(p.weights,j,i)
        _cellular_mutate_genome!(rng,p.weights,j,q.mutation,q.mutation_probability)
        _cellular_copy_slot!(p.hidden,j,i); _cellular_copy_slot!(p.adaptation,j,i); _cellular_copy_slot!(p.plastic_marks,j,i)
        # Symmetric partition perturbation preserves mean developmental state.
        for k in 1:h
            shift=Float32(q.developmental_partition)*randn(rng,Float32)
            p.hidden[k,j]=clamp(p.hidden[k,j]+shift,-1f0,1f0); p.hidden[k,i]=clamp(p.hidden[k,i]-shift,-1f0,1f0)
        end
        _cellular_copy_slot!(p.last_actions,j,i)
        direction=randn(rng,3)
        w.grid.size[3]==1 && (direction[3]=0)
        norm(direction)==0 && (direction[1]=1)
        direction./=norm(direction)
        # Daughter placement is a small symmetry break; passive expansion then
        # creates space through physical contact rather than teleporting a cell.
        separation=.1q.exchange_radius
        for d in 1:3
            center=p.positions[d,i]
            p.positions[d,i]=d==3 && w.grid.size[3]==1 ? 0. : mod(center-separation*direction[d]/2,w.grid.extent[d])
            p.positions[d,j]=d==3 && w.grid.size[3]==1 ? 0. : mod(center+separation*direction[d]/2,w.grid.extent[d])
        end
        volume=(p.biomass[i]+sum(view(p.metabolites,:,i)))/q.density
        radius=w.grid.size[3]>1 ? cbrt(3volume/(4pi)) : sqrt(volume/pi)
        degree=count(edge->p.ids[i] in edge,p.bonds.endpoints)
        q.max_degree>degree && p.last_actions[8,i]>=0 && add_cell_bond!(p.bonds,p.ids[i],id,2radius)
        p.births+=1
    end
    p
end

advance!(w::World,p::CellularEcology,dt)=ecology_step!(w,p,dt)
function cellular_statistics(p::CellularEcology)
    n=length(p.ids)
    (population=n,births=p.births,deaths=p.deaths,blocked_births=p.blocked_births,
     max_generation=maximum(p.generations;init=0),surviving_lineages=length(unique(p.lineages)),
     biomass=sum(p.biomass),internal_material=sum(p.metabolites),stored_energy=sum(p.energy),
     work_spent=p.work_spent,chemical_heat=p.chemical_heat,
     internal_potential=sum(p.metabolites[s,i]*p.resources.potentials[s] for s in 1:2,i in 1:n;init=0.),
     bonds=length(p.bonds.endpoints),junction_transferred=p.junction_transferred,
     synthesized=p.synthesized,oxidized=p.oxidized,shed=p.shed,inert_residual=p.inert_residual,
     mean_age=n==0 ? 0. : sum(p.age)/n)
end
ecology_statistics(p::CellularEcology)=cellular_statistics(p)
function record_process_fields(p::CellularEcology)
    names=(:ids,:parents,:lineages,:generations,:positions,:biomass,:metabolites,:energy,:hidden,:adaptation,:weights,:plastic_marks,:last_actions,:age)
    merge(NamedTuple{names}(Tuple(copy(getfield(p,n)) for n in names)),
        (bond_endpoints=isempty(p.bonds.endpoints) ? zeros(Int,2,0) : hcat(collect.(p.bonds.endpoints)...),
         bond_rest_lengths=copy(p.bonds.rest_lengths),parameters=Tuple(values(p.parameters)),
         parameter_names=Tuple(string.(keys(p.parameters))),statistics=Tuple(values(cellular_statistics(p))),
         statistic_names=Tuple(string.(keys(cellular_statistics(p)))),
         position_coordinates="physical periodic [0,extent); centered rendering subtracts extent/2"))
end
