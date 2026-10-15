"""Cell-built extracellular catalysts turn inaccessible feedstock into shared food.

The two enzyme materials are immobile, persistent and inert in the chemical-work
ledger. Cells construct them from equal masses of both intracellular precursors
and pay work. Each enzyme catalyses release of its matching feedstock; it does
not create matter or chemical potential. Enzyme decay returns inert waste.

A separate inherited readout of cellular developmental state and local chemistry
regulates investment and allocation. This adds an evolvable effector without
changing the cellular controller or checkpoint layout. No ancestry gate, colony
reward or prescribed producer/consumer classes are present.
"""
mutable struct CatalyticHabitat{P}
    cells::CellularEcology
    substrates::ResourcePool
    enzymes::Vector{Array{Float64,3}}
    release::Vector{Array{Float64,3}}
    ids::Vector{Int}
    weights::Array{Float32,3}
    expression::Matrix{Float64}
    last_secretion::Matrix{Float64}
    secreted::Vector{Float64}
    solubilized::Vector{Float64}
    degraded::Vector{Float64}
    work_spent::Float64
    inert_residual::Float64
    seed::Int
    parameters::P
end

function CatalyticHabitat(w::World,cells::CellularEcology,substrates::ResourcePool;
        secretion_rate=.08,secretion_work=.2,allocation_exponent=2.,catalytic_rate=.35,
        half_saturation=.2,decay_rate=.015,investment_bias=-1.,allocation_bias_scale=1.2,
        mutation=.08,mutation_probability=.03,regulation=true,seed=w.seed)
    q=(;secretion_rate=Float64(secretion_rate),secretion_work=Float64(secretion_work),
        allocation_exponent=Float64(allocation_exponent),catalytic_rate=Float64(catalytic_rate),
        half_saturation=Float64(half_saturation),decay_rate=Float64(decay_rate),
        investment_bias=Float64(investment_bias),allocation_bias_scale=Float64(allocation_bias_scale),
        mutation=Float64(mutation),mutation_probability=Float64(mutation_probability),regulation=Bool(regulation))
    n=length(cells.ids); h=size(cells.hidden,1); rng=MersenneTwister(seed+831)
    weights=.12f0.*randn(rng,Float32,3,h+10,n)
    for i in 1:n
        weights[1,end,i]+=Float32(investment_bias)
        preference=Float32(allocation_bias_scale)*randn(rng,Float32)
        weights[2,end,i]+=preference; weights[3,end,i]-=preference
    end
    p=CatalyticHabitat(cells,substrates,[zeros(Float64,w.grid.size) for _ in 1:2],
        [zeros(Float64,w.grid.size) for _ in 1:2],copy(cells.ids),weights,zeros(3,n),zeros(2,n),
        zeros(2),zeros(2),zeros(2),0.,0.,Int(seed),q)
    validate(w,p);p
end

function validate(w::World,p::CatalyticHabitat)
    q=p.parameters;c=p.cells;f=p.substrates
    all(x->isfinite(x)&&x>=0,(q.secretion_rate,q.secretion_work,q.allocation_exponent,q.catalytic_rate,
        q.half_saturation,q.decay_rate,q.allocation_bias_scale,q.mutation,q.mutation_probability)) &&
        q.half_saturation>0 && q.allocation_exponent>=1 && q.mutation_probability<=1 &&
        isfinite(q.investment_bias) || throw(ArgumentError("invalid catalytic habitat parameters"))
    f!==c.resources && length(f.fields)==2 && all(i->f.potentials[i]>=c.resources.potentials[i],1:2) ||
        throw(ArgumentError("two distinct downhill feedstocks required"))
    validate(w,f);validate(w,c)
    length(p.enzymes)==length(p.release)==2 &&
        all(a->size(a)==w.grid.size && all(x->isfinite(x)&&x>=0,a),vcat(p.enzymes,p.release)) || throw(ArgumentError("invalid enzyme material or release field"))
    # An Any container preserves the actual Float32/Float64 array identities;
    # numeric array promotion during vcat would copy fields and hide aliases.
    inventories=Any[]
    for group in (f.fields,c.resources.fields,p.enzymes,p.release,values(w.fields));append!(inventories,group);end
    any(Base.mightalias(inventories[i],inventories[j]) for i in eachindex(inventories) for j in 1:i-1) &&
        throw(ArgumentError("catalytic reservoirs and observation buffers must not alias"))
    length(unique(p.ids))==length(p.ids) && size(p.expression)==(3,length(p.ids)) &&
        size(p.last_secretion)==(2,length(p.ids)) || throw(DimensionMismatch("catalytic identity/readout arrays"))
    size(p.weights)==(3,size(c.hidden,1)+10,length(p.ids)) && all(isfinite,p.weights) || throw(ArgumentError("invalid secretion genotype"))
    isfinite(p.inert_residual)&&p.inert_residual>=0 || throw(ArgumentError("invalid catalytic residual"))
    nothing
end
stability_limit(w,p::CatalyticHabitat)=Inf

function _catalytic_sync!(p)
    c=p.cells; p.ids==c.ids && return p
    old=Dict(id=>i for (i,id) in enumerate(p.ids)); n=length(c.ids)
    weights=zeros(Float32,size(p.weights,1),size(p.weights,2),n)
    for (j,id) in enumerate(c.ids)
        inherited=get(old,id,0); newborn=inherited==0
        if newborn
            inherited=get(old,c.parents[j],0)
            inherited>0 || throw(ArgumentError("secretion genotype missing for cell $id; schedule habitat every cellular step"))
        end
        weights[:,:,j].=view(p.weights,:,:,inherited)
        if newborn
            rng=MersenneTwister(string(p.seed,":catalytic-division:",id,":",c.parents[j]))
            for k in axes(weights,2),l in axes(weights,1)
                rand(rng)<p.parameters.mutation_probability && (weights[l,k,j]+=Float32(p.parameters.mutation)*randn(rng,Float32))
            end
        end
    end
    p.ids=copy(c.ids);p.weights=weights;p.expression=zeros(3,n);p.last_secretion=zeros(2,n)
    p
end

# Deposit only representable material and never charge a cell for a requested
# increment that the owned field did not receive.
function _catalytic_deposit!(field,stencil,i,amount)
    actual=0.
    for e in stencil.offsets[i]:stencil.offsets[i+1]-1
        site=stencil.sites[e]; want=min(amount-actual,amount*stencil.weights[e]);before=field[site]
        after=before+want;after-before>want && (after=prevfloat(after))
        gained=max(0.,after-before)
        while actual+gained>amount
            after=max(before,prevfloat(after));gained=after-before
        end
        field[site]=after;actual+=gained
    end
    actual
end

advance!(w::World,p::CatalyticHabitat,dt::Real)=_advance_catalytic!(w,p,dt,p.parameters,p.cells.parameters)
function _advance_catalytic!(w,p,dt,q,cq)
    check_dt(w,p,dt); delta=Float64(Float32(dt)); _catalytic_sync!(p)
    c=p.cells;f=p.substrates;r=c.resources;s=c.stencil;n=length(c.ids);h=size(c.hidden,1)
    fill!(p.last_secretion,0);foreach(a->fill!(a,0),p.release)
    input=zeros(Float64,h+10)
    for i in 1:n
        input[1:h].=view(c.hidden,:,i)
        fill!(view(input,h+1:h+4),0)
        for e in s.offsets[i]:s.offsets[i+1]-1
            site=s.sites[e];weight=s.weights[e]/s.cell_volume
            for k in 1:2
                input[h+k]+=weight*f.fields[k][site]
                input[h+2+k]+=weight*p.enzymes[k][site]
            end
        end
        for j in h+1:h+4;input[j]=input[j]/(1+input[j]);end
        input[h+5]=c.metabolites[1,i]/max(c.biomass[i]*cq.store_capacity,eps())
        input[h+6]=c.metabolites[2,i]/max(c.biomass[i]*cq.store_capacity,eps())
        input[h+7]=c.energy[i]/(c.biomass[i]+c.energy[i])
        input[h+8]=c.biomass[i]/(cq.division_mass+c.biomass[i])
        input[h+9]=c.age[i]/(20+c.age[i]);input[h+10]=1
        if q.regulation
            # Three bounded readouts; subtracting max keeps exponentials finite.
            output=[sum(Float64(p.weights[j,k,i])*input[k] for k in eachindex(input)) for j in 1:3]
            investment=1/(1+exp(-clamp(output[1],-40.,40.)))
            difference=clamp(output[2]-output[3],-40.,40.)
            allocation=1/(1+exp(-difference))
        else
            investment=1/(1+exp(-clamp(q.investment_bias,-40.,40.)));allocation=.5
        end
        p.expression[:,i].=(investment,allocation,1-allocation)
        for k in 1:2
            fraction=k==1 ? allocation : 1-allocation
            requested=delta*q.secretion_rate*c.biomass[i]*investment*fraction^q.allocation_exponent
            amount=min(requested,2minimum(view(c.metabolites,:,i)))
            if q.secretion_work>0
                affordable=c.energy[i]/q.secretion_work
                affordable*q.secretion_work>c.energy[i] && (affordable=prevfloat(affordable))
                amount=min(amount,max(0.,affordable))
            end
            actual=_catalytic_deposit!(p.enzymes[k],s,i,amount)
            c.metabolites[:,i].-=actual/2
            work=actual*q.secretion_work;c.energy[i]-=work;c.work_spent+=work;p.work_spent+=work
            c.chemical_heat+=actual/2*sum(Float64,view(r.potentials,1:2))
            p.secreted[k]+=actual;p.last_secretion[k,i]=actual
        end
    end
    for k in 1:2,site in eachindex(p.enzymes[k])
        old=f.fields[k][site]
        requested=delta*q.catalytic_rate*p.enzymes[k][site]*Float64(old)/(q.half_saturation*s.cell_volume+Float64(old))
        after=Float32(max(0.,Float64(old)-min(Float64(old),requested)))
        taken=Float64(old)-Float64(after);f.fields[k][site]=after;f.withdrawn[k]+=taken
        delivered=resource_add!(r,k,site,taken)
        p.inert_residual+=taken-delivered
        c.chemical_heat+=taken*Float64(f.potentials[k])-delivered*Float64(r.potentials[k])
        p.solubilized[k]+=delivered;p.release[k][site]=delivered/delta/s.cell_volume
        before=p.enzymes[k][site];next=before*exp(-q.decay_rate*delta)
        decayed=before-next;p.enzymes[k][site]=next
        waste=resource_add!(r,3,site,decayed)
        p.inert_residual+=decayed-waste;p.degraded[k]+=decayed
    end
    # Secreted material left the cell body; the next resource stage sees the new
    # cellular porosity. Enzymes are catalytic material, not an invented wall.
    cellular_porosity!(w,c)
    p
end

owned_material_mass(p::CatalyticHabitat)=sum(a->sum(a),p.enzymes)+p.inert_residual
composition_references(p::CatalyticHabitat)=(p.cells,p.substrates)
process_ports(::CatalyticHabitat)=ProcessPorts(reads=(:cellular_state,:feedstock,:resources),
    writes=(:feedstock,:resources,:cellular_inventory,:permeability),owned=(:enzymes,:weights,:inert_residual),
    memory="inherited secretion readout and persistent spatial catalyst material",
    law="paid complementary secretion; downhill solubilization; conservative enzyme decay into waste")
function catalytic_statistics(p::CatalyticHabitat)
    (enzyme_mass=map(sum,p.enzymes),secreted=copy(p.secreted),solubilized=copy(p.solubilized),
     degraded=copy(p.degraded),work_spent=p.work_spent,inert_residual=p.inert_residual,
     secretion_genotypes=size(p.weights,3),regulation=p.parameters.regulation)
end
function record_process_fields(p::CatalyticHabitat)
    (ids=copy(p.ids),weights=copy(p.weights),expression=copy(p.expression),last_secretion=copy(p.last_secretion),
     enzyme_a=copy(p.enzymes[1]),enzyme_b=copy(p.enzymes[2]),release_a=copy(p.release[1]),release_b=copy(p.release[2]),
     statistics=catalytic_statistics(p),parameters=p.parameters)
end
