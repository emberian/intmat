"""Paid catalyst doubles as a porous, reversibly sorptive material.

Food bound to the immobile scaffold is unavailable to cellular import and to the
resource pool's diffusion/outflow. Each food has a finite class of binding sites;
site capacity follows the actual enzyme inventory. Bound and free food carry the
same *bookkeeping* potential: this is a kinetic partition law, not a molecular
binding free-energy model or a source of work. No extra scaffold is invented.
"""
mutable struct SorptiveMatrix{P}
    habitat::CatalyticHabitat
    bound::Vector{Array{Float64,3}}
    permeability::Array{Float32,3}
    adsorbed::Vector{Float64}
    desorbed::Vector{Float64}
    displaced::Vector{Float64}
    inert_residual::Float64
    parameters::P
end

function SorptiveMatrix(w::World,h::CatalyticHabitat;capacity_per_enzyme=4.,
        adsorption_rate=2.,desorption_rate=.1,cross_binding=0.,transport_barrier=10.)
    q=(;capacity_per_enzyme=Float64(capacity_per_enzyme),adsorption_rate=Float64(adsorption_rate),
        desorption_rate=Float64(desorption_rate),cross_binding=Float64(cross_binding),
        transport_barrier=Float64(transport_barrier))
    p=SorptiveMatrix(h,[zeros(Float64,w.grid.size) for _ in 1:2],ones(Float32,w.grid.size),
        zeros(2),zeros(2),zeros(2),0.,q)
    validate(w,p);matrix_porosity!(w,p);p
end

function validate(w::World,p::SorptiveMatrix)
    validate(w,p.habitat);q=p.parameters;r=p.habitat.cells.resources
    all(x->isfinite(x)&&x>=0,values(q)) && q.cross_binding<=1 || throw(ArgumentError("invalid sorptive matrix parameters"))
    r.medium===nothing || throw(ArgumentError("matrix requires the cellular pool's combined permeability, not a bypassing medium reference"))
    length(p.bound)==2 && all(a->size(a)==w.grid.size && all(x->isfinite(x)&&x>=0,a),p.bound) || throw(ArgumentError("invalid bound inventories"))
    size(p.permeability)==w.grid.size && all(x->isfinite(x)&&0<=x<=1,p.permeability) || throw(ArgumentError("invalid matrix permeability"))
    isfinite(p.inert_residual)&&p.inert_residual>=0 || throw(ArgumentError("invalid sorptive residual"))
    f=p.habitat.substrates
    arrays=Any[p.bound...,p.permeability,r.permeability,r.scratch,f.permeability,f.scratch]
    for group in (r.fields,f.fields,p.habitat.enzymes,p.habitat.release,values(w.fields),values(w.scratch));append!(arrays,group);end
    any(Base.mightalias(arrays[i],arrays[j]) for i in eachindex(arrays) for j in 1:i-1) && throw(ArgumentError("sorptive inventories and buffers must not alias"))
    nothing
end
stability_limit(w,p::SorptiveMatrix)=Inf
composition_references(p::SorptiveMatrix)=(p.habitat,)
owned_material_mass(p::SorptiveMatrix)=sum(sum,p.bound)+p.inert_residual
matrix_potential(p::SorptiveMatrix)=sum(k->Float64(p.habitat.cells.resources.potentials[k])*sum(p.bound[k]),1:2)

# Exact finite-capacity reversible kinetics for fixed local total M and capacity C:
# dB/dt = (kon/V)(M-B)(C-B) - koff B. The stable small root avoids
# cancellation at low concentration; the repeated-root branch covers koff=0,M=C.
function sorption_target(B,M,C,kon,koff,volume,dt)
    C==0 && return 0.
    kon==0 && return B*exp(-koff*dt)
    a=kon/volume
    d=hypot(a*(M-C),sqrt(koff*(2a*(M+C)+koff)))
    d==0 && return C-(C-B)/(1+a*(C-B)*dt)
    small=2a*M*C/(a*(M+C)+koff+d)
    v=B-small;decay=exp(-d*dt)
    clamp(small+v*decay/(1-a*v*(-expm1(-d*dt))/d),0.,min(M,C))
end

# Actual representable field changes pay for transfers. Sub-ulp mass that cannot
# inhabit either endpoint is retained as inert material, with its potential heat.
function _sorptive_exchange!(p,k,site,target)
    r=p.habitat.cells.resources;old=p.bound[k][site];free=r.fields[k][site]
    residual=0.;moved=0.
    if target>old
        take=min(Float64(free),target-old)
        after=Float32(Float64(free)-take)
        Float64(free)-Float64(after)>take && (after=nextfloat(after))
        taken=Float64(free)-Float64(after)
        next=old+taken;next-old>taken && (next=prevfloat(next))
        gained=next-old;r.fields[k][site]=after;p.bound[k][site]=next
        r.withdrawn[k]+=taken;p.adsorbed[k]+=gained;residual=taken-gained;moved=gained
    elseif target<old
        delivered=resource_add!(r,k,site,old-target)
        next=max(0.,old-delivered)
        old-next<delivered && (next=max(0.,prevfloat(next)))
        lost=old-next;p.bound[k][site]=next;p.desorbed[k]+=delivered
        residual=lost-delivered;moved=-delivered
    end
    p.inert_residual+=residual
    p.habitat.cells.chemical_heat+=residual*Float64(r.potentials[k])
    moved
end

function advance!(w::World,p::SorptiveMatrix,dt::Real)
    check_dt(w,p,dt);q=p.parameters;h=p.habitat;r=h.cells.resources
    volume=h.cells.stencil.cell_volume;delta=Float64(Float32(dt))
    for k in 1:2,site in eachindex(p.bound[k])
        capacity=q.capacity_per_enzyme*((1-q.cross_binding)*h.enzymes[k][site]+q.cross_binding*h.enzymes[3-k][site])
        B=p.bound[k][site]
        if B>capacity
            p.displaced[k]-=_sorptive_exchange!(p,k,site,capacity)
            B=p.bound[k][site]
            if B>capacity
                # Scaffold loss is compulsory, unlike ordinary slow desorption.
                # If the soluble Float32 endpoint cannot receive the remainder,
                # preserve its mass as inert residue and dissipate its potential.
                lost=B-capacity;p.bound[k][site]=capacity;B=capacity
                p.inert_residual+=lost
                h.cells.chemical_heat+=lost*Float64(r.potentials[k])
            end
        end
        M=B+Float64(r.fields[k][site])
        target=sorption_target(B,M,capacity,q.adsorption_rate,q.desorption_rate,volume,delta)
        _sorptive_exchange!(p,k,site,target)
    end
    matrix_porosity!(w,p)
    p
end

function matrix_porosity!(w,p::SorptiveMatrix)
    h=p.habitat;q=p.parameters;r=h.cells.resources;volume=h.cells.stencil.cell_volume
    # Reconstruct body/background porosity before multiplication: repeated calls
    # never compound yesterday's resistance. This stage follows habitat secretion.
    cellular_porosity!(w,h.cells)
    for site in eachindex(p.permeability)
        resistance=1+q.transport_barrier*(h.enzymes[1][site]+h.enzymes[2][site])/volume
        p.permeability[site]=Float32(inv(resistance))
        r.permeability[site]*=p.permeability[site]
    end
    p
end

matrix_statistics(p::SorptiveMatrix)=(bound_mass=map(sum,p.bound),bound_potential=matrix_potential(p),
    adsorbed=copy(p.adsorbed),desorbed=copy(p.desorbed),capacity_released=copy(p.displaced),
    minimum_permeability=minimum(p.permeability),mean_permeability=sum(Float64,p.permeability)/length(p.permeability),
    inert_residual=p.inert_residual)
process_ports(::SorptiveMatrix)=ProcessPorts(reads=(:enzymes,:free_precursors,:cellular_porosity),
    writes=(:resources,:permeability),owned=(:bound,:permeability,:inert_residual),
    memory="immobile bound food and paid catalyst geometry",law="finite reversible sorption and concentration-dependent transport resistance")
record_process_fields(p::SorptiveMatrix)=(bound_a=copy(p.bound[1]),bound_b=copy(p.bound[2]),
    permeability=copy(p.permeability),statistics=matrix_statistics(p),parameters=p.parameters)

function sorptive_arena_snapshot!(fields,metrics,p,w,xs,ys,k)
    metrics["matrix"]=Dict(string(key)=>value for (key,value) in pairs(matrix_statistics(p)))
    for (name,arr,factor,units) in (("retained_a",p.bound[1],1/p.habitat.cells.stencil.cell_volume,"bound precursor mass / physical volume"),
            ("retained_b",p.bound[2],1/p.habitat.cells.stencil.cell_volume,"bound precursor mass / physical volume"),
            ("matrix_permeability",p.permeability,1.,"relative solute permeability"))
        fields[name]=Dict("shape"=>[length(xs),length(ys)],"values"=>vec(arr[xs,ys,k:k]).*factor,"units"=>units)
    end
    metrics["cellular_model"]="growing cells / inherited development / catalytic sorptive matrix"
    nothing
end
