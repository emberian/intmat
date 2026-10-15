"""Cells and their construction effector share one anabolic material/work budget.

Unlike the legacy `(cells, habitat)` schedule, secretion is not restricted to
whatever precursor growth leaves behind. Desired growth and both construction
fluxes are scaled proportionally by the same limiting precursor and work. The
normal cellular path remains available with unchanged arithmetic.
"""
mutable struct ConstructiveEcology{C,H,M}
    cells::C
    habitat::H
    mechanics::M
    growth_requested::Float64
    construction_requested::Vector{Float64}
    growth_planned::Float64
    construction_planned::Vector{Float64}
end
function ConstructiveEcology(w::World,c::CellularEcology,h::CatalyticHabitat;mechanics=nothing)
    mechanics===nothing || mechanics isa MatrixDrag || throw(ArgumentError("unsupported constructive mechanics law"))
    p=ConstructiveEcology(c,h,mechanics,0.,zeros(2),0.,zeros(2))
    validate(w,p);p
end
function validate(w::World,p::ConstructiveEcology)
    p.habitat.cells===p.cells || throw(ArgumentError("construction and cellular owners must share the same cells"))
    validate(w,p.cells);validate(w,p.habitat)
    _validate_matrix_mechanics(w,p.habitat,p.mechanics)
    all(x->isfinite(x)&&x>=0,(p.growth_requested,p.growth_planned,p.construction_requested...,p.construction_planned...)) ||
        throw(ArgumentError("invalid anabolic diagnostics"))
    nothing
end
stability_limit(w,p::ConstructiveEcology)=min(stability_limit(w,p.cells),stability_limit(w,p.habitat))
process_children(p::ConstructiveEcology)=(p.cells,p.habitat)

struct CatalyticAnabolism{P,D}
    owner::P
    input::Vector{Float64}
    environment::Matrix{Float64}
    requests::Matrix{Float64}
    drag::D
    drag_ids::Vector{Int}
end

function CatalyticAnabolism(p)
    h=p.habitat;c=p.cells;s=c.stencil;n=length(c.ids)
    environment=zeros(4,n)
    # All cells observe the same pre-construction field, not additions by earlier
    # IDs during this anabolic sweep. Intracellular state is read at planning.
    for i in 1:n,e in s.offsets[i]:s.offsets[i+1]-1
        site=s.sites[e];weight=s.weights[e]/s.cell_volume
        for k in 1:2
            environment[k,i]+=weight*h.substrates.fields[k][site]
            environment[k+2,i]+=weight*h.enzymes[k][site]
        end
    end
    drag=_matrix_drag_coefficients(p.mechanics,environment)
    for j in eachindex(environment);environment[j]/=1+environment[j];end
    CatalyticAnabolism(p,zeros(Float64,size(c.hidden,1)+10),environment,zeros(2,n),drag,drag===nothing ? Int[] : copy(c.ids))
end

function _cellular_planned_work(planned,w,c,i,delta,cq,a::CatalyticAnabolism)
    p=a.owner;h=p.habitat;q=h.parameters
    investment,allocation=_catalytic_readout!(h,i,a.input,q,cq,a.environment)
    request_a=delta*q.secretion_rate*c.biomass[i]*investment*allocation^q.allocation_exponent
    request_b=delta*q.secretion_rate*c.biomass[i]*investment*(1-allocation)^q.allocation_exponent
    a.requests[:,i].=(request_a,request_b)
    _matrix_drag_forecast(planned+q.secretion_work*(request_a+request_b),w,c,i,delta,cq,a.drag)
end

function _cellular_mechanical_drag(w,c,a::CatalyticAnabolism)
    a.drag===nothing && return nothing
    c.ids==a.drag_ids && return a.drag
    # Terminal cells leave before mechanics; stable IDs preserve each survivor's
    # earlier material observation without assigning a dead cell's drag to it.
    index=Dict(id=>i for (i,id) in enumerate(a.drag_ids))
    [a.drag[index[id]] for id in c.ids]
end

function _cellular_anabolic_build!(w,c,i,delta,growth_request,cq,a::CatalyticAnabolism)
    p=a.owner;h=p.habitat;q=h.parameters
    request_a=a.requests[1,i];request_b=a.requests[2,i]
    total=growth_request+request_a+request_b
    work=cq.growth_work*growth_request+q.secretion_work*(request_a+request_b)
    available=2minimum(view(c.metabolites,:,i))
    all(x->isfinite(x)&&x>=0,(growth_request,request_a,request_b,total,work,available,c.energy[i])) ||
        throw(ArgumentError("shared anabolism requires finite nonnegative requests and donors"))
    scale=min(1.,total>0 ? available/total : 1.,work>0 ? c.energy[i]/work : 1.)
    growth=growth_request*scale;first=request_a*scale;second=request_b*scale
    # Round the shared plan toward both donor budgets before any effector acts.
    for _ in 1:8
        growth+first+second<=available && cq.growth_work*growth+q.secretion_work*(first+second)<=c.energy[i] && break
        scale=max(0.,prevfloat(scale));growth=growth_request*scale;first=request_a*scale;second=request_b*scale
    end
    if growth+first+second>available || cq.growth_work*growth+q.secretion_work*(first+second)>c.energy[i]
        growth=first=second=0.
    end
    p.growth_requested+=growth_request;p.construction_requested[1]+=request_a;p.construction_requested[2]+=request_b
    p.growth_planned+=growth;p.construction_planned[1]+=first;p.construction_planned[2]+=second
    _cellular_anabolic_build!(w,c,i,delta,growth,cq,nothing)
    # These transfers retain ordinary donor-bounded representable deposition;
    # any unrepresentable portion stays intracellular, not a fabricated refund.
    _catalytic_construct!(h,i,1,first,q)
    _catalytic_construct!(h,i,2,second,q)
    nothing
end

function advance!(w::World,p::ConstructiveEcology,dt::Real)
    check_dt(w,p,dt);c=p.cells;h=p.habitat
    issorted(c.ids) || reorder_ecology!(c,sortperm(c.ids))
    _catalytic_sync!(h);fill!(h.last_secretion,0);foreach(a->fill!(a,0),h.release)
    allocator=CatalyticAnabolism(p)
    ecology_step!(w,c,dt;anabolism=allocator)
    _catalytic_sync!(h;preserve_observation=true)
    _catalytic_react!(w,h,Float64(Float32(dt)),h.parameters)
    p
end

anabolic_statistics(p::ConstructiveEcology)=(mode="shared finite material/work budget",
    growth_requested=p.growth_requested,construction_requested=copy(p.construction_requested),
    growth_planned=p.growth_planned,construction_planned=copy(p.construction_planned))
process_ports(::ConstructiveEcology)=ProcessPorts(reads=(:resources,:feedstock,:enzymes,:cellular_state),
    writes=(:resources,:feedstock,:cellular_inventory,:permeability),
    memory="owned cellular development and inherited construction readout",
    law="maintenance followed by proportional complementary-material/work allocation to growth and construction")
record_process_fields(p::ConstructiveEcology)=(cells=record_process_fields(p.cells),
    habitat=record_process_fields(p.habitat),anabolism=anabolic_statistics(p))
