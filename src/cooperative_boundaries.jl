"""A funded, decaying local film coupled through one `ActiveMedium` stage.

Film is distinct from mineral deposit. Maintenance callers must debit exactly
what `maintain_boundary!` returns from their owned material and work reservoirs.
Advance this wrapper instead of also advancing its referenced medium.
"""
mutable struct BoundaryMaintenance{M,A}
    medium::M
    film::A
    work_per_mass::Float32
    decay::Float32
    contributed_material::Float64
    spent_work::Float64
    recycled_material::Float64
end
function BoundaryMaintenance(w::World,medium::ActiveMedium; work_per_mass=1,decay=.1)
    b=BoundaryMaintenance(medium,zeros(Float32,w.grid.size),Float32(work_per_mass),
        Float32(decay),0.0,0.0,0.0)
    validate(w,b)
    b
end
function validate(w::World,b::BoundaryMaintenance)
    validate(w,b.medium)
    size(b.film)==w.grid.size || throw(DimensionMismatch("film/world grid mismatch"))
    isfinite(b.work_per_mass) && b.work_per_mass>0 || throw(ArgumentError("work_per_mass must be positive and finite"))
    isfinite(b.decay) && b.decay>=0 || throw(ArgumentError("decay must be finite and nonnegative"))
    all(x->isfinite(x)&&x>=0,b.film) || throw(ArgumentError("film must be finite and nonnegative"))
    nothing
end
"""Allocate affordable film and return exact caller debits (material, work).

No resource is minted: the caller supplies available budgets, then transfers the
returned material from its reservoir. Zero requests implement withholding.
"""
function maintain_boundary!(b::BoundaryMaintenance,cell::NTuple{3,Int},material_available::Real,
        work_available::Real,requested_mass::Real)
    all(x->isfinite(x)&&x>=0,(material_available,work_available,requested_mass)) ||
        throw(ArgumentError("maintenance budgets must be finite and nonnegative"))
    checkbounds(b.film,cell...)
    # Round downward so Float32 conversion cannot overdraw either budget.
    cap=min(Float64(material_available),Float64(work_available)/b.work_per_mass,Float64(requested_mass))
    amount=Float32(cap)
    Float64(amount)>cap && (amount=prevfloat(amount))
    old=b.film[cell...]
    updated=old+amount
    isfinite(updated) || throw(ArgumentError("film overflow"))
    actual=max(0f0,updated-old)
    Float64(actual)>cap && (updated=prevfloat(updated); actual=max(0f0,updated-old))
    work=Float64(actual)*b.work_per_mass
    b.film[cell...]=updated
    b.contributed_material+=actual
    b.spent_work+=work
    (material=Float64(actual),work=work)
end
boundary_material_mass(b::BoundaryMaintenance)=sum(Float64,b.film)
"""Local public-good exposure, including adjacent beneficiaries across periodic edges.

This is a physical exposure diagnostic, not evidence of cooperative fitness.
Transport itself uses symmetric faces; a film cell slows both incident faces.
"""
function local_boundary_benefit(b::BoundaryMaintenance,cell::NTuple{3,Int})
    checkbounds(b.film,cell...)
    exposure=Float64(b.film[cell...]); count=1
    for axis in 1:3
        if size(b.film,axis)>1
            for side in (-1,1)
                exposure+=b.film[phase_neighbor(cell,axis,side,size(b.film))...]
                count+=1
            end
        end
    end
    exposure/count
end
function advance!(w::World,b::BoundaryMaintenance,dt::Real)
    check_dt(w,b,dt)
    survival=exp(-Float64(b.decay)*Float64(dt))
    for idx in eachindex(b.film)
        old=b.film[idx]
        b.film[idx]=Float32(old*survival)
        recycled=old-b.film[idx]
        w.fields.nutrient[idx]+=recycled
        b.recycled_material+=recycled
    end
    # Transactional projection lets the existing medium operator use total
    # boundary material for transport and pinning. Restore mineral deposit
    # exactly (not by subtracting floats) before any external observer runs.
    deposit=copy(b.medium.membrane)
    b.medium.membrane .+= b.film
    try
        advance!(w,b.medium,dt)
    finally
        copyto!(b.medium.membrane,deposit)
    end
    w
end
function record_process_fields(b::BoundaryMaintenance)
    medium=record_process_fields(b.medium)
    names=Tuple(Symbol("medium_",n) for n in keys(medium))
    recorded_values=Tuple(v isa AbstractArray ? copy(v) : v for v in Base.values(medium))
    merge(NamedTuple{names}(recorded_values),(film=copy(b.film),permeability=copy(b.medium.permeability),signal=copy(b.medium.signal),
     work_per_mass=b.work_per_mass,decay=b.decay,
     contributed_material=b.contributed_material,spent_work=b.spent_work,
     recycled_material=b.recycled_material,
     medium_reference="wrapped ActiveMedium; advance only BoundaryMaintenance"))
end

# The wrapper owns advancement of its medium; the phase is a borrowed reference.
process_children(b::BoundaryMaintenance)=(b.medium,)
