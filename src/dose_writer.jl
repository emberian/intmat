"""Opt-in isolated thermal writing with dose as the sole material state.

The thermal update and exposure quadrature are identical to ThermalWriter.
Uniform initial material is reconstructed as nutrient=n₀ exp(-conversion*dose)
and deposit=d₀+n₀*(-expm1(-conversion*dose)). This preserves the analytic
reaction invariant, with different Float32 rounding from incremental exchange.
This process owns its material; it does not read or update World nutrient/deposit.
External feeding, erosion and material feedback need a different carrier.

Temperature uses ping-pong buffers: `p.writer` and its temperature/scratch
arrays are borrowed step-local handles. Reacquire them after each advance;
copy a field to retain a snapshot. Dose keeps its array identity. The original
ThermalWriter retains its stable-buffer contract.
"""
mutable struct DoseThermalWriter{W}
    writer::W
    alternate::W
    initial_nutrient::Float32
    initial_deposit::Float32
end
function DoseThermalWriter(w::World;initial_nutrient=1f0,initial_deposit=0f0,kwargs...)
    q=ThermalWriter(w;kwargs...)
    alternate=ThermalWriter(q.scratch,q.temperature,q.dose,q.path,q.center,q.radius,q.height,q.period,q.width,q.heat,q.diffusion,q.cooling,q.threshold,q.conversion,q.pulse_rate,q.start,q.stop)
    p=DoseThermalWriter(q,alternate,Float32(initial_nutrient),Float32(initial_deposit))
    validate(w,p);p
end
process_children(p::DoseThermalWriter)=(p.writer,p.alternate)
function validate(w::World,p::DoseThermalWriter)
    validate(w,p.writer)
    all(x->isfinite(x)&&x>=0,(p.initial_nutrient,p.initial_deposit)) || throw(ArgumentError("initial dose-writer inventory must be finite and nonnegative"))
end
stability_limit(w,p::DoseThermalWriter)=stability_limit(w,p.writer)
@kernel function write_dose_kernel!(out,temp,dose,h,extent,dt,position,width,power,D,cooling,threshold)
    i,j,k=@index(Global,NTuple)
    T,exposure=thermal_exposure(temp,i,j,k,h,extent,dt,position,width,power,D,cooling,threshold)
    dose[i,j,k]+=exposure
    out[i,j,k]=T
end
function advance!(w::World,p::DoseThermalWriter,dt::Real)
    check_dt(w,p,dt);q=p.writer;t=w.time+Float64(dt)/2
    active=q.start<=t<q.stop
    envelope=q.pulse_rate==0 ? 1f0 : Float32(.5*(1+cos(2pi*Float64(q.pulse_rate)*t)))
    power=active ? q.heat*envelope : 0f0
    launch!(w,write_dose_kernel!,q.scratch,q.temperature,q.dose,w.grid.spacing,w.grid.extent,Float32(dt),beam_position(q,t),q.width,power,q.diffusion,q.cooling,q.threshold)
    p.writer,p.alternate=p.alternate,p.writer
    w
end
@kernel function dose_material_kernel!(nutrient,deposit,dose,n0,d0,conversion)
    i=@index(Global,Linear)
    fraction=-expm1(-conversion*dose[i])
    deposit[i]=d0+n0*fraction
    nutrient[i]=n0*exp(-conversion*dose[i])
end
"""Reconstruct material into caller-owned backend arrays, e.g. bounded output tiles."""
function writer_material!(nutrient,deposit,p::DoseThermalWriter;dose=p.writer.dose)
    q=p.writer
    size(nutrient)==size(deposit)==size(dose) || throw(DimensionMismatch("dose material output grid"))
    (Base.mightalias(nutrient,deposit) || any(a->Base.mightalias(nutrient,a)||Base.mightalias(deposit,a),(dose,q.temperature,q.scratch,q.dose))) && throw(ArgumentError("material readout arrays must not alias each other or any owned thermal/dose state"))
    backend=KernelAbstractions.get_backend(dose)
    all(a->KernelAbstractions.get_backend(a)==backend,(nutrient,deposit,dose)) || throw(ArgumentError("dose material output backend"))
    dose_material_kernel!(backend)(nutrient,deposit,dose,p.initial_nutrient,p.initial_deposit,q.conversion;ndrange=length(dose))
    KernelAbstractions.synchronize(backend);(;nutrient,deposit)
end
record_process_fields(p::DoseThermalWriter)=merge(record_process_fields(p.writer),(initial_nutrient=p.initial_nutrient,initial_deposit=p.initial_deposit,material_representation="analytic uniform-stock reaction from dose"))

# The private reservoir never exchanges material with World fields.
owned_material_mass(p::DoseThermalWriter)=(Float64(p.initial_nutrient)+Float64(p.initial_deposit))*length(p.writer.dose)
process_ports(::DoseThermalWriter)=ProcessPorts(owned=(:writer,),drives="prescribed Gaussian heat source; private fixed uniform material reservoir",memory="temperature and dose",law="shared explicit thermal step; analytic private inventory readout")
