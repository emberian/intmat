"""Stabilized semi-implicit Cahn–Hilliard process with constant mobility.
The periodic finite-difference Laplacian is diagonalized by FFT, preserving the
explicit reference's spatial discretization. Mineral pinning is supported for
uniform deposits; varying deposits require the explicit reference process.
`chemical_forcing` is an externally writable chemical-potential field, never a
composition source: its Laplacian conserves the signed composition.
"""
struct ImplicitPhaseSeparation{P,F,R,C}
    reference::P
    chemical_forcing::F
    stabilization::Float64
    eigenvalues::F
    real_buffer::R
    spectrum::C
    nonlinear_spectrum::C
end
function Base.getproperty(p::ImplicitPhaseSeparation,n::Symbol)
    n in fieldnames(typeof(p)) ? getfield(p,n) : getproperty(getfield(p,:reference),n)
end
Base.propertynames(p::ImplicitPhaseSeparation)= (fieldnames(typeof(p))...,propertynames(p.reference)...)
function ImplicitPhaseSeparation(w::World; stabilization=2, chemical_forcing=nothing, kwargs...)
    KernelAbstractions.supports_float64(w.backend) || throw(ArgumentError("ImplicitPhaseSeparation requires Float64 FFT support on its backend; reduced precision is not implicit"))
    isfinite(stabilization) && stabilization>=0 || throw(ArgumentError("stabilization must be finite and nonnegative"))
    chemical_forcing===nothing || size(chemical_forcing)==w.grid.size || throw(DimensionMismatch("chemical forcing/world grid mismatch"))
    reference=PhaseSeparation(w;kwargs...)
    forcing=KernelAbstractions.zeros(w.backend,Float32,w.grid.size...)
    chemical_forcing===nothing || copyto!(forcing,Float32.(chemical_forcing))
    size(forcing)==w.grid.size || throw(DimensionMismatch("chemical forcing/world grid mismatch"))
    spectral_size=(w.grid.size[1]÷2+1,w.grid.size[2],w.grid.size[3])
    eigenvalues=zeros(Float32,spectral_size)
    for I in CartesianIndices(eigenvalues)
        eigenvalues[I]=sum(4sinpi((I[d]-1)/w.grid.size[d])^2/w.grid.spacing[d]^2 for d in 1:3 if w.grid.size[d]>1; init=0.0)
    end
    if !(w.backend isa CPU)
        host=eigenvalues
        eigenvalues=KernelAbstractions.zeros(w.backend,Float32,spectral_size...)
        copyto!(eigenvalues,host)
    end
    real_buffer=KernelAbstractions.zeros(w.backend,Float64,w.grid.size...)
    spectrum=KernelAbstractions.zeros(w.backend,ComplexF64,spectral_size...); nonlinear=similar(spectrum)
    p=ImplicitPhaseSeparation(reference,forcing,Float64(stabilization),eigenvalues,real_buffer,spectrum,nonlinear)
    validate(w,p)
    p
end
function implicit_mobility(w,p)
    p.mineral_pinning==0 && return Float64(p.mobility)
    lo,hi=extrema(w.fields.deposit)
    lo==hi || throw(ArgumentError("implicit phase requires uniform deposit with mineral_pinning>0; use PhaseSeparation for spatially varying pinned mobility or explicitly choose mineral_pinning=0"))
    isfinite(lo) && lo>=0 || throw(ArgumentError("deposit must be finite and nonnegative"))
    Float64(p.mobility)/(1+Float64(p.mineral_pinning)*lo)
end
function validate(w::World,p::ImplicitPhaseSeparation)
    validate(w,p.reference)
    all(isfinite,p.phi) || throw(ArgumentError("phase composition must be finite"))
    size(p.chemical_forcing)==w.grid.size || throw(DimensionMismatch("chemical forcing/world grid mismatch"))
    all(isfinite,p.chemical_forcing) || throw(ArgumentError("chemical forcing must be finite"))
    implicit_mobility(w,p)
    nothing
end
# Only explicit nutrient transport/activity feedback retains a timestep bound.
function stability_limit(w::World,p::ImplicitPhaseSeparation)
    min(diffusion_limit(w,p.nutrient_diffusion),p.activity_relaxation==0 ? Inf : inv(Float64(p.activity_relaxation)))
end
# Optional worker-owned native resources stay outside the physical state.
abstract type AbstractImplicitPhaseRuntime end
implicit_phase_layout(p)=(backend=KernelAbstractions.get_backend(p.real_buffer),
    buffers=(objectid(p.real_buffer),objectid(p.spectrum)),
    real=(size(p.real_buffer),strides(p.real_buffer),eltype(p.real_buffer)),
    complex=(size(p.spectrum),strides(p.spectrum),eltype(p.spectrum)))
struct ImplicitPhaseRuntime{F,C,R,L} <: AbstractImplicitPhaseRuntime
    forward::F
    inverse_complex::C
    inverse_real::R
    threads::Int
    layout::L
end
struct DeviceImplicitPhaseRuntime{F,R,L} <: AbstractImplicitPhaseRuntime
    forward::F
    inverse_real::R
    normalization::Float64
    layout::L
end
ImplicitPhaseRuntime(p::ImplicitPhaseSeparation)=implicit_phase_runtime(p)
# Shared resource factory also supports explicitly separate experimental states.
function implicit_phase_runtime(p)
    backend=KernelAbstractions.get_backend(p.real_buffer)
    if backend isa CPU
        forward=FFTW.plan_rfft(p.real_buffer;flags=FFTW.ESTIMATE)
        # The separable CPU inverse avoids a full-grid planning copy.
        inverse_complex=FFTW.plan_ifft!(p.spectrum,(2,3);flags=FFTW.ESTIMATE)
        inverse_real=FFTW.plan_irfft(p.spectrum,size(p.real_buffer,1),(1,);flags=FFTW.ESTIMATE)
        return ImplicitPhaseRuntime(forward,inverse_complex,inverse_real,Int(FFTW.get_num_threads()),implicit_phase_layout(p))
    end
    KernelAbstractions.supports_float64(backend) || throw(ArgumentError("implicit phase runtime requires Float64 support"))
    forward=FFTW.plan_rfft(p.real_buffer)
    # Unnormalized inverse plus Julia scaling does not require a device BLAS.
    inverse_real=FFTW.plan_brfft(p.spectrum,size(p.real_buffer,1))
    DeviceImplicitPhaseRuntime(forward,inverse_real,inv(Float64(length(p.real_buffer))),implicit_phase_layout(p))
end
@inline implicit_phase_nonlinear(u,forcing,bulk)=bulk*(u^3-u)+forcing
@inline function implicit_phase_mode(value,nonlinear,q,dt,m,s,kappa)
    a=Float64(dt)*m*q
    ((1+a*s)*value-a*nonlinear)/(1+a*(s+kappa*q))
end
@kernel function implicit_phase_real_kernel!(out,phi)
    i=@index(Global,Linear)
    out[i]=Float64(phi[i])
end
@kernel function implicit_phase_nonlinear_kernel!(out,phi,forcing,bulk)
    i=@index(Global,Linear)
    out[i]=implicit_phase_nonlinear(Float64(phi[i]),forcing[i],bulk)
end
@kernel function implicit_phase_mode_kernel!(spectrum,nonlinear,eigenvalues,dt,m,s,kappa)
    i=@index(Global,Linear)
    spectrum[i]=implicit_phase_mode(spectrum[i],nonlinear[i],Float64(eigenvalues[i]),dt,m,s,kappa)
end
@kernel function implicit_phase_finish_kernel!(phi,real_buffer,normalization)
    i=@index(Global,Linear)
    phi[i]=real_buffer[i]*normalization
end
function implicit_phase_real!(w,p,::ImplicitPhaseRuntime)
    for I in eachindex(p.phi);p.real_buffer[I]=p.phi[I];end
end
function implicit_phase_nonlinear!(w,p,::ImplicitPhaseRuntime)
    for I in eachindex(p.phi)
        p.real_buffer[I]=implicit_phase_nonlinear(Float64(p.phi[I]),p.chemical_forcing[I],p.bulk)
    end
end
function implicit_phase_modes!(w,p,dt,m,s,::ImplicitPhaseRuntime)
    for I in eachindex(p.spectrum)
        p.spectrum[I]=implicit_phase_mode(p.spectrum[I],p.nonlinear_spectrum[I],Float64(p.eigenvalues[I]),dt,m,s,p.kappa)
    end
end
function implicit_phase_inverse!(w,p,r::ImplicitPhaseRuntime)
    r.inverse_complex*p.spectrum
    mul!(p.real_buffer,r.inverse_real,p.spectrum)
    for I in eachindex(p.phi);p.phi[I]=p.real_buffer[I];end
end
function implicit_phase_launch!(w,kernel,args...;ndrange)
    kernel(w.backend,256)(args...;ndrange)
    nothing
end
implicit_phase_real!(w,p,::DeviceImplicitPhaseRuntime)=implicit_phase_launch!(w,implicit_phase_real_kernel!,p.real_buffer,p.phi;ndrange=length(p.phi))
implicit_phase_nonlinear!(w,p,::DeviceImplicitPhaseRuntime)=implicit_phase_launch!(w,implicit_phase_nonlinear_kernel!,p.real_buffer,p.phi,p.chemical_forcing,p.bulk;ndrange=length(p.phi))
implicit_phase_modes!(w,p,dt,m,s,::DeviceImplicitPhaseRuntime)=implicit_phase_launch!(w,implicit_phase_mode_kernel!,p.spectrum,p.nonlinear_spectrum,p.eigenvalues,dt,m,s,p.kappa;ndrange=length(p.spectrum))
function implicit_phase_inverse!(w,p,r::DeviceImplicitPhaseRuntime)
    mul!(p.real_buffer,r.inverse_real,p.spectrum)
    implicit_phase_launch!(w,implicit_phase_finish_kernel!,p.phi,p.real_buffer,r.normalization;ndrange=length(p.phi))
end
function advance!(w::World,p::ImplicitPhaseSeparation,dt::Real)
    check_dt(w,p,dt) # Reject invalid steps before allocating native plans.
    implicit_phase_advance!(w,p,dt,ImplicitPhaseRuntime(p))
end
function advance!(w::World,p::ImplicitPhaseSeparation,dt::Real,runtime::AbstractImplicitPhaseRuntime)
    check_dt(w,p,dt)
    implicit_phase_advance!(w,p,dt,runtime)
end
function implicit_phase_advance!(w,p,dt,runtime)
    runtime.layout==implicit_phase_layout(p) || throw(ArgumentError("implicit phase runtime backend/array layout differs"))
    m=implicit_mobility(w,p);s=p.stabilization
    implicit_phase_real!(w,p,runtime)
    mul!(p.spectrum,runtime.forward,p.real_buffer)
    implicit_phase_nonlinear!(w,p,runtime)
    mul!(p.nonlinear_spectrum,runtime.forward,p.real_buffer)
    implicit_phase_modes!(w,p,dt,m,s,runtime)
    implicit_phase_inverse!(w,p,runtime)
    phase_observe!(w,p)
    launch!(w,phase_couple_kernel!,w.scratch.nutrient,w.fields.nutrient,w.fields.activity,p.phi,p.interface,w.grid.spacing,Float32(dt),p.nutrient_diffusion,p.phase_barrier,p.interface_activity,p.activity_relaxation)
    copyto!(w.fields.nutrient,w.scratch.nutrient)
    w
end
phase_energy(w::World,p::ImplicitPhaseSeparation)=phase_energy(w,p.reference)
phase_optical_fields(w::World,p::ImplicitPhaseSeparation)=phase_optical_fields(w,p.reference)
record_process_fields(p::ImplicitPhaseSeparation)=merge(record_process_fields(p.reference),(chemical_forcing=copy(p.chemical_forcing),stabilization=p.stabilization))

# Physical transfer is explicit because raw GPUArrays serialization does not
# preserve repeated device-array aliases. One memo spans world and process.
function implicit_phase_transfer(w::World,p::ImplicitPhaseSeparation,backend)
    KernelAbstractions.supports_float64(backend) || throw(ArgumentError("implicit phase transfer requires Float64 support"))
    memo=IdDict()
    transfer(a)=get!(memo,a) do
        out=KernelAbstractions.zeros(backend,eltype(a),size(a)...)
        copyto!(out,a)
        out
    end
    world=World(w.grid,backend,map(transfer,w.fields),map(transfer,w.scratch),w.time,w.seed)
    r=p.reference
    reference=PhaseSeparation(transfer(r.phi),transfer(r.chemical_potential),transfer(r.interface),transfer(r.scratch),
        r.bulk,r.kappa,r.mobility,r.mineral_pinning,r.nutrient_diffusion,r.phase_barrier,r.interface_activity,r.activity_relaxation)
    phase=ImplicitPhaseSeparation(reference,transfer(p.chemical_forcing),p.stabilization,transfer(p.eigenvalues),
        transfer(p.real_buffer),transfer(p.spectrum),transfer(p.nonlinear_spectrum))
    validate(world,phase)
    world,phase
end
const IMPLICIT_PHASE_CHECKPOINT_SOURCE=(implicit=bytes2hex(sha256(read(@__FILE__))),
    phase=bytes2hex(sha256(read(joinpath(@__DIR__,"phase.jl")))),core=bytes2hex(sha256(read(joinpath(@__DIR__,"core.jl")))))
function implicit_phase_checkpoint(w::World,p::ImplicitPhaseSeparation)
    validate(w,p)
    world,phase=implicit_phase_transfer(w,p,CPU())
    (;format="Interior.ImplicitPhase",version=1,julia_version=VERSION,source=IMPLICIT_PHASE_CHECKPOINT_SOURCE,world,phase)
end
function implicit_phase_restore(payload;backend=CPU())
    payload.format=="Interior.ImplicitPhase"&&payload.version==1&&payload.julia_version==VERSION&&payload.source==IMPLICIT_PHASE_CHECKPOINT_SOURCE ||
        throw(ArgumentError("implicit phase checkpoint source/schema differs"))
    implicit_phase_transfer(payload.world,payload.phase,backend)
end
