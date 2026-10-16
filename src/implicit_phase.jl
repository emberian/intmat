"""CPU stabilized semi-implicit Cahn–Hilliard process with constant mobility.
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
    w.backend isa CPU || throw(ArgumentError("ImplicitPhaseSeparation currently requires CPU"))
    isfinite(stabilization) && stabilization>=0 || throw(ArgumentError("stabilization must be finite and nonnegative"))
    reference=PhaseSeparation(w;kwargs...)
    forcing=chemical_forcing===nothing ? zeros(Float32,w.grid.size) : Float32.(chemical_forcing)
    size(forcing)==w.grid.size || throw(DimensionMismatch("chemical forcing/world grid mismatch"))
    spectral_size=(w.grid.size[1]÷2+1,w.grid.size[2],w.grid.size[3])
    eigenvalues=zeros(Float32,spectral_size)
    for I in CartesianIndices(eigenvalues)
        eigenvalues[I]=sum(4sinpi((I[d]-1)/w.grid.size[d])^2/w.grid.spacing[d]^2 for d in 1:3 if w.grid.size[d]>1; init=0.0)
    end
    real_buffer=zeros(Float64,w.grid.size)
    spectrum=zeros(ComplexF64,spectral_size); nonlinear=similar(spectrum)
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
# FFTW plans retain dimensions/strides/native handles, not these Julia arrays.
struct ImplicitPhaseRuntime{F,C,R}
    forward::F
    inverse_complex::C
    inverse_real::R
    threads::Int
end
function ImplicitPhaseRuntime(p::ImplicitPhaseSeparation)
    forward=FFTW.plan_rfft(p.real_buffer;flags=FFTW.ESTIMATE)
    inverse_complex=FFTW.plan_ifft!(p.spectrum,(2,3);flags=FFTW.ESTIMATE)
    inverse_real=FFTW.plan_irfft(p.spectrum,size(p.real_buffer,1),(1,);flags=FFTW.ESTIMATE)
    ImplicitPhaseRuntime(forward,inverse_complex,inverse_real,Int(FFTW.get_num_threads()))
end
advance!(w::World,p::ImplicitPhaseSeparation,dt::Real)=advance!(w,p,dt,ImplicitPhaseRuntime(p))
function advance!(w::World,p::ImplicitPhaseSeparation,dt::Real,runtime::ImplicitPhaseRuntime)
    check_dt(w,p,dt)
    m=implicit_mobility(w,p); s=p.stabilization
    forward=runtime.forward; inverse_complex=runtime.inverse_complex; inverse_real=runtime.inverse_real
    for I in eachindex(p.phi); p.real_buffer[I]=p.phi[I]; end
    mul!(p.spectrum,forward,p.real_buffer)
    for I in eachindex(p.phi)
        u=Float64(p.phi[I])
        p.real_buffer[I]=p.bulk*(u^3-u)+p.chemical_forcing[I]
    end
    mul!(p.nonlinear_spectrum,forward,p.real_buffer)
    for I in eachindex(p.spectrum)
        q=Float64(p.eigenvalues[I]); a=Float64(dt)*m*q
        p.spectrum[I]=((1+a*s)*p.spectrum[I]-a*p.nonlinear_spectrum[I])/(1+a*(s+p.kappa*q))
    end
    inverse_complex*p.spectrum
    mul!(p.real_buffer,inverse_real,p.spectrum)
    for I in eachindex(p.phi); p.phi[I]=p.real_buffer[I]; end
    phase_observe!(w,p)
    launch!(w,phase_couple_kernel!,w.scratch.nutrient,w.fields.nutrient,w.fields.activity,p.phi,p.interface,w.grid.spacing,Float32(dt),p.nutrient_diffusion,p.phase_barrier,p.interface_activity,p.activity_relaxation)
    copyto!(w.fields.nutrient,w.scratch.nutrient)
    w
end
phase_energy(w::World,p::ImplicitPhaseSeparation)=phase_energy(w,p.reference)
phase_optical_fields(w::World,p::ImplicitPhaseSeparation)=phase_optical_fields(w,p.reference)
record_process_fields(p::ImplicitPhaseSeparation)=merge(record_process_fields(p.reference),(chemical_forcing=copy(p.chemical_forcing),stabilization=p.stabilization))
