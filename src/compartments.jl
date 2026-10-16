"""An active, periodic compartment environment referencing a separately advanced phase.

Membrane material aliases `world.fields.deposit`: inhabitants must pay for builds
from their cargo, and recover eroded material themselves. `signal_source` and
`remodeling` are per-stage increments, consumed by `advance!`. Nutrient exchange
uses symmetric face fluxes and internally substeps to preserve positivity.
"""
struct ActiveMedium{P,A}
    phase::P
    membrane::A
    permeability::A
    signal::A
    signal_source::A
    remodeling::A
    forcing::A
    scratch::A
    nutrient_diffusion::Float32
    signal_diffusion::Float32
    signal_decay::Float32
    interface_barrier::Float32
    membrane_barrier::Float32
    signal_opening::Float32
    remodeling_gain::Float32
    forcing_relaxation::Float32
    membrane_pinning::Float32
    signal_bias::Float32
end
function ActiveMedium(w::World,phase; nutrient_diffusion=.03, signal_diffusion=.02,
        signal_decay=.4, interface_barrier=8, membrane_barrier=6, signal_opening=2,
        remodeling_gain=.3, forcing_relaxation=1, membrane_pinning=.5, signal_bias=.1)
    w.backend isa CPU || throw(ArgumentError("ActiveMedium currently requires CPU"))
    hasproperty(phase,:chemical_forcing) || throw(ArgumentError("phase must expose chemical_forcing"))
    arrays=ntuple(_->zeros(Float32,w.grid.size),6)
    params=Float32.((nutrient_diffusion,signal_diffusion,signal_decay,interface_barrier,
        membrane_barrier,signal_opening,remodeling_gain,forcing_relaxation,membrane_pinning,signal_bias))
    m=ActiveMedium(phase,w.fields.deposit,arrays...,params...)
    validate(w,m)
    medium_observe!(m)
    m
end
function validate(w::World,m::ActiveMedium)
    w.backend isa CPU || throw(ArgumentError("ActiveMedium currently requires CPU"))
    all(x->isfinite(x)&&x>=0,(m.nutrient_diffusion,m.signal_diffusion,m.signal_decay,
        m.interface_barrier,m.membrane_barrier,m.signal_opening,m.remodeling_gain,
        m.forcing_relaxation,m.membrane_pinning,m.signal_bias)) || throw(ArgumentError("medium parameters must be finite and nonnegative"))
    all(a->size(a)==w.grid.size,(m.membrane,m.permeability,m.signal,m.signal_source,
        m.remodeling,m.forcing,m.scratch,m.phase.phi,m.phase.interface,m.phase.chemical_forcing)) || throw(DimensionMismatch("medium/world grid mismatch"))
    m.membrane === w.fields.deposit || throw(ArgumentError("medium membrane must alias world deposit"))
    all(x->isfinite(x)&&x>=0,m.membrane) || throw(ArgumentError("membrane material must be finite and nonnegative"))
    all(x->isfinite(x)&&x>=0,m.signal_source) || throw(ArgumentError("signal sources must be finite and nonnegative"))
    all(x->isfinite(x)&&x>=0,m.signal) || throw(ArgumentError("signals must be finite and nonnegative"))
    all(x->isfinite(x)&&x>=0,w.fields.nutrient) || throw(ArgumentError("nutrients must be finite and nonnegative"))
    all(isfinite,m.forcing) || throw(ArgumentError("forcing must be finite"))
    all(isfinite,m.phase.phi) && all(x->isfinite(x)&&x>=0,m.phase.interface) || throw(ArgumentError("phase environment must be finite with nonnegative interface"))
    all(isfinite,m.remodeling) || throw(ArgumentError("remodeling cues must be finite"))
    nothing
end
"""Refresh transport permeability from current interfaces, funded deposits and signals."""
function medium_observe!(m::ActiveMedium)
    @. m.permeability=1f0/(1f0+(m.interface_barrier*m.phase.interface+
        m.membrane_barrier*m.membrane)/(1f0+m.signal_opening*m.signal))
    m
end
# Each directed face has the same coefficient at either end. Updating into a
# scratch field makes exchanges equal and opposite, including periodic faces.
function medium_diffuse!(out,field,permeability,grid,D,dt)
    backend=KernelAbstractions.get_backend(field)
    if !(backend isa CPU)
        size(out)==size(field)==size(permeability)==grid.size || throw(DimensionMismatch("medium diffusion grid"))
        KernelAbstractions.get_backend(out)==backend==KernelAbstractions.get_backend(permeability) || throw(ArgumentError("medium diffusion backends differ"))
        (Base.mightalias(out,field)||Base.mightalias(out,permeability)) && throw(ArgumentError("device medium diffusion requires scratch distinct from both inputs"))
        medium_diffuse_kernel!(backend,256)(out,field,permeability,grid.spacing,D,dt;ndrange=grid.size)
        copyto!(field,out)
        KernelAbstractions.synchronize(backend)
        return field
    end
    dims=grid.size;h=grid.spacing
    for k in 1:dims[3],j in 1:dims[2],i in 1:dims[1]
        out[i,j,k]=_medium_diffused_value(field,permeability,(i,j,k),dims,h,D,dt)
    end
    copyto!(field,out)
end

@inline function _medium_face_flux(center,neighbor,a,b,D,h2)
    face=(a+b)==0f0 ? 0f0 : 2f0*a*b/(a+b)
    D*face*(neighbor-center)/h2
end

@inline function _medium_diffused_value(field,permeability,pos,dims,h,D,dt)
    flux=0f0
    for axis in 1:3
        if dims[axis]>1
            for side in (-1,1)
                q=phase_neighbor(pos,axis,side,dims)
                flux+=_medium_face_flux(field[pos...],field[q...],permeability[pos...],permeability[q...],D,h[axis]^2)
            end
        end
    end
    field[pos...]+dt*flux
end
@kernel function medium_diffuse_kernel!(out,field,permeability,h,D,dt)
    i,j,k=@index(Global,NTuple)
    out[i,j,k]=_medium_diffused_value(field,permeability,(i,j,k),size(field),h,D,dt)
end

# Concrete CPU fields use direct periodic neighbors. The order of the six face
# contributions is identical to the generic axis/side loop above.
function medium_diffuse!(out::Array{Float32,3},field::Array{Float32,3},
        permeability::Array{Float32,3},grid::Grid,D::Float32,dt::Float32)
    dims=grid.size
    size(out)==size(field)==size(permeability)==dims || throw(DimensionMismatch("medium diffusion grid"))
    nx,ny,nz=dims;h2=grid.spacing.^2
    @inbounds for k in 1:nz,j in 1:ny,i in 1:nx
        center=field[i,j,k];a=permeability[i,j,k];flux=0f0
        if nx>1
            im=i==1 ? nx : i-1;ip=i==nx ? 1 : i+1
            flux+=_medium_face_flux(center,field[im,j,k],a,permeability[im,j,k],D,h2[1])
            flux+=_medium_face_flux(center,field[ip,j,k],a,permeability[ip,j,k],D,h2[1])
        end
        if ny>1
            jm=j==1 ? ny : j-1;jp=j==ny ? 1 : j+1
            flux+=_medium_face_flux(center,field[i,jm,k],a,permeability[i,jm,k],D,h2[2])
            flux+=_medium_face_flux(center,field[i,jp,k],a,permeability[i,jp,k],D,h2[2])
        end
        if nz>1
            km=k==1 ? nz : k-1;kp=k==nz ? 1 : k+1
            flux+=_medium_face_flux(center,field[i,j,km],a,permeability[i,j,km],D,h2[3])
            flux+=_medium_face_flux(center,field[i,j,kp],a,permeability[i,j,kp],D,h2[3])
        end
        out[i,j,k]=center+dt*flux
    end
    copyto!(field,out)
end
function advance!(w::World,m::ActiveMedium,dt::Real)
    check_dt(w,m,dt)
    @. m.signal+=m.signal_source
    fill!(m.signal_source,0f0)
    # Cues accumulate into a decaying material-independent chemical bias.
    @. m.forcing+=m.remodeling_gain*m.remodeling
    fill!(m.remodeling,0f0)
    medium_observe!(m)
    D=max(m.nutrient_diffusion,m.signal_diffusion)
    invh2=sum((inv(Float64(w.grid.spacing[d])^2) for d in 1:3 if w.grid.size[d]>1);init=0.0)
    steps=max(1,ceil(Int,Float64(dt)*2D*invh2/.9))
    subdt=Float32(dt/steps)
    for _ in 1:steps
        medium_diffuse!(m.scratch,w.fields.nutrient,m.permeability,w.grid,m.nutrient_diffusion,subdt)
        medium_diffuse!(m.scratch,m.signal,m.permeability,w.grid,m.signal_diffusion,subdt)
    end
    m.signal .*= exp(-m.signal_decay*Float32(dt))
    m.forcing .*= exp(-m.forcing_relaxation*Float32(dt))
    medium_observe!(m)
    # Chemical-potential forcing acts under the phase solver's Laplacian:
    # boundaries deform without creating or removing phase composition.
    @. m.phase.chemical_forcing=m.forcing+m.membrane_pinning*m.membrane*m.phase.phi-
        m.signal_bias*m.signal
    w
end

"""Record medium state; `membrane` is an alias of world deposit. Record the
referenced phase as its own process stage to capture phase state and rules."""
function record_process_fields(m::ActiveMedium)
    names=Tuple(n for n in fieldnames(typeof(m)) if n!=:phase)
    merge(NamedTuple{names}(Tuple(getfield(m,n) for n in names)),
        (membrane_alias="world.fields.deposit",phase_reference="separate phase process stage"))
end
