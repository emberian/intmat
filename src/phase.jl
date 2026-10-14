"""Conserved periodic Cahn–Hilliard material with mineral pinning and interface feedback.
`phi` is a signed composition, not a clipped volume fraction. All persistent state
is process-owned. The constructor seed affects only the initial perturbation.
"""
struct PhaseSeparation{A}
    phi::A
    chemical_potential::A
    interface::A
    scratch::A
    bulk::Float32
    kappa::Float32
    mobility::Float32
    mineral_pinning::Float32
    nutrient_diffusion::Float32
    phase_barrier::Float32
    interface_activity::Float32
    activity_relaxation::Float32
end
function PhaseSeparation(w::World; mean=0f0, noise=0.12f0, seed=w.seed,
        bulk=1f0, kappa=0.7f0, mobility=1f0, mineral_pinning=2f0,
        nutrient_diffusion=0.04f0, phase_barrier=8f0,
        interface_activity=0.5f0, activity_relaxation=0.3f0)
    all(isfinite,(mean,noise)) && noise>=0 || throw(ArgumentError("finite mean and nonnegative noise required"))
    rng=MersenneTwister(seed)
    host=Float32(mean) .+ Float32(noise).*(rand(rng,Float32,w.grid.size...).-0.5f0)
    host .+= Float32(mean)-sum(host)/length(host)
    arrays=ntuple(_->KernelAbstractions.zeros(w.backend,Float32,w.grid.size...),4)
    copyto!(arrays[1],host)
    p=PhaseSeparation(arrays...,Float32.((bulk,kappa,mobility,mineral_pinning,nutrient_diffusion,phase_barrier,interface_activity,activity_relaxation))...)
    validate(w,p)
    phase_observe!(w,p)
    p
end
function validate(w::World,p::PhaseSeparation)
    parameters=(p.bulk,p.kappa,p.mobility,p.mineral_pinning,p.nutrient_diffusion,p.phase_barrier,p.interface_activity,p.activity_relaxation)
    all(x->isfinite(x)&&x>=0,parameters) || throw(ArgumentError("phase parameters must be finite and nonnegative"))
    p.bulk>0 && p.kappa>0 || throw(ArgumentError("phase bulk and kappa must be positive"))
    all(a->size(a)==w.grid.size,(p.phi,p.chemical_potential,p.interface,p.scratch)) || throw(DimensionMismatch("phase/world grid mismatch"))
    lo,hi=extrema(p.phi); dlo,dhi=extrema(w.fields.deposit)
    all(isfinite,(lo,hi,dlo,dhi)) && dlo>=0 || throw(ArgumentError("phase must be finite and deposit nonnegative"))
end
function stability_limit(w::World,p::PhaseSeparation)
    s=sum(inv(Float64(w.grid.spacing[d])^2) for d in 1:(w.grid.size[3]==1 ? 2 : 3))
    lo,hi=extrema(p.phi); r=max(1.5,abs(lo),abs(hi))
    # Conservative guard on the fourth-order explicit operator. Cubic bound uses
    # a margin beyond the current state; this is not unconditional stability.
    ch=p.mobility*4s*(p.bulk*(3r*r+1)+4p.kappa*s)
    coupling=2p.nutrient_diffusion*s+p.activity_relaxation
    0.7/max(ch,coupling,eps(Float64))
end
@inline function phase_neighbor(pos,axis,side,dims)
    q=pos[axis]+side
    q=q<1 ? dims[axis] : (q>dims[axis] ? 1 : q)
    axis==1 ? (q,pos[2],pos[3]) : axis==2 ? (pos[1],q,pos[3]) : (pos[1],pos[2],q)
end
@kernel function phase_observe_kernel!(mu,interface,phi,h,bulk,kappa)
    i,j,k=@index(Global,NTuple)
    u=phi[i,j,k]
    mu[i,j,k]=bulk*(u*u*u-u)-kappa*laplacian(phi,i,j,k,h)
    grad2=0f0; pos=(i,j,k); dims=size(phi)
    for axis in 1:3
        if dims[axis]>1
            plus=phase_neighbor(pos,axis,1,dims); minus=phase_neighbor(pos,axis,-1,dims)
            grad=(phi[plus...]-phi[minus...])/(2h[axis])
            grad2+=grad*grad
        end
    end
    interface[i,j,k]=sqrt(grad2)
end
function phase_observe!(w,p)
    launch!(w,phase_observe_kernel!,p.chemical_potential,p.interface,p.phi,w.grid.spacing,p.bulk,p.kappa)
end
@kernel function phase_flux_kernel!(out,phi,mu,deposit,h,dt,mobility,pinning)
    i,j,k=@index(Global,NTuple)
    pos=(i,j,k); dims=size(phi); flux=0f0
    for axis in 1:3
        if dims[axis]>1
            for side in (-1,1)
                q=phase_neighbor(pos,axis,side,dims)
                face=mobility/(1f0+pinning*(deposit[pos...]+deposit[q...])/2f0)
                flux+=face*(mu[q...]-mu[pos...])/h[axis]^2
            end
        end
    end
    out[pos...]=phi[pos...]+dt*flux
end
@kernel function phase_couple_kernel!(nout,n,a,phi,interface,h,dt,D,barrier,gain,relax)
    i,j,k=@index(Global,NTuple)
    pos=(i,j,k); dims=size(phi); flux=0f0
    for axis in 1:3
        if dims[axis]>1
            for side in (-1,1)
                q=phase_neighbor(pos,axis,side,dims)
                # Dense positive phase has smaller permeability. Face arithmetic
                # mean keeps equal and opposite nutrient transfers.
                left=1f0/(1f0+barrier*(1f0+tanh(phi[pos...]))/2f0)
                right=1f0/(1f0+barrier*(1f0+tanh(phi[q...]))/2f0)
                flux+=D*(left+right)/2f0*(n[q...]-n[pos...])/h[axis]^2
            end
        end
    end
    nout[pos...]=n[pos...]+dt*flux
    target=gain*interface[pos...]/(1f0+interface[pos...])
    a[pos...]+=dt*relax*(target-a[pos...])
end
function advance!(w::World,p::PhaseSeparation,dt::Real)
    check_dt(w,p,dt)
    phase_observe!(w,p)
    launch!(w,phase_flux_kernel!,p.scratch,p.phi,p.chemical_potential,w.fields.deposit,w.grid.spacing,Float32(dt),p.mobility,p.mineral_pinning)
    copyto!(p.phi,p.scratch)
    phase_observe!(w,p)
    launch!(w,phase_couple_kernel!,w.scratch.nutrient,w.fields.nutrient,w.fields.activity,p.phi,p.interface,w.grid.spacing,Float32(dt),p.nutrient_diffusion,p.phase_barrier,p.interface_activity,p.activity_relaxation)
    copyto!(w.fields.nutrient,w.scratch.nutrient)
    w
end
"""Discrete free energy, using the same forward-edge gradient as the chemical potential."""
function phase_energy(w::World,p::PhaseSeparation)
    phi=Array(p.phi); energy=0.0; dims=size(phi)
    for k in axes(phi,3),j in axes(phi,2),i in axes(phi,1)
        pos=(i,j,k); u=Float64(phi[pos...])
        energy+=p.bulk*(u*u-1)^2/4
        for axis in 1:3
            if dims[axis]>1
                q=phase_neighbor(pos,axis,1,dims)
                energy+=p.kappa/2*((Float64(phi[q...])-u)/w.grid.spacing[axis])^2
            end
        end
    end
    energy*prod(w.grid.spacing[d] for d in 1:(dims[3]==1 ? 2 : 3))
end
record_process_fields(p::PhaseSeparation)=NamedTuple{fieldnames(typeof(p))}(Tuple(getfield(p,n) for n in fieldnames(typeof(p))))
"""Optical observation from composition and real interfaces; it never mutates state."""
function phase_optical_fields(w::World,p::PhaseSeparation)
    phi=Array(p.phi); interface=Array(p.interface)
    fraction=(1f0 .+ tanh.(2f0.*phi))./2f0
    (density=0.04f0 .+ 0.7f0.*fraction, structure=fraction,
     activity=0.8f0.*interface .+ 0.2f0.*max.(Array(w.fields.activity),0f0))
end
