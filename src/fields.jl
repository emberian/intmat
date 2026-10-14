"""FitzHugh–Nagumo excitable field; accumulated deposit suppresses transport and shifts threshold."""
Base.@kwdef struct Excitable
    diffusion::Float32=0.008f0
    threshold::Float32=0.12f0
    recovery_rate::Float32=0.08f0
    recovery_coupling::Float32=0.8f0
    gain::Float32=3f0
    material_feedback::Float32=0.5f0
end
function stability_limit(w,p::Excitable)
    # Bound the local reaction Jacobian over current activity and every admissible
    # mineral-shifted threshold. This is a step-size guard, not a global guarantee
    # for the unbounded cubic after arbitrary preceding process updates.
    ulo,uhi=extrema(w.fields.activity)
    isfinite(ulo) && isfinite(uhi) || return 0f0
    derivative=0.0
    for theta in (Float64(p.threshold),Float64(p.threshold+p.material_feedback))
        vertex=clamp((1+theta)/3,Float64(ulo),Float64(uhi))
        for u in (Float64(ulo),Float64(uhi),vertex)
            derivative=max(derivative,abs(Float64(p.gain)*(-3u*u+2*(1+theta)*u-theta)))
        end
    end
    reaction=max(derivative+1,Float64(p.recovery_rate)*(1+Float64(p.recovery_coupling)))
    axes=w.grid.size[3]==1 ? (1:2) : (1:3)
    diffusion=2Float64(p.diffusion)*sum(inv(Float64(w.grid.spacing[d])^2) for d in axes)
    min(0.05f0,diffusion_limit(w,p.diffusion),0.8/(reaction+diffusion))
end
function validate(w::World,p::Excitable)
    validate(p)
    lo,hi=extrema(w.fields.activity); rlo,rhi=extrema(w.fields.recovery)
    dlo,dhi=extrema(w.fields.deposit)
    all(isfinite,(lo,hi,rlo,rhi,dlo,dhi)) && dlo>=0 ||
        throw(ArgumentError("excitable state must be finite with nonnegative deposit"))
end
@kernel function excite_kernel!(ao,ro,a,r,d,h,dt,p)
    i,j,k=@index(Global,NTuple)
    u=a[i,j,k]; v=r[i,j,k]; mineral=d[i,j,k]
    threshold=p.threshold+p.material_feedback*mineral/(1f0+mineral)
    ao[i,j,k]=u+dt*(p.diffusion/(1f0+mineral)*laplacian(a,i,j,k,h)+p.gain*u*(1f0-u)*(u-threshold)-v)
    ro[i,j,k]=v+dt*p.recovery_rate*(u-p.recovery_coupling*v)
end
function advance!(w::World,p::Excitable,dt::Real)
    check_dt(w,p,dt)
    f=w.fields; s=w.scratch
    launch!(w,excite_kernel!,s.activity,s.recovery,f.activity,f.recovery,f.deposit,w.grid.spacing,Float32(dt),p)
    copyto!(f.activity,s.activity); copyto!(f.recovery,s.recovery)
    w
end
"""Conservative nutrient-to-solid reaction, gated by activity. Dissolution returns feedstock."""
Base.@kwdef struct Deposition
    rate::Float32=0.6f0
    threshold::Float32=0.25f0
    dissolution::Float32=0.002f0
end
@kernel function deposit_kernel!(n,d,a,dt,p)
    i,j,k=@index(Global,NTuple)
    solid=d[i,j,k]; food=n[i,j,k]
    grow=min(food,max(0f0,p.rate*dt*max(0f0,a[i,j,k]-p.threshold)*food))
    dissolve=min(solid,max(0f0,p.dissolution*dt*solid))
    flux=grow-dissolve
    n[i,j,k]=food-flux; d[i,j,k]=solid+flux
end
function advance!(w::World,p::Deposition,dt::Real)
    check_dt(w,p,dt)
    launch!(w,deposit_kernel!,w.fields.nutrient,w.fields.deposit,w.fields.activity,Float32(dt),p)
    w
end
"""Conservative periodic diffusion of a selected field."""
struct Diffusion
    field::Symbol
    rate::Float32
    material_feedback::Float32
    function Diffusion(field::Symbol,rate::Real; material_feedback=0f0)
        isfinite(rate)&&rate>=0 || throw(ArgumentError("diffusion rate must be finite and nonnegative"))
        isfinite(material_feedback)&&material_feedback>=0 || throw(ArgumentError("material feedback must be finite and nonnegative"))
        new(field,Float32(rate),Float32(material_feedback))
    end
end
stability_limit(w,p::Diffusion)=diffusion_limit(w,p.rate)
function validate(w::World,p::Diffusion)
    hasproperty(w.fields,p.field) && hasproperty(w.scratch,p.field) ||
        throw(ArgumentError("unknown diffusion field $(p.field)"))
    isfinite(p.rate) && p.rate>=0 && isfinite(p.material_feedback) && p.material_feedback>=0 ||
        throw(ArgumentError("diffusion parameters must be finite and nonnegative"))
end
@kernel function diffuse_kernel!(out,a,deposit,h,dt,D,feedback)
    i,j,k=@index(Global,NTuple)
    v=0f0
    dims=size(a)
    for axis in 1:3
        if dims[axis]>1
            for side in (-1,1)
                pos=(i,j,k)
                q=pos[axis]+side
                q=q<1 ? dims[axis] : (q>dims[axis] ? 1 : q)
                neighbor=axis==1 ? (q,j,k) : (axis==2 ? (i,q,k) : (i,j,q))
                mobility=1f0/(1f0+feedback*(deposit[i,j,k]+deposit[neighbor...])/2f0)
                v+=mobility*(a[neighbor...]-a[i,j,k])/h[axis]^2
            end
        end
    end
    out[i,j,k]=a[i,j,k]+dt*D*v
end
function advance!(w::World,p::Diffusion,dt::Real)
    check_dt(w,p,dt)
    a=getproperty(w.fields,p.field); out=getproperty(w.scratch,p.field)
    launch!(w,diffuse_kernel!,out,a,w.fields.deposit,w.grid.spacing,Float32(dt),p.rate,p.material_feedback)
    copyto!(a,out)
    w
end
"""External standing pressure wave and its velocity; an explicitly driven artistic medium."""
Base.@kwdef struct WaveForcing
    amplitude::Float32=0.12f0
    frequency::Float32=0.7f0
    wavenumber::Float32=3f0
    activity_gain::Float32=0.08f0
end
@kernel function wave_kernel!(pressure,vx,vy,vz,a,d,h,extent,t,dt,p)
    i,j,k=@index(Global,NTuple)
    x=(i-0.5f0)*h[1]-extent[1]/2; y=(j-0.5f0)*h[2]-extent[2]/2
    phase=2f0*Float32(pi)*p.frequency*t
    q=p.wavenumber
    amp=p.amplitude/(1f0+d[i,j,k])
    pressure[i,j,k]=amp*sin(q*x)*sin(q*y)*cos(phase)
    vx[i,j,k]=-amp*q*cos(q*x)*sin(q*y)*sin(phase)
    vy[i,j,k]=amp*q*sin(q*x)*cos(q*y)*sin(phase)
    vz[i,j,k]=0f0
    a[i,j,k]+=dt*p.activity_gain*pressure[i,j,k]
end
function advance!(w::World,p::WaveForcing,dt::Real)
    check_dt(w,p,dt)
    f=w.fields
    launch!(w,wave_kernel!,f.pressure,f.vx,f.vy,f.vz,f.activity,f.deposit,w.grid.spacing,w.grid.extent,Float32(p.frequency!=0 ? mod(w.time,1.0/abs(Float64(p.frequency))) : 0.0),Float32(dt),p)
    w
end

function validate(p::Excitable)
    all(isfinite,(p.diffusion,p.threshold,p.recovery_rate,p.recovery_coupling,p.gain,p.material_feedback)) || throw(ArgumentError("excitable parameters must be finite"))
    all(>=(0),(p.diffusion,p.recovery_rate,p.recovery_coupling,p.gain,p.material_feedback)) || throw(ArgumentError("excitable rates must be nonnegative"))
end
function validate(p::Deposition)
    all(isfinite,(p.rate,p.threshold,p.dissolution)) && p.rate>=0 && p.dissolution>=0 || throw(ArgumentError("deposition rates must be finite and nonnegative"))
end
function validate(p::WaveForcing)
    all(isfinite,(p.amplitude,p.frequency,p.wavenumber,p.activity_gain)) || throw(ArgumentError("wave parameters must be finite"))
end
