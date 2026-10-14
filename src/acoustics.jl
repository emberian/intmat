"""Damped scalar wave propagation with mineral-dependent face coefficients.
The pressure-like scalar obeys p_tt + 2γ p_t = div(c(d)^2 grad p) + source.
This is a nondimensional artistic acoustic analogue, not a calibrated solid model.
"""
mutable struct AcousticField{A}
    pressure::A
    previous::A
    scratch::A
    envelope::A
    speed::Float32
    material_feedback::Float32
    damping::Float32
    source_amplitude::Float32
    frequency::Float32
    source_radius::Float32
    source_center::NTuple{3,Float32}
    activity_gain::Float32
    intensity_scale::Float32
    envelope_time::Float32
    timestep::Float32
end
function AcousticField(w::World; speed=0.3f0,material_feedback=8f0,damping=0.08f0,
                       source_amplitude=18f0,frequency=0.8f0,source_radius=0.08f0,
                       source_center=(-0.45f0,0f0,0f0),activity_gain=0.35f0,
                       intensity_scale=0.001f0,envelope_time=0.4f0)
    alloc()=KernelAbstractions.zeros(w.backend,Float32,w.grid.size...)
    p=AcousticField(alloc(),alloc(),alloc(),alloc(),Float32(speed),Float32(material_feedback),Float32(damping),
        Float32(source_amplitude),Float32(frequency),Float32(source_radius),Float32.(Tuple(source_center)),
        Float32(activity_gain),Float32(intensity_scale),Float32(envelope_time),0f0)
    validate(w,p); p
end
function validate(w::World,p::AcousticField)
    vals=(p.speed,p.material_feedback,p.damping,p.source_amplitude,p.frequency,p.source_radius,p.activity_gain,p.intensity_scale,p.envelope_time)
    all(isfinite,vals) && all(isfinite,p.source_center) || throw(ArgumentError("acoustic parameters must be finite"))
    p.speed>0 && p.source_radius>0 && p.intensity_scale>0 && p.envelope_time>0 || throw(ArgumentError("acoustic speed and scales must be positive"))
    all(>=(0),(p.material_feedback,p.damping,p.frequency,p.activity_gain)) || throw(ArgumentError("acoustic coupling and damping must be nonnegative"))
    all(a->size(a)==w.grid.size,(p.pressure,p.previous,p.scratch,p.envelope)) || throw(ArgumentError("acoustic field grid mismatch"))
    all(a->typeof(KernelAbstractions.get_backend(a))==typeof(w.backend),(p.pressure,p.previous,p.scratch,p.envelope)) || throw(ArgumentError("acoustic backend mismatch"))
    minimum(w.fields.deposit)>=0 || throw(ArgumentError("acoustic mineral must be nonnegative"))
end
stability_limit(w,p::AcousticField)=0.9f0/(p.speed*sqrt(sum(inv(w.grid.spacing[d]^2) for d in 1:(w.grid.size[3]==1 ? 2 : 3))))
@kernel function acoustic_kernel!(out,pressure,previous,envelope,activity,deposit,h,extent,t,dt,
                                 speed,feedback,damping,amplitude,frequency,radius,center,gain,scale,tau)
    i,j,k=@index(Global,NTuple)
    dims=size(pressure); u=pressure[i,j,k]; divergence=0f0
    for axis in 1:3
        if dims[axis]>1
            for side in (-1,1)
                index=(i,j,k)[axis]+side
                index=index<1 ? dims[axis] : (index>dims[axis] ? 1 : index)
                neighbor=axis==1 ? (index,j,k) : (axis==2 ? (i,index,k) : (i,j,index))
                c2=speed^2/(1f0+feedback*(deposit[i,j,k]+deposit[neighbor...])/2f0)
                divergence+=c2*(pressure[neighbor...]-u)/h[axis]^2
            end
        end
    end
    x=(i-0.5f0)*h[1]-extent[1]/2-center[1]
    y=(j-0.5f0)*h[2]-extent[2]/2-center[2]
    z=dims[3]==1 ? 0f0 : (k-0.5f0)*h[3]-extent[3]/2-center[3]
    source=amplitude*exp(-(x*x+y*y+z*z)/(2radius^2))*sin(2f0*Float32(pi)*frequency*t)
    next=(2u-(1f0-damping*dt)*previous[i,j,k]+dt^2*(divergence+source))/(1f0+damping*dt)
    out[i,j,k]=next
    intensity=envelope[i,j,k]+(1f0-exp(-dt/tau))*(next^2-envelope[i,j,k])
    envelope[i,j,k]=intensity
    activity[i,j,k]+=dt*gain*intensity/(scale+intensity)
end
function advance!(w::World,p::AcousticField,dt::Real)
    check_dt(w,p,dt)
    p.timestep==0 || Float32(dt)==p.timestep || throw(ArgumentError("AcousticField requires a fixed timestep; rebuild state to change dt"))
    p.timestep=Float32(dt)
    launch!(w,acoustic_kernel!,p.scratch,p.pressure,p.previous,p.envelope,w.fields.activity,w.fields.deposit,
        w.grid.spacing,w.grid.extent,Float32(p.frequency>0 ? mod(w.time,1.0/Float64(p.frequency)) : 0.0),Float32(dt),p.speed,p.material_feedback,p.damping,p.source_amplitude,
        p.frequency,p.source_radius,p.source_center,p.activity_gain,p.intensity_scale,p.envelope_time)
    p.previous,p.pressure,p.scratch=p.pressure,p.scratch,p.previous
    copyto!(w.fields.pressure,p.pressure)
    w
end
"""Initialize a pulse at zero initial velocity, respecting the leapfrog half-acceleration."""
function seed_acoustic!(w::World,p::AcousticField; amplitude=1f0,radius=0.12f0,center=(0f0,0f0,0f0),dt=0.01f0)
    validate(w,p); check_dt(w,p,dt)
    radius>0 && isfinite(radius) && isfinite(amplitude) || throw(ArgumentError("invalid acoustic pulse"))
    host=zeros(Float32,w.grid.size)
    for k in axes(host,3),j in axes(host,2),i in axes(host,1)
        distance=sum((((i,j,k)[d]-0.5f0)*w.grid.spacing[d]-w.grid.extent[d]/2-center[d])^2 for d in 1:(w.grid.size[3]==1 ? 2 : 3))
        host[i,j,k]=amplitude*exp(-distance/(2radius^2))
    end
    prev=copy(host); deposit=Array(w.fields.deposit)
    for k in axes(host,3),j in axes(host,2),i in axes(host,1)
        divergence=0f0
        for axis in 1:3
            w.grid.size[axis]==1 && continue
            for side in (-1,1)
                q=mod1((i,j,k)[axis]+side,w.grid.size[axis])
                neighbor=axis==1 ? (q,j,k) : axis==2 ? (i,q,k) : (i,j,q)
                c2=p.speed^2/(1f0+p.material_feedback*(deposit[i,j,k]+deposit[neighbor...])/2f0)
                divergence+=c2*(host[neighbor...]-host[i,j,k])/w.grid.spacing[axis]^2
            end
        end
        prev[i,j,k]+=0.5f0*Float32(dt)^2*divergence
    end
    copyto!(p.pressure,host); copyto!(p.previous,prev)
    copyto!(w.fields.pressure,host)
    p.timestep=Float32(dt)
    p
end

function check_dt(w,p::AcousticField,dt)
    invoke(check_dt,Tuple{Any,Any,Any},w,p,dt)
    p.timestep==0 || Float32(dt)==p.timestep || throw(ArgumentError("AcousticField requires a fixed timestep"))
end
