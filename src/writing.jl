"""Focused thermal writing with persistent temperature and irreversible dose.

A deliberate parametric toolpath deposits a Gaussian heat source. Temperature
then diffuses and cools; exposure above `threshold` converts available nutrient
to deposit conservatively. This is a nondimensional phenomenological writing
model, not optical breakdown, pulse-scale laser physics, or sonoluminescence.
`pulse_rate` controls a smooth intensity envelope (zero means continuous).
"""
struct ThermalWriter{A}
    temperature::A
    scratch::A
    dose::A
    path::Symbol
    center::NTuple{3,Float32}
    radius::Float32
    height::Float32
    period::Float32
    width::Float32
    heat::Float32
    diffusion::Float32
    cooling::Float32
    threshold::Float32
    conversion::Float32
    pulse_rate::Float32
    start::Float32
    stop::Float32
end
function ThermalWriter(w::World;path=:lissajous,center=(0,0,0),radius=0.55,height=0.4,
                       period=4,width=0.045,heat=16,diffusion=0.0003,cooling=4,
                       threshold=0.1,conversion=12,pulse_rate=0,start=0,stop=4)
    arrays=ntuple(_->KernelAbstractions.zeros(w.backend,Float32,w.grid.size...),3)
    p=ThermalWriter(arrays...,Symbol(path),Tuple(Float32.(center)),Float32(radius),Float32(height),Float32(period),Float32(width),Float32(heat),Float32(diffusion),Float32(cooling),Float32(threshold),Float32(conversion),Float32(pulse_rate),Float32(start),Float32(stop))
    validate(p);p
end
function validate(p::ThermalWriter)
    p.path in (:stationary,:helix,:lissajous) || throw(ArgumentError("writer path must be stationary, helix, or lissajous"))
    all(isfinite,(p.center...,p.radius,p.height,p.period,p.width,p.heat,p.diffusion,p.cooling,p.threshold,p.conversion,p.pulse_rate,p.start,p.stop)) || throw(ArgumentError("writer parameters must be finite"))
    p.period>0 && p.width>0 && p.stop>=p.start || throw(ArgumentError("writer period and width must be positive, stop must follow start"))
    all(>=(0),(p.radius,p.height,p.heat,p.diffusion,p.cooling,p.threshold,p.conversion,p.pulse_rate)) || throw(ArgumentError("writer rates and dimensions must be nonnegative"))
end
function validate(w::World,p::ThermalWriter)
    validate(p)
    for a in (p.temperature,p.scratch,p.dose)
        size(a)==w.grid.size || throw(ArgumentError("writer grid differs from world"))
        KernelAbstractions.get_backend(a)==w.backend || throw(ArgumentError("writer backend differs from world"))
    end
end
"""Beam center in physical grid coordinates at time t; not a solver trajectory."""
function beam_position(p::ThermalWriter,t::Real)
    phase=2f0*Float32(pi)*Float32(t-p.start)/p.period
    offset=p.path==:stationary ? (0f0,0f0,0f0) : p.path==:helix ?
        (p.radius*cos(3phase),p.radius*sin(3phase),p.height*(2f0*Float32(t-p.start)/p.period-1f0)) :
        (p.radius*sin(2phase),p.radius*sin(3phase+0.4f0),p.height*cos(phase))
    ntuple(d->p.center[d]+offset[d],3)
end
function stability_limit(w,p::ThermalWriter)
    rate=p.cooling+2f0*p.diffusion*sum(inv(w.grid.spacing[d]^2) for d in 1:(w.grid.size[3]==1 ? 2 : 3))
    rate==0 ? Inf32 : 0.9f0/rate
end
@kernel function write_heat_kernel!(out,temp,dose,nutrient,deposit,h,extent,dt,position,width,power,D,cooling,threshold,conversion)
    i,j,k=@index(Global,NTuple)
    x=(i-.5f0)*h[1]-extent[1]/2-position[1]
    y=(j-.5f0)*h[2]-extent[2]/2-position[2]
    z=size(temp,3)==1 ? 0f0 : (k-.5f0)*h[3]-extent[3]/2-position[3]
    source=power*exp(-(x*x+y*y+z*z)/(2f0*width*width))
    T=max(0f0,temp[i,j,k]+dt*(D*laplacian(temp,i,j,k,h)-cooling*temp[i,j,k]+source))
    exposure=dt*max(0f0,T-threshold)
    dose[i,j,k]+=exposure
    food=nutrient[i,j,k]
    flux=food*(-expm1(-conversion*exposure))
    nutrient[i,j,k]=food-flux; deposit[i,j,k]+=flux
    out[i,j,k]=T
end
function advance!(w::World,p::ThermalWriter,dt::Real)
    check_dt(w,p,dt)
    t=w.time+Float32(dt)/2
    active=p.start<=t<p.stop
    envelope=p.pulse_rate==0 ? 1f0 : .5f0*(1f0+cos(2f0*Float32(pi)*p.pulse_rate*t))
    power=active ? p.heat*envelope : 0f0
    launch!(w,write_heat_kernel!,p.scratch,p.temperature,p.dose,w.fields.nutrient,w.fields.deposit,w.grid.spacing,w.grid.extent,Float32(dt),beam_position(p,t),p.width,power,p.diffusion,p.cooling,p.threshold,p.conversion)
    copyto!(p.temperature,p.scratch)
    w
end

function record_process_fields(p::ThermalWriter)
    names=fieldnames(typeof(p))
    NamedTuple{names}(Tuple(getfield(p,n) for n in names))
end
