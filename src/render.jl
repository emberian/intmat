"""Rendering controls in nondimensional specimen coordinates; orbit is in radians."""
Base.@kwdef struct RenderSettings
    width::Int = 640
    height::Int = 640
    steps::Int = 160
    orbit::Float64 = 0.55
    elevation::Float64 = 0.23
    absorption::Float64 = 1.1
    emission::Float64 = 3.5
    scattering::Float64 = 0.18
    exposure::Float64 = 1.4
    cutaway::Float64 = 1.0
end

@inline _vadd(a,b) = ntuple(i->a[i]+b[i],3)
@inline _vmul(a,s) = ntuple(i->a[i]*s,3)
@inline _dot(a,b) = sum(a[i]*b[i] for i in 1:3)
@inline _cross(a,b) = (a[2]*b[3]-a[3]*b[2],a[3]*b[1]-a[1]*b[3],a[1]*b[2]-a[2]*b[1])
@inline _unit(a) = _vmul(a,1/sqrt(_dot(a,a)))

"""Trilinear sampling of a cell-centered field on [-1,1]^3; zero outside."""
@inline function volume_sample(a::AbstractArray{<:Real,3}, p)
    any(x->abs(x)>1,p) && return 0.0
    n=size(a)
    q=ntuple(d->clamp((p[d]+1)*n[d]/2+0.5,1.0,Float64(n[d])),3)
    lo=ntuple(d->floor(Int,q[d]),3)
    hi=ntuple(d->min(lo[d]+1,n[d]),3)
    w=ntuple(d->q[d]-lo[d],3)
    value=0.0
    @inbounds for k in 0:1, j in 0:1, i in 0:1
        value += a[i==0 ? lo[1] : hi[1],j==0 ? lo[2] : hi[2],k==0 ? lo[3] : hi[3]] *
            (i==0 ? 1-w[1] : w[1])*(j==0 ? 1-w[2] : w[2])*(k==0 ? 1-w[3] : w[3])
    end
    value
end

"""Exact constant-cell Beer–Lambert update, with source per unit length."""
@inline function transport_segment(transmittance, radiance, extinction, source, ds)
    atten=exp(-extinction*ds)
    weight=extinction>1e-12 ? -expm1(-extinction*ds)/extinction : ds
    (transmittance*atten, radiance+transmittance*source*weight)
end

function _box_interval(origin, direction)
    near=-Inf; far=Inf
    for d in 1:3
        if abs(direction[d])<1e-12
            abs(origin[d])>1 && return nothing
        else
            a=(-1-origin[d])/direction[d]; b=(1-origin[d])/direction[d]
            near=max(near,min(a,b)); far=min(far,max(a,b))
        end
    end
    far>max(near,0) ? (max(near,0),far) : nothing
end

"""Render real 3D fields to linear RGB. Required: density; optional mask, activity,
structure. Positive activity emits cyan light; structure increases amber absorption.
A geometric mask confines the carrier. Single scattering uses a bounded directional
light with twelve-step shadow integration. This is straight-ray transport, not refraction.
"""
function render_volume(fields::NamedTuple; settings=RenderSettings())
    s=settings
    s.width>0 && s.height>0 && s.steps>0 || throw(ArgumentError("positive render dimensions and steps required"))
    density=fields.density
    ndims(density)==3 || throw(ArgumentError("density must be a 3D field"))
    for name in (:mask,:activity,:structure)
        hasproperty(fields,name) && size(getproperty(fields,name))!=size(density) && throw(DimensionMismatch("$name shape differs"))
    end
    eye=(3.8*cos(s.elevation)*sin(s.orbit),3.8*sin(s.elevation),3.8*cos(s.elevation)*cos(s.orbit))
    forward=_unit(_vmul(eye,-1)); right=_unit(_cross(forward,(0.,1.,0.))); up=_cross(right,forward)
    light=_unit((-0.6,0.8,0.5)); rgb=zeros(Float32,3,s.width,s.height)
    sample(name,p,default)=hasproperty(fields,name) ? volume_sample(getproperty(fields,name),p) : default
    Threads.@threads for y in 1:s.height
        for x in 1:s.width
            u=(2*(x-0.5)/s.width-1)*s.width/s.height*0.39
            v=(1-2*(y-0.5)/s.height)*0.39
            ray=_unit(_vadd(forward,_vadd(_vmul(right,u),_vmul(up,v))))
            interval=_box_interval(eye,ray)
            background=(0.025,0.035,0.055)
            if interval===nothing
                rgb[:,x,y].=background; continue
            end
            t0,t1=interval; ds=(t1-t0)/s.steps
            T=(1.,1.,1.); L=(0.,0.,0.)
            for step in 1:s.steps
                p=_vadd(eye,_vmul(ray,t0+(step-0.5)*ds))
                p[3]>s.cutaway && continue
                m=clamp(sample(:mask,p,1.),0.,1.)
                rho=max(0.,volume_sample(density,p))*m
                rho<1e-7 && continue
                activity=max(0.,sample(:activity,p,0.))*m
                structure=max(0.,sample(:structure,p,0.))*m
                sigma=ntuple(c->s.absorption*rho*(0.26,0.72,1.8)[c]+structure*(0.16,0.42,0.85)[c]+s.scattering*rho,3)
                illumination=1.0
                if s.scattering>0
                    opticaldepth=0.0
                    for shadow in 1:12
                        q=_vadd(p,_vmul(light,(shadow-0.5)*0.16))
                        q[3]>s.cutaway && continue
                        opticaldepth+=max(0.,volume_sample(density,q))*clamp(sample(:mask,q,1.),0.,1.)*0.16
                    end
                    illumination=exp(-(s.absorption*0.65+s.scattering)*opticaldepth)
                end
                source=ntuple(c->s.emission*activity*(0.06,0.8,1.0)[c]+s.scattering*rho*illumination*(2.2,1.4,0.65)[c],3)
                result=ntuple(c->transport_segment(T[c],L[c],sigma[c],source[c],ds),3)
                T=ntuple(c->result[c][1],3); L=ntuple(c->result[c][2],3)
                maximum(T)<0.001 && break
            end
            @inbounds for c in 1:3
                rgb[c,x,y]=L[c]+T[c]*background[c]
            end
        end
    end
    rgb
end

"""Write a binary PPM with exponential exposure mapping and sRGB encoding."""
function write_ppm(path::AbstractString, rgb; exposure=1.4)
    size(rgb,1)==3 || throw(DimensionMismatch("RGB channel dimension must be 3"))
    mkpath(dirname(abspath(path)))
    open(path,"w") do io
        write(io,"P6\n$(size(rgb,2)) $(size(rgb,3))\n255\n")
        for y in axes(rgb,3), x in axes(rgb,2), c in 1:3
            linear=1-exp(-max(0.,rgb[c,x,y])*exposure)
            srgb=linear<=0.0031308 ? 12.92linear : 1.055linear^(1/2.4)-0.055
            write(io,UInt8(round(Int,255clamp(srgb,0.,1.))))
        end
    end
    path
end
render_volume(fields::NamedTuple,path::AbstractString;settings=RenderSettings()) = write_ppm(path,render_volume(fields;settings);exposure=settings.exposure)

"""Write Mitsuba VOL v3, float32 scalar grid, x axis fastest, bounds [-1,1]^3.
The format is little endian, as supported on Julia's standard desktop platforms."""
function write_vol(path::AbstractString, data::AbstractArray{<:Real,3})
    Base.ENDIAN_BOM==0x04030201 || error("VOL writer currently requires a little-endian host")
    mkpath(dirname(abspath(path)))
    open(path,"w") do io
        write(io,codeunits("VOL")); write(io,UInt8(3)); write(io,Int32(1))
        write(io,Int32[size(data)...]); write(io,Int32(1))
        write(io,Float32[-1,-1,-1,1,1,1]); write(io,Float32.(vec(data)))
    end
    path
end

"""Export each scalar field to Mitsuba VOL and float32 raw, with a JSON manifest.
The raw layout is Julia column-major (x fastest), explicitly recorded for conversion."""
function export_volume(fields::NamedTuple, directory::AbstractString)
    mkpath(directory)
    names=Symbol[]
    shape=nothing
    for name in keys(fields)
        data=getproperty(fields,name)
        data isa AbstractArray{<:Real,3} || continue
        if shape===nothing
            shape=size(data)
        elseif size(data)!=shape
            throw(DimensionMismatch("all export fields must share one shape"))
        end
        push!(names,name)
        write_vol(joinpath(directory,"$name.vol"),data)
        open(joinpath(directory,"$name.f32"),"w") do io
            write(io,Float32.(vec(data)))
        end
    end
    isempty(names) && throw(ArgumentError("no scalar volume fields"))
    dims=size(getproperty(fields,first(names)))
    open(joinpath(directory,"volume.json"),"w") do io
        print(io,"{\"shape\": [$(join(dims,","))], \"bounds\": [-1,-1,-1,1,1,1], \"dtype\": \"little-endian float32\", \"order\": \"F\", \"fields\": [")
        print(io,join(["\"$name\"" for name in names],",")); print(io,"]}\n")
    end
    joinpath(directory,"volume.json")
end
