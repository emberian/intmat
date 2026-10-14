"""Periodic cell-centered Cartesian grid; nz=1 selects two spatial dimensions."""
struct Grid
    size::NTuple{3,Int}
    extent::NTuple{3,Float32}
    spacing::NTuple{3,Float32}
    function Grid(nx::Integer, ny::Integer, nz::Integer=1; extent=(2f0,2f0,2f0))
        dims = (Int(nx),Int(ny),Int(nz))
        all(>(0), dims) || throw(ArgumentError("grid dimensions must be positive"))
        ex = Float32.(Tuple(extent))
        length(ex)==3 && all(x->isfinite(x)&&x>0,ex) || throw(ArgumentError("extent must contain three positive finite lengths"))
        new(dims,ex,ex ./ dims)
    end
end
mutable struct World{B,F,S}
    grid::Grid
    backend::B
    fields::F
    scratch::S
    time::Float64
    seed::Int
end
function World(grid::Grid; seed=17, backend=CPU())
    names = (:nutrient,:deposit,:activity,:recovery,:pressure,:vx,:vy,:vz)
    arrays = map(names) do name
        a = KernelAbstractions.zeros(backend,Float32,grid.size...)
        name == :nutrient && fill!(a,1f0)
        a
    end
    scratch = map(_->KernelAbstractions.zeros(backend,Float32,grid.size...),names)
    World(grid,backend,NamedTuple{names}(arrays),NamedTuple{names}(scratch),0.0,Int(seed))
end
"""Advance processes in tuple order. Time advances once, after all processes.
The Float64 clock accumulates the actual Float32 integration timestep, preserving
long-run scheduling precision without changing field or kernel precision.
"""
function step!(w::World, processes::Tuple, dt::Real)
    isfinite(dt) && dt>0 || throw(ArgumentError("dt must be positive and finite"))
    isfinite(Float32(dt)) && Float32(dt)>0 || throw(ArgumentError("dt must remain positive and finite in Float32"))
    for p in processes
        check_dt(w,p,dt)
    end
    for p in processes
        advance!(w,p,Float32(dt))
    end
    w.time += Float64(Float32(dt))
    w
end
step!(w::World,p,dt::Real) = step!(w,(p,),dt)
material_mass(w::World) = sum(w.fields.nutrient) + sum(w.fields.deposit)
"""Seed a smooth, reproducible activity patch on host, then transfer to the backend."""
function seed_activity!(w::World; amplitude=0.9f0, radius=0.22f0, center=(0f0,0f0,0f0), noise=0.02f0)
    radius>0 || throw(ArgumentError("radius must be positive"))
    rng=MersenneTwister(w.seed)
    a=Array{Float32}(undef,w.grid.size)
    for k in axes(a,3), j in axes(a,2), i in axes(a,1)
        xyz=ntuple(d->(( (i,j,k)[d]-0.5f0)*w.grid.spacing[d]-w.grid.extent[d]/2),3)
        r2=sum((xyz[d]-center[d])^2 for d in 1:(w.grid.size[3]==1 ? 2 : 3))
        a[i,j,k]=amplitude*exp(-r2/(2radius^2))+noise*(rand(rng,Float32)-0.5f0)
    end
    copyto!(w.fields.activity,a)
    w
end
@inline function laplacian(a,i,j,k,h)
    nx,ny,nz=size(a)
    im=i==1 ? nx : i-1; ip=i==nx ? 1 : i+1
    jm=j==1 ? ny : j-1; jp=j==ny ? 1 : j+1
    v=(a[im,j,k]+a[ip,j,k]-2a[i,j,k])/h[1]^2 + (a[i,jm,k]+a[i,jp,k]-2a[i,j,k])/h[2]^2
    if nz>1
        km=k==1 ? nz : k-1; kp=k==nz ? 1 : k+1
        v+=(a[i,j,km]+a[i,j,kp]-2a[i,j,k])/h[3]^2
    end
    v
end
function diffusion_limit(w,D)
    D<=0 && return Inf32
    0.45f0/(D*sum(inv(w.grid.spacing[d]^2) for d in 1:(w.grid.size[3]==1 ? 2 : 3)))
end
stability_limit(w,p)=Inf32
function check_dt(w,p,dt)
    isfinite(dt) && dt>0 || throw(ArgumentError("dt must be positive and finite"))
    isfinite(Float32(dt)) && Float32(dt)>0 || throw(ArgumentError("dt must remain positive and finite in Float32"))
    validate(w,p)
    Float32(dt)<=stability_limit(w,p) || throw(ArgumentError("dt=$(Float32(dt)) exceeds explicit stability limit $(stability_limit(w,p)) for $(typeof(p))"))
end
function launch!(w,kernel,args...)
    kernel(w.backend,256)(args...; ndrange=w.grid.size)
    KernelAbstractions.synchronize(w.backend)
end

validate(p)=nothing
validate(w,p)=validate(p)
