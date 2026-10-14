"""CPU observer for smooth graded-index fields. No interface Snell/Fresnel law.
Coordinates use Grid.extent; sampling is cell-centered trilinear with its analytic
piecewise gradient. Outside the box the index is 1 and gradient zero.
"""
struct GradedIndex{A}
    index::A
    extent::NTuple{3,Float64}
    function GradedIndex(index::AbstractArray{<:Real,3};extent=(2.,2.,2.))
        index isa StridedArray || throw(ArgumentError("GRIN observer requires a CPU strided array; transfer GPU state explicitly with Array"))
        all(>(0),size(index)) || throw(ArgumentError("index dimensions must be nonempty"))
        all(x->isfinite(x)&&x>0,index) || throw(ArgumentError("index must be positive and finite"))
        ex=Tuple(Float64.(extent))
        length(ex)==3 && all(x->isfinite(x)&&x>0,ex) || throw(ArgumentError("extent must be positive and finite"))
        new{typeof(index)}(index,ex)
    end
end
function validate_grin(f::GradedIndex)
    all(x->isfinite(x)&&x>0,f.index) || throw(ArgumentError("index must remain positive and finite"))
    nothing
end
@inline function index_gradient(f::GradedIndex,p)
    any(d->abs(p[d])>f.extent[d]/2,1:3) && return (1.,(0.,0.,0.))
    dims=size(f.index)
    q=ntuple(d->clamp((p[d]/f.extent[d]+.5)*dims[d]+.5,1.,Float64(dims[d])),3)
    lo=ntuple(d->floor(Int,q[d]),3); hi=ntuple(d->min(lo[d]+1,dims[d]),3)
    w=ntuple(d->q[d]-lo[d],3)
    scale=ntuple(d->(q[d]>1 && q[d]<dims[d]) ? dims[d]/f.extent[d] : 0.,3)
    n=0.; gx=0.;gy=0.;gz=0.
    @inbounds for k in 0:1,j in 0:1,i in 0:1
        bx=i==0 ? 1-w[1] : w[1]; by=j==0 ? 1-w[2] : w[2]; bz=k==0 ? 1-w[3] : w[3]
        a=Float64(f.index[i==0 ? lo[1] : hi[1],j==0 ? lo[2] : hi[2],k==0 ? lo[3] : hi[3]])
        n+=a*bx*by*bz
        gx+=a*(i==0 ? -1 : 1)*scale[1]*by*bz
        gy+=a*bx*(j==0 ? -1 : 1)*scale[2]*bz
        gz+=a*bx*by*(k==0 ? -1 : 1)*scale[3]
    end
    n,(gx,gy,gz)
end
@inline function _grin_acceleration(f,p,v)
    n,g=index_gradient(f,p); parallel=_dot(v,g)
    ntuple(d->(g[d]-v[d]*parallel)/n,3)
end
"""Arc-length midpoint integration of dr/ds=v, dv/ds=(∇n-v(v⋅∇n))/n.
Normalizes direction after each step. Returns a bounded path; trapped/backward
rays are status=:unresolved, not silently continued as straight rays.
"""
function trace_grin(f::GradedIndex,origin,direction;ds=.01,target_z=f.extent[3]/2,max_length=8.,save_path=true,_validated=false)
    _validated || validate_grin(f)
    all(isfinite,origin) && all(isfinite,direction) && _dot(direction,direction)>0 || throw(ArgumentError("finite origin and nonzero direction required"))
    isfinite(ds) && ds>0 && isfinite(max_length) && max_length>0 && isfinite(target_z) || throw(ArgumentError("positive finite step and length required"))
    p=Tuple(Float64.(origin));v=_unit(Tuple(Float64.(direction)));length=0.;opl=0.
    path=NTuple{3,Float64}[]; save_path && push!(path,p)
    p[3]>=target_z && return (position=p,direction=v,arc_length=length,optical_length=opl,status=:reached,path=path)
    for _ in 1:ceil(Int,max_length/ds)
        h=min(ds,max_length-length);h<=0 && break
        a=_grin_acceleration(f,p,v)
        pm=_vadd(p,_vmul(v,h/2));vm=_unit(_vadd(v,_vmul(a,h/2)))
        am=_grin_acceleration(f,pm,vm)
        pn=_vadd(p,_vmul(vm,h));vn=_unit(_vadd(v,_vmul(am,h)))
        if pn[3]>=target_z && pn[3]>p[3]
            fraction=(target_z-p[3])/(pn[3]-p[3])
            h*=fraction;pn=_vadd(p,_vmul(_vadd(pn,_vmul(p,-1)),fraction))
            vn=_unit(_vadd(v,_vmul(_vadd(vn,_vmul(v,-1)),fraction)))
        end
        nm,_=index_gradient(f,_vmul(_vadd(p,pn),.5));opl+=nm*h;length+=h;p=pn;v=vn
        save_path && push!(path,p)
        p[3]>=target_z-1e-12 && return (position=p,direction=v,arc_length=length,optical_length=opl,status=:reached,path=path)
    end
    (position=p,direction=v,arc_length=length,optical_length=opl,status=:unresolved,path=path)
end
"""Orthographic background chart observed through the frozen index field.
Camera is at z=-extent_z/2; rays run toward +z. A background plane at background_z
lies outside the index volume. No scattering, emission, polarization or interfaces.
Returns linear RGB, hit coordinates, and unresolved ray count.
"""
function render_grin(f::GradedIndex;width=320,height=320,ds=.012,background_z=1.7,span=2.4)
    validate_grin(f)
    width>0 && height>0 && span>0 || throw(ArgumentError("positive image dimensions/span required"))
    background_z>f.extent[3]/2 || throw(ArgumentError("background must lie behind volume"))
    rgb=zeros(Float32,3,width,height);hits=zeros(Float64,2,width,height);failed=zeros(Int,height)
    Threads.@threads for y in 1:height
        for x in 1:width
            u=((x-.5)/width-.5)*span*width/height;v=(.5-(y-.5)/height)*span
            ray=trace_grin(f,(u,v,-f.extent[3]/2),(0.,0.,1.);ds=ds,save_path=false,_validated=true)
            if ray.status!=:reached || ray.direction[3]<=1e-8
                rgb[:,x,y].=(1f0,0f0,1f0);failed[y]+=1;continue
            end
            t=(background_z-ray.position[3])/ray.direction[3]
            point=_vadd(ray.position,_vmul(ray.direction,t));hits[:,x,y].=point[1:2]
            # A fixed luminous chart: fine grid over warm/cool checks and circles.
            a,b=point[1],point[2]; pitch=.12
            line=min(abs(mod(a+pitch/2,pitch)-pitch/2),abs(mod(b+pitch/2,pitch)-pitch/2))<.006
            check=isodd(floor(Int,a/.36)+floor(Int,b/.36))
            circle=abs(sqrt(a*a+b*b)-.68)<.012
            color=line ? (.04,.055,.07) : circle ? (1.,.25,.08) : check ? (.72,.82,.88) : (.9,.72,.42)
            rgb[:,x,y].=color
        end
    end
    (rgb=rgb,hits=hits,unresolved=sum(failed))
end
export GradedIndex,index_gradient,trace_grin,render_grin
