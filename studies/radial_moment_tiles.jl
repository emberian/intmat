module RadialMomentTiles
using Interior
export MomentGeometry, LazyMomentPyramid, MomentPyramid, SupportPlan, compile_support, refresh!, sample, storage_bytes, memory_forecast, scatter, allocate_requests

"""Experimental hierarchical polynomial gathers at unchanged physical radius.

Fields are borrowed read-only until refresh! rebuilds their moments. The radial
polynomial is exact in real arithmetic; moment reduction changes floating-point
order and is NOT a bitwise replacement for the production Float32 iterator.
Boundary leaves retain explicit periodic sphere tests and nearest-site fallback.
"""
abstract type AbstractMomentGeometry end
abstract type AbstractMomentField <: AbstractMomentGeometry end
struct MomentPyramid{T,A,G} <: AbstractMomentField
    field::A
    grid::G
    radius::Float64
    leaf_width::Int
    levels::Vector{Array{T,4}} # moment, tile x, tile y, tile z
end
struct MomentLevelShape
    shape::NTuple{4,Int}
end
Base.size(s::MomentLevelShape,d::Int)=s.shape[d]
struct MomentGeometry{G} <: AbstractMomentGeometry
    grid::G
    radius::Float64
    leaf_width::Int
    levels::Vector{MomentLevelShape}
end
function MomentGeometry(grid,radius;leaf_width=4)
    isfinite(radius)&&radius>0&&leaf_width>0 || throw(ArgumentError("positive finite geometry required"))
    dims=ntuple(d->cld(grid.size[d],leaf_width),3);levels=MomentLevelShape[]
    while true
        push!(levels,MomentLevelShape((14,dims...)))
        dims==(1,1,1) && break
        dims=ntuple(d->cld(dims[d],2),3)
    end
    MomentGeometry(grid,Float64(radius),Int(leaf_width),levels)
end
struct LazyMomentPyramid{T,A,G} <: AbstractMomentField
    field::A
    grid::G
    radius::Float64
    leaf_width::Int
    levels::Vector{MomentLevelShape}
    index::Vector{Array{Int32,3}}
    moments::Matrix{T}
end
@inline node_moments(p::MomentPyramid,l,i,j,k)=ntuple(a->Float64(p.levels[l][a,i,j,k]),Val(14))
@inline function node_moments(p::LazyMomentPyramid,l,i,j,k)
    slot=p.index[l][i,j,k]
    slot>0 || throw(ArgumentError("node not covered by lazy support plan set"))
    ntuple(a->Float64(p.moments[a,slot]),Val(14))
end
@inline m2(m,i,j)=m[i==j ? 4+i : i+j==3 ? 8 : i+j==4 ? 9 : 10]
@inline function basis(x,y,z)
    q=(x*x+y*y)+z*z
    (one(x),x,y,z,x*x,y*y,z*z,x*y,x*z,y*z,x*q,y*q,z*q,q*q)
end
@inline function translated(m,d)
    x,y,z=d;dd=(x*x+y*y)+z*z;tr=m[5]+m[6]+m[7]
    dot1=x*m[2]+y*m[3]+z*m[4]
    v=ntuple(i->m[1+i]+d[i]*m[1],3)
    second=ntuple(6) do k
        i,j=k<=3 ? (k,k) : k==4 ? (1,2) : k==5 ? (1,3) : (2,3)
        m2(m,i,j)+d[i]*m[1+j]+d[j]*m[1+i]+d[i]*d[j]*m[1]
    end
    third=ntuple(3) do i
        cross=sum(d[j]*m2(m,i,j) for j in 1:3)
        m[10+i]+d[i]*tr+2cross+2d[i]*dot1+dd*m[1+i]+d[i]*dd*m[1]
    end
    quadratic=x*x*m[5]+y*y*m[6]+z*z*m[7]+2(x*y*m[8]+x*z*m[9]+y*z*m[10])
    fourth=m[14]+4(x*m[11]+y*m[12]+z*m[13])+2dd*tr+4quadratic+4dd*dot1+dd*dd*m[1]
    (m[1],v...,second...,third...,fourth)
end
@inline function bounds(p,level,i,j,k)
    width=p.leaf_width*(1<<(level-1));dims=p.grid.size
    lo=((i-1)*width+1,(j-1)*width+1,(k-1)*width+1)
    hi=ntuple(d->min(dims[d],lo[d]+width-1),3)
    lo,hi
end
@inline function center(p,lo,hi)
    ntuple(d->((lo[d]+hi[d])/2-.5)*Float64(p.grid.extent[d])/p.grid.size[d],3)
end
function MomentPyramid(field,grid,radius;leaf_width=4,precision=Float64)
    size(field)==grid.size || throw(DimensionMismatch("moment field/grid differ"))
    isfinite(radius)&&radius>0&&leaf_width>0 || throw(ArgumentError("positive finite radius and tile width required"))
    precision in (Float32,Float64) || throw(ArgumentError("Float32 or Float64 moments required"))
    all(isfinite,field) || throw(ArgumentError("finite moment field required"))
    dims=ntuple(d->cld(grid.size[d],leaf_width),3);levels=Array{precision,4}[]
    while true
        push!(levels,zeros(precision,14,dims...))
        dims==(1,1,1) && break
        dims=ntuple(d->cld(dims[d],2),3)
    end
    p=MomentPyramid(field,grid,Float64(radius),Int(leaf_width),levels)
    refresh!(p)
end
function refresh!(p::MomentPyramid{T}) where T
    all(isfinite,p.field) || throw(ArgumentError("finite moment field required"))
    nd=p.grid.size[3]==1 ? 2 : 3;h=ntuple(d->Float64(p.grid.extent[d])/p.grid.size[d],3)
    leaf=p.levels[1];fill!(leaf,0)
    for k in axes(leaf,4),j in axes(leaf,3),i in axes(leaf,2)
        lo,hi=bounds(p,1,i,j,k);c=center(p,lo,hi)
        for z in lo[3]:hi[3],y in lo[2]:hi[2],x in lo[1]:hi[1]
            t=ntuple(d->d>nd ? zero(T) : T((( (x,y,z)[d]-.5)*h[d]-c[d])/p.radius),3)
            monomials=basis(t...);value=T(p.field[x,y,z])
            for a in 1:14;leaf[a,i,j,k]+=value*monomials[a];end
        end
    end
    for level in 2:length(p.levels)
        dst=p.levels[level];src=p.levels[level-1];fill!(dst,0)
        for k in axes(dst,4),j in axes(dst,3),i in axes(dst,2)
            lo,hi=bounds(p,level,i,j,k);c=center(p,lo,hi)
            for z in 2k-1:min(2k,size(src,4)),y in 2j-1:min(2j,size(src,3)),x in 2i-1:min(2i,size(src,2))
                clo,chi=bounds(p,level-1,x,y,z);cc=center(p,clo,chi)
                shift=ntuple(d->d>nd ? zero(T) : T((cc[d]-c[d])/p.radius),3)
                values=translated(ntuple(a->src[a,x,y,z],Val(14)),shift)
                for a in 1:14;dst[a,i,j,k]+=values[a];end
            end
        end
    end
    p
end
function LazyMomentPyramid(field,geometry::MomentGeometry,plans;precision=Float64)
    size(field)==geometry.grid.size || throw(DimensionMismatch("lazy field/grid differ"))
    precision in (Float32,Float64) || throw(ArgumentError("Float32 or Float64 moments required"))
    index=[zeros(Int32,s.shape[2:4]...) for s in geometry.levels]
    function mark(l,i,j,k)
        index[l][i,j,k]!=0 && return
        index[l][i,j,k]=-1
        if l>1
            child=index[l-1]
            for z in 2k-1:min(2k,size(child,3)),y in 2j-1:min(2j,size(child,2)),x in 2i-1:min(2i,size(child,1))
                mark(l-1,x,y,z)
            end
        end
    end
    for plan in plans
        plan.grid_size==geometry.grid.size && plan.extent==Float64.(geometry.grid.extent) &&
            plan.radius==geometry.radius && plan.leaf_width==geometry.leaf_width || throw(ArgumentError("lazy plan geometry differs"))
        for (l,i,j,k) in plan.blocks;l>0 && mark(Int(l),i,j,k);end
    end
    slots=0
    for level in index,I in eachindex(level)
        if level[I]!=0
            slots+=1;slots<=typemax(Int32) || throw(ArgumentError("lazy index overflow"));level[I]=slots
        end
    end
    p=LazyMomentPyramid(field,geometry.grid,geometry.radius,geometry.leaf_width,geometry.levels,index,zeros(precision,14,slots))
    refresh!(p)
end
function refresh!(p::LazyMomentPyramid{T}) where T
    # Only touched leaves are inspected; boundary values are read directly when sampling.
    fill!(p.moments,0);nd=p.grid.size[3]==1 ? 2 : 3
    h=ntuple(d->Float64(p.grid.extent[d])/p.grid.size[d],3)
    for l in eachindex(p.index),I in CartesianIndices(p.index[l])
        slot=p.index[l][I];slot==0 && continue
        i,j,k=Tuple(I);lo,hi=bounds(p,l,i,j,k);c=center(p,lo,hi)
        if l==1
            for z in lo[3]:hi[3],y in lo[2]:hi[2],x in lo[1]:hi[1]
                value=T(p.field[x,y,z]);isfinite(value) || throw(ArgumentError("finite lazy field required"))
                t=ntuple(d->d>nd ? zero(T) : T((((x,y,z)[d]-.5)*h[d]-c[d])/p.radius),3)
                b=basis(t...)
                for a in 1:14;p.moments[a,slot]+=value*b[a];end
            end
        else
            src=p.index[l-1]
            for z in 2k-1:min(2k,size(src,3)),y in 2j-1:min(2j,size(src,2)),x in 2i-1:min(2i,size(src,1))
                clo,chi=bounds(p,l-1,x,y,z);cc=center(p,clo,chi)
                shift=ntuple(d->d>nd ? zero(T) : T((cc[d]-c[d])/p.radius),3)
                childslot=src[x,y,z]
                values=translated(ntuple(a->p.moments[a,childslot],Val(14)),shift)
                for a in 1:14;p.moments[a,slot]+=values[a];end
            end
        end
    end
    p
end
storage_bytes(p::LazyMomentPyramid)=sizeof(p.moments)+sum(sizeof,p.index)
@inline function polynomial(m,d)
    x,y,z=d;A=1-((x*x+y*y)+z*z)
    quadratic=x*x*m[5]+y*y*m[6]+z*z*m[7]+2(x*y*m[8]+x*z*m[9]+y*z*m[10])
    A*A*m[1]-4A*(x*m[2]+y*m[3]+z*m[4])-2A*(m[5]+m[6]+m[7])+4quadratic+
        4(x*m[11]+y*m[12]+z*m[13])+m[14]
end
@inline function geometry_moments(p,lo,hi)
    nd=p.grid.size[3]==1 ? 2 : 3
    n=ntuple(d->Float64(hi[d]-lo[d]+1),3)
    step=ntuple(d->d>nd ? 0. : Float64(p.grid.extent[d])/p.grid.size[d]/p.radius,3)
    s2=ntuple(d->step[d]^2*n[d]*(n[d]^2-1)/12,3)
    s4=ntuple(d->step[d]^4*n[d]*(n[d]^2-1)*(3n[d]^2-7)/240,3)
    xx=s2[1]*n[2]*n[3];yy=s2[2]*n[1]*n[3];zz=s2[3]*n[1]*n[2]
    fourth=s4[1]*n[2]*n[3]+s4[2]*n[1]*n[3]+s4[3]*n[1]*n[2]+
        2(s2[1]*s2[2]*n[3]+s2[1]*s2[3]*n[2]+s2[2]*s2[3]*n[1])
    (prod(n),0.,0.,0.,xx,yy,zz,0.,0.,0.,0.,0.,0.,fourth)
end
mutable struct Accumulator
    numerator::Float64
    normalization::Float64
    nodes::Int
    full_nodes::Int
    direct_voxels::Int
    accepted_voxels::Int
    blocks::Union{Nothing,Vector{NTuple{4,Int32}}}
end
struct SupportPlan
    grid_size::NTuple{3,Int}
    extent::NTuple{3,Float64}
    radius::Float64
    leaf_width::Int
    position::NTuple{3,Float64}
    blocks::Vector{NTuple{4,Int32}} # positive level=interior; negative=boundary leaf
    normalization::Float64
    accepted_voxels::Int
    fallback::Bool
    axis_start::NTuple{3,Int}
    axis_distance2::NTuple{3,Vector{Float64}}
end
@inline periodic_delta(x,L)=mod(x+L/2,L)-L/2
@inline function point_distance2(p,position,site,::Nothing)
    nd=p.grid.size[3]==1 ? 2 : 3
    sum(periodic_delta((site[d]-.5)*(Float64(p.grid.extent[d])/p.grid.size[d])-position[d],Float64(p.grid.extent[d]))^2 for d in 1:nd)
end
@inline function point_distance2(p,position,site,cache)
    distance=ntuple(d->cache.axis_distance2[d][mod(site[d]-cache.axis_start[d],p.grid.size[d])+1],3)
    (distance[1]+distance[2])+distance[3]
end
function visit!(acc,p,level,i,j,k,position,gather::Val{G},cache=nothing) where G
    acc.nodes+=1
    lo,hi=bounds(p,level,i,j,k);c=center(p,lo,hi)
    nd=p.grid.size[3]==1 ? 2 : 3
    delta=ntuple(d->d>nd ? 0. : periodic_delta(c[d]-position[d],Float64(p.grid.extent[d])),3)
    span=ntuple(d->d>nd ? 0. : (hi[d]-lo[d])*Float64(p.grid.extent[d])/p.grid.size[d]/2,3)
    padding=64eps(Float64)*max(maximum(Float64.(p.grid.extent)),p.radius)
    lower=sum(max(0.,abs(delta[d])-span[d]-padding)^2 for d in 1:nd)
    lower>=p.radius^2 && return nothing
    upper=sum((abs(delta[d])+span[d]+padding)^2 for d in 1:nd)
    continuous=all(abs(delta[d])+span[d]+padding<Float64(p.grid.extent[d])/2 for d in 1:nd)
    if continuous && upper<p.radius^2
        offset=ntuple(d->delta[d]/p.radius,3)
        if G
            m=node_moments(p,level,i,j,k)
            acc.numerator+=polynomial(m,offset)
        end
        acc.blocks===nothing || push!(acc.blocks,Int32.((level,i,j,k)))
        acc.normalization+=polynomial(geometry_moments(p,lo,hi),offset)
        acc.full_nodes+=1;acc.accepted_voxels+=prod(hi[d]-lo[d]+1 for d in 1:3)
    elseif level==1
        h=ntuple(d->Float64(p.grid.extent[d])/p.grid.size[d],3)
        count_before=acc.accepted_voxels
        for z in lo[3]:hi[3],y in lo[2]:hi[2],x in lo[1]:hi[1]
            acc.direct_voxels+=1
            r2=point_distance2(p,position,(x,y,z),cache)
            if r2<p.radius^2
                weight=(1-r2/p.radius^2)^2
                G && (acc.numerator+=weight*Float64(p.field[x,y,z]))
                acc.normalization+=weight;acc.accepted_voxels+=1
            end
        end
        if acc.blocks!==nothing && acc.accepted_voxels>count_before
            push!(acc.blocks,Int32.((-level,i,j,k)))
        end
    else
        src=p.levels[level-1]
        for z in 2k-1:min(2k,size(src,4)),y in 2j-1:min(2j,size(src,3)),x in 2i-1:min(2i,size(src,2))
            visit!(acc,p,level-1,x,y,z,position,gather,cache)
        end
    end
    nothing
end
function sample(p::AbstractMomentField,position)
    length(position)==3 && all(isfinite,position) || throw(ArgumentError("finite physical xyz required"))
    xyz=ntuple(d->mod(Float64(position[d]),Float64(p.grid.extent[d])),3)
    acc=Accumulator(0.,0.,0,0,0,0,nothing)
    visit!(acc,p,length(p.levels),1,1,1,xyz,Val(true))
    fallback=acc.normalization==0
    if fallback
        site=ntuple(d->min(p.grid.size[d],floor(Int,xyz[d]/Float64(p.grid.extent[d])*p.grid.size[d])+1),3)
        acc.numerator=Float64(p.field[site...]);acc.normalization=1.;acc.accepted_voxels=1
    end
    acc.normalization>0 && isfinite(acc.numerator) || throw(ArgumentError("unresolved moment cancellation/nonfinite query"))
    (;value=acc.numerator/acc.normalization,normalization=acc.normalization,
        nodes=acc.nodes,full_nodes=acc.full_nodes,direct_voxels=acc.direct_voxels,
        accepted_voxels=acc.accepted_voxels,fallback)
end
"""Compile field-independent interior/boundary tiles once per held position.

The plan owns its position snapshot; movement, birth, grid/radius or tiling changes
require a new plan. Field mutations require refreshed moments but not a new plan.
No production support rebuild or lifecycle schedule is changed by this helper.
"""
function compile_support(p::AbstractMomentGeometry,position)
    length(position)==3 && all(isfinite,position) || throw(ArgumentError("finite physical xyz required"))
    xyz=ntuple(d->mod(Float64(position[d]),Float64(p.grid.extent[d])),3)
    nd=p.grid.size[3]==1 ? 2 : 3
    starts=ntuple(3) do d
        h=Float64(p.grid.extent[d])/p.grid.size[d]
        d>nd || 2p.radius>=Float64(p.grid.extent[d]) ? 1 : floor(Int,(xyz[d]-p.radius)/h+.5)-p.leaf_width
    end
    distances=ntuple(3) do d
        h=Float64(p.grid.extent[d])/p.grid.size[d];L=Float64(p.grid.extent[d])
        count=d>nd ? 1 : 2p.radius>=L ? p.grid.size[d] : min(p.grid.size[d],2ceil(Int,p.radius/h)+2+2p.leaf_width)
        [d>nd ? 0. : periodic_delta((mod1(starts[d]+j-1,p.grid.size[d])-.5)*h-xyz[d],L)^2 for j in 1:count]
    end
    cache=(axis_start=starts,axis_distance2=distances)
    acc=Accumulator(0.,0.,0,0,0,0,NTuple{4,Int32}[])
    visit!(acc,p,length(p.levels),1,1,1,xyz,Val(false),cache)
    fallback=acc.normalization==0
    fallback && (acc.normalization=1.;acc.accepted_voxels=1)
    acc.normalization>0 || throw(ArgumentError("unresolved moment normalization"))
    SupportPlan(p.grid.size,Float64.(p.grid.extent),p.radius,p.leaf_width,xyz,acc.blocks,
        acc.normalization,acc.accepted_voxels,fallback,starts,distances)
end
@inline function boundary_value(p,plan,lo,hi)
    total=0.;r2max=p.radius^2
    nz=p.grid.size[3];ny=p.grid.size[2];nx=p.grid.size[1]
    iz=mod(lo[3]-plan.axis_start[3],nz)+1
    @inbounds for z in lo[3]:hi[3]
        dz=plan.axis_distance2[3][iz]
        iy=mod(lo[2]-plan.axis_start[2],ny)+1
        for y in lo[2]:hi[2]
            dy=plan.axis_distance2[2][iy]
            ix=mod(lo[1]-plan.axis_start[1],nx)+1
            for x in lo[1]:hi[1]
                dx=plan.axis_distance2[1][ix];r2=(dx+dy)+dz
                r2<r2max && (total+=(1-r2/r2max)^2*Float64(p.field[x,y,z]))
                ix=ix==nx ? 1 : ix+1
            end
            iy=iy==ny ? 1 : iy+1
        end
        iz=iz==nz ? 1 : iz+1
    end
    total
end
@inline function block_value(p,level,i,j,k,position)
    lo,hi=bounds(p,level,i,j,k);c=center(p,lo,hi)
    nd=p.grid.size[3]==1 ? 2 : 3
    offset=ntuple(d->d>nd ? 0. : periodic_delta(c[d]-position[d],Float64(p.grid.extent[d]))/p.radius,3)
    polynomial(node_moments(p,level,i,j,k),offset)
end
function sample(p::AbstractMomentField,plan::SupportPlan)
    p.grid.size==plan.grid_size && Float64.(p.grid.extent)==plan.extent &&
        p.radius==plan.radius && p.leaf_width==plan.leaf_width || throw(ArgumentError("moment support plan geometry differs"))
    if plan.fallback
        site=ntuple(d->min(p.grid.size[d],floor(Int,plan.position[d]/plan.extent[d]*p.grid.size[d])+1),3)
        return Float64(p.field[site...])
    end
    numerator=0.;nd=p.grid.size[3]==1 ? 2 : 3
    h=ntuple(d->plan.extent[d]/p.grid.size[d],3)
    for (signed_level,i,j,k) in plan.blocks
        level=abs(Int(signed_level))
        if signed_level>0
            numerator+=block_value(p,level,i,j,k,plan.position)
        else
            lo,hi=bounds(p,level,i,j,k)
            numerator+=boundary_value(p,plan,lo,hi)
        end
    end
    numerator/plan.normalization
end
@inline function kernel_coefficients(d)
    x,y,z=d;A=1-((x*x+y*y)+z*z)
    (A*A,-4A*x,-4A*y,-4A*z,-2A+4x*x,-2A+4y*y,-2A+4z*z,
        8x*y,8x*z,8y*z,4x,4y,4z,1.)
end
@inline cmatrix(c,i,j)=i==j ? c[4+i] : m2(c,i,j)/2
@inline function shifted_coefficients(c,delta)
    x,y,z=delta;q=(x*x+y*y)+z*z
    dotd=x*c[11]+y*c[12]+z*c[13]
    quadratic=x*x*c[5]+y*y*c[6]+z*z*c[7]+x*y*c[8]+x*z*c[9]+y*z*c[10]
    constant=c[1]+x*c[2]+y*c[3]+z*c[4]+quadratic+dotd*q+c[14]*q*q
    linear=ntuple(3) do i
        c[1+i]+2sum(cmatrix(c,i,j)*delta[j] for j in 1:3)+c[10+i]*q+
            2dotd*delta[i]+4c[14]*q*delta[i]
    end
    second=ntuple(6) do k
        if k<=3
            c[4+k]+dotd+2c[10+k]*delta[k]+2c[14]*q+4c[14]*delta[k]^2
        else
            i,j=k==4 ? (1,2) : k==5 ? (1,3) : (2,3)
            c[4+k]+2(c[10+i]*delta[j]+c[10+j]*delta[i])+8c[14]*delta[i]*delta[j]
        end
    end
    third=ntuple(i->c[10+i]+4c[14]*delta[i],3)
    (constant,linear...,second...,third...,c[14])
end

"""Mathematical transpose of normalized gathers, with changed Float64 order.

Interior cell contributions accumulate fourteen polynomial coefficients per node;
coarse coefficients propagate to children, then each leaf voxel has one writer.
Boundary tiles keep direct sphere tests. This is an additive REQUEST operator,
not paid catalyst delivery: production CAS rounding/credit attribution is not
silently replaced. Negative roundoff is recomputed directly, never clipped.
"""
function scatter(p::MomentPyramid,plans,amounts)
    length(plans)==length(amounts) && all(a->isfinite(a)&&a>=0,amounts) || throw(ArgumentError("finite nonnegative amount per support plan required"))
    all(plan->plan.grid_size==p.grid.size && plan.extent==Float64.(p.grid.extent) &&
        plan.radius==p.radius && plan.leaf_width==p.leaf_width,plans) || throw(ArgumentError("scatter plan geometry differs"))
    coefficients=[zeros(Float64,size(level)) for level in p.levels]
    output=zeros(Float64,p.grid.size);nd=p.grid.size[3]==1 ? 2 : 3
    for (plan,amount) in zip(plans,amounts)
        amount==0 && continue
        if plan.fallback
            site=ntuple(d->min(p.grid.size[d],floor(Int,plan.position[d]/plan.extent[d]*p.grid.size[d])+1),3)
            output[site...]+=amount;continue
        end
        for (signed_level,i,j,k) in plan.blocks
            level=abs(Int(signed_level));lo,hi=bounds(p,level,i,j,k)
            if signed_level>0
                centerpoint=center(p,lo,hi)
                delta=ntuple(d->d>nd ? 0. : periodic_delta(centerpoint[d]-plan.position[d],plan.extent[d])/p.radius,3)
                c=kernel_coefficients(delta)
                for a in 1:14;coefficients[level][a,i,j,k]+=amount/plan.normalization*c[a];end
            else
                for z in lo[3]:hi[3],y in lo[2]:hi[2],x in lo[1]:hi[1]
                    r2=point_distance2(p,plan.position,(x,y,z),plan)
                    r2<p.radius^2 && (output[x,y,z]+=amount/plan.normalization*(1-r2/p.radius^2)^2)
                end
            end
        end
    end
    for level in length(coefficients):-1:2
        src=coefficients[level];dst=coefficients[level-1]
        for k in axes(src,4),j in axes(src,3),i in axes(src,2)
            c=ntuple(a->src[a,i,j,k],Val(14));all(iszero,c) && continue
            lo,hi=bounds(p,level,i,j,k);parent_center=center(p,lo,hi)
            for z in 2k-1:min(2k,size(dst,4)),y in 2j-1:min(2j,size(dst,3)),x in 2i-1:min(2i,size(dst,2))
                clo,chi=bounds(p,level-1,x,y,z);child_center=center(p,clo,chi)
                delta=ntuple(d->d>nd ? 0. : (child_center[d]-parent_center[d])/p.radius,3)
                shifted=shifted_coefficients(c,delta)
                for a in 1:14;dst[a,x,y,z]+=shifted[a];end
            end
        end
    end
    leaf=coefficients[1];h=ntuple(d->Float64(p.grid.extent[d])/p.grid.size[d],3)
    for k in axes(leaf,4),j in axes(leaf,3),i in axes(leaf,2)
        c=ntuple(a->leaf[a,i,j,k],Val(14));all(iszero,c) && continue
        lo,hi=bounds(p,1,i,j,k);tile_center=center(p,lo,hi)
        for z in lo[3]:hi[3],y in lo[2]:hi[2],x in lo[1]:hi[1]
            t=ntuple(d->d>nd ? 0. : (((x,y,z)[d]-.5)*h[d]-tile_center[d])/p.radius,3)
            b=basis(t...);output[x,y,z]+=sum(c[a]*b[a] for a in 1:14)
        end
    end
    fallback_sites=0
    for site in CartesianIndices(output)
        isfinite(output[site]) || throw(ArgumentError("nonfinite polynomial request accumulation"))
        if output[site]<0
            fallback_sites+=1;actual=0.
            for (plan,amount) in zip(plans,amounts)
                r2=point_distance2(p,plan.position,Tuple(site),nothing)
                r2<p.radius^2 && (actual+=amount/plan.normalization*(1-r2/p.radius^2)^2)
            end
            output[site]=actual
        end
    end
    (;field=output,roundoff_fallback_sites=fallback_sites,workspace_bytes=sum(sizeof,coefficients))
end
"""CPU accounting proof for the synchronous request operator, not a cell step.

Withdraw actual representable Float32 field endpoints first. A finite escrow
funds returned Float64 grants; unassigned withdrawal is explicit inert mass and
potential heat. Inputs are not mutated. Production cell-store endpoint credits,
shared work budgets, CAS catalyst exports and lifecycle are deliberately absent.
"""
function allocate_requests(field::AbstractArray{Float32,3},p,plans,requests;potential=1.)
    size(field)==p.grid.size && all(x->isfinite(x)&&x>=0,field) || throw(ArgumentError("finite nonnegative matching Float32 source required"))
    isfinite(potential)&&potential>=0 || throw(ArgumentError("finite nonnegative chemical potential required"))
    demand=scatter(p,plans,requests)
    output=copy(field);scale=zeros(Float64,p.grid.size);actual=zeros(Float64,p.grid.size)
    for i in eachindex(output)
        before=Float64(field[i]);wanted=min(before,demand.field[i])
        after=Float32(max(0.,before-wanted))
        before-Float64(after)>wanted && (after=nextfloat(after))
        output[i]=after;actual[i]=before-Float64(after)
        scale[i]=demand.field[i]>0 ? actual[i]/demand.field[i] : 0.
    end
    ratio=MomentPyramid(scale,p.grid,p.radius;leaf_width=p.leaf_width)
    withdrawn=sum(actual);escrow=withdrawn;grants=zeros(Float64,length(plans));request_caps=0
    for i in eachindex(plans)
        want=requests[i]*sample(ratio,plans[i])
        isfinite(want)&&want>=0 || throw(ArgumentError("unresolved negative/nonfinite allocation moment"))
        want>requests[i] && (request_caps+=1)
        wanted=min(want,requests[i],escrow)
        after=max(0.,escrow-wanted)
        escrow-after>wanted && (after=nextfloat(after))
        grants[i]=escrow-after;escrow=after
    end
    (;field=output,grants,withdrawn,inert_residual=escrow,chemical_heat=escrow*potential,
        request_cap_count=request_caps,roundoff_fallback_sites=demand.roundoff_fallback_sites)
end
storage_bytes(p::SupportPlan)=sizeof(p.blocks)+sum(sizeof,p.axis_distance2)
storage_bytes(p::MomentPyramid)=sum(sizeof,p.levels)
function memory_forecast(dims;leaf_width=4,precision_bytes=4)
    tile=ntuple(d->cld(dims[d],leaf_width),3);nodes=0
    while true
        nodes+=prod(tile);tile==(1,1,1) && break;tile=ntuple(d->cld(tile[d],2),3)
    end
    (;moment_bytes=14nodes*precision_bytes,nodes,
        base_grid_bytes=84prod(dims),coupled_grid_bytes=268prod(dims))
end
end
