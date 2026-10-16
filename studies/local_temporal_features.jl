"""Small, device-shaped temporal features for native study policies.

No policy, reward, physical inventory or production controller is changed here.
Call `observe!` inside a decision kernel to avoid a separate launch. Inputs are
the native bounded material observations, in [-1,1]. State follows resident IDs.
"""
module LocalTemporalFeatures
using KernelAbstractions

export coefficients, state, observe!, contrast, observe_kernel!, inherit_kernel!

"Exact exponential coefficients for one physical observation interval."
function coefficients(interval::Real; fast::Real=.4, slow::Real=2.)
    all(isfinite,(interval,fast,slow)) || throw(ArgumentError("finite clocks required"))
    interval>0 && 0<fast<slow || throw(ArgumentError("require interval>0 and 0<fast<slow"))
    (fast=Float32(-expm1(-interval/fast)),slow=Float32(-expm1(-interval/slow)))
end

"Allocate with the same device/storage family as the native input array."
function state(prototype,capacity::Integer)
    capacity>0 || throw(ArgumentError("positive capacity required"))
    (ids=fill!(similar(prototype,Int32,capacity),0),
     fast=fill!(similar(prototype,Float32,capacity,3),0f0),
     slow=fill!(similar(prototype,Float32,capacity,3),0f0))
end

"Observe once at a decision boundary; first observation creates no false trend."
@inline function observe!(z,i,id,s1,s2,s3,af,as)
    changed=z.ids[i]!=id
    for k in 1:3
        x=k==1 ? s1 : k==2 ? s2 : s3
        if changed
            z.fast[i,k]=x;z.slow[i,k]=x
        else
            z.fast[i,k]+=af*(x-z.fast[i,k])
            z.slow[i,k]+=as*(x-z.slow[i,k])
        end
    end
    z.ids[i]=id
    nothing
end

"Bounded signed temporal contrast, in [-1,1] for bounded source observations."
@inline contrast(z,i,k)=(z.fast[i,k]-z.slow[i,k])*.5f0

@kernel function observe_kernel!(z,active,ids,inputs,firstcol,af,as)
    i=@index(Global,Linear)
    if active[i]==1
        observe!(z,i,ids[i],inputs[i,firstcol],inputs[i,firstcol+1],inputs[i,firstcol+2],af,as)
    end
end

"""Run after native division and before the child's first observation.

copy_parent=false explicitly resets the newborn's temporal state. When true,
copy normalized information (not inventory) only from a matching parent ID.
Only the child is written, so the parent has no partition/race side effect.
Native division assigns at most one child per parent in a step; no child may
also be a parent in this same inheritance pass. Non-division identity changes
are handled independently by observe!.
"""
@kernel function inherit_kernel!(z,active,ids,new_parent,copy_parent)
    child=@index(Global,Linear)
    parent=new_parent[child]
    if active[child]==1 && parent>0
        valid=copy_parent && active[parent]==1 && z.ids[parent]==ids[parent]
        for k in 1:3
            z.fast[child,k]=valid ? z.fast[parent,k] : 0f0
            z.slow[child,k]=valid ? z.slow[parent,k] : 0f0
        end
        z.ids[child]=valid ? ids[child] : Int32(0)
    end
end
end
