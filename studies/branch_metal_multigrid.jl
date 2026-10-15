# Isolated KA/Metal multigrid support; not loaded by the live CPU producer.
using Interior
using Interior: KernelAbstractions
using .KernelAbstractions: @kernel, @index, @localmem, @synchronize
struct DeviceMG{A,M}
 u::A;rhs::A;residual::A;occupied::M
end
@kernel function branch_rb_kernel!(u,b,mask,n,color,omega)
 q,j,k=@index(Global,NTuple)
 i=2q+mod(color-j-k,2)
 @inbounds begin
 if 1<i<n && 1<j<n && 1<k<n && (i+j+k)%2==color && mask[i,j,k]==0
  avg=(u[i-1,j,k]+u[i+1,j,k]+u[i,j-1,k]+u[i,j+1,k]+u[i,j,k-1]+u[i,j,k+1]+b[i,j,k])*(1f0/6f0)
  u[i,j,k]=muladd(omega,avg-u[i,j,k],u[i,j,k])
 end
 end
end
@kernel function branch_residual_kernel!(r,u,b,mask,n)
 i,j,k=@index(Global,NTuple)
 @inbounds begin
 if 1<i<n && 1<j<n && 1<k<n && mask[i,j,k]==0
  r[i,j,k]=b[i,j,k]+u[i-1,j,k]+u[i+1,j,k]+u[i,j-1,k]+u[i,j+1,k]+u[i,j,k-1]+u[i,j,k+1]-6f0*u[i,j,k]
 else
  r[i,j,k]=0f0
 end
 end
end
@kernel function branch_restrict_kernel!(b,r,mask,n)
 i,j,k=@index(Global,NTuple)
 @inbounds begin
 if 1<i<n && 1<j<n && 1<k<n && mask[i,j,k]==0
  s=0f0
  for dz in -1:1,dy in -1:1,dx in -1:1
   weight=(dx==0 ? 2f0 : 1f0)*(dy==0 ? 2f0 : 1f0)*(dz==0 ? 2f0 : 1f0)
   s+=weight*r[2i-1+dx,2j-1+dy,2k-1+dz]
  end
  b[i,j,k]=s/16f0
 else
  b[i,j,k]=0f0
 end
 end
end
@kernel function branch_prolong_kernel!(u,c,mask,n)
 i,j,k=@index(Global,NTuple)
 @inbounds begin
 if 1<i<n && 1<j<n && 1<k<n && mask[i,j,k]==0
  ci=(i+1)÷2;cj=(j+1)÷2;ck=(k+1)÷2
  tx=iseven(i) ? .5f0 : 0f0;ty=iseven(j) ? .5f0 : 0f0;tz=iseven(k) ? .5f0 : 0f0
  s=0f0
  for dz in 0:1,dy in 0:1,dx in 0:1
   s+=(dx==0 ? 1f0-tx : tx)*(dy==0 ? 1f0-ty : ty)*(dz==0 ? 1f0-tz : tz)*c[ci+dx,cj+dy,ck+dz]
  end
  u[i,j,k]+=s
 end
 end
end
@kernel function branch_coarse_kernel!(u,b,mask,n)
 q=@index(Global,Linear)
 if q==1
  @inbounds for sweep in 1:40,color in 0:1,k in 2:n-1,j in 2:n-1,i in 2:n-1
   if (i+j+k)%2==color && mask[i,j,k]==0
    u[i,j,k]=(u[i-1,j,k]+u[i+1,j,k]+u[i,j-1,k]+u[i,j+1,k]+u[i,j,k-1]+u[i,j,k+1]+b[i,j,k])*(1f0/6f0)
   end
  end
 end
end
@kernel function branch_mark_kernel!(u,mask,i,j,k,value)
 q=@index(Global,Linear)
 if q==1;mask[i,j,k]=0x01;u[i,j,k]=value;end
end
@kernel function branch_gather_kernel!(values,u,indices)
 q=@index(Global,Linear)
 values[q]=u[indices[q]]
end
# The convergence residual needs only its maximum, unlike the residual field
# consumed by restriction. This optional scratch avoids a full-grid write/read.
struct DeviceMGNorm{V}
 maxima::V
end
DeviceMGNorm(level::DeviceMG)=DeviceMGNorm(similar(level.u,Float32,cld(length(level.u),256)))
@kernel function branch_residual_norm_kernel!(maxima,u,b,mask,n)
 q=@index(Global,Linear);lane=@index(Local,Linear);block=@index(Group,Linear)
 scratch=@localmem Float32 (256,)
 value=0f0
 if q<=length(u)
  i=mod(q-1,n)+1;j=mod((q-1)÷n,n)+1;k=(q-1)÷(n*n)+1
  @inbounds if 1<i<n && 1<j<n && 1<k<n && mask[i,j,k]==0
   value=abs(b[i,j,k]+u[i-1,j,k]+u[i+1,j,k]+u[i,j-1,k]+u[i,j+1,k]+u[i,j,k-1]+u[i,j,k+1]-6f0*u[i,j,k])
  end
 end
 @inbounds scratch[lane]=value
 @synchronize
 for offset in (128,64,32,16,8,4,2,1)
  @inbounds if lane<=offset;scratch[lane]=max(scratch[lane],scratch[lane+offset]);end
  @synchronize
 end
 @inbounds if lane==1;maxima[block]=scratch[1];end
end
function device_residual_norm!(workspace::DeviceMGNorm,level)
 n=size(level.u,1);blocks=cld(length(level.u),256)
 length(workspace.maxima)==blocks || throw(DimensionMismatch("residual norm workspace/grid mismatch"))
 backend=KernelAbstractions.get_backend(level.u)
 KernelAbstractions.get_backend(workspace.maxima)==backend || throw(ArgumentError("residual norm workspace backend mismatch"))
 branch_residual_norm_kernel!(backend,256)(workspace.maxima,level.u,level.rhs,level.occupied,n;ndrange=blocks*256)
 Float64(maximum(workspace.maxima))/6
end

function device_levels(cpu,upload)
 [DeviceMG(upload(l.u),upload(l.rhs),upload(l.residual),upload(UInt8.(l.occupied))) for l in cpu]
end
function device_rb!(level,sweeps;omega=1f0,workgroup=(8,8,4))
 backend=KernelAbstractions.get_backend(level.u);n=size(level.u,1)
 for sweep in 1:sweeps,color in 0:1
  branch_rb_kernel!(backend,workgroup)(level.u,level.rhs,level.occupied,n,color,omega;ndrange=((n-1)÷2,n,n))
 end
end
function device_residual!(level;workgroup=(8,8,4))
 backend=KernelAbstractions.get_backend(level.u);n=size(level.u,1)
 branch_residual_kernel!(backend,workgroup)(level.residual,level.u,level.rhs,level.occupied,n;ndrange=(n,n,n))
end
function device_vcycle!(levels,l=1;workgroup=(8,8,4))
 f=levels[l];backend=KernelAbstractions.get_backend(f.u)
 if l==length(levels)
  branch_coarse_kernel!(backend,1)(f.u,f.rhs,f.occupied,size(f.u,1);ndrange=1);return
 end
 device_rb!(f,3;workgroup);device_residual!(f;workgroup);c=levels[l+1];fill!(c.u,0);nc=size(c.u,1)
 branch_restrict_kernel!(backend,workgroup)(c.rhs,f.residual,c.occupied,nc;ndrange=(nc,nc,nc))
 device_vcycle!(levels,l+1;workgroup);n=size(f.u,1)
 branch_prolong_kernel!(backend,workgroup)(f.u,c.u,f.occupied,n;ndrange=(n,n,n))
 device_rb!(f,3;workgroup)
end
"""Converged harmonic solve. `workgroup` changes execution layout only.
The historical (8,8,4) default preserves existing callers; (32,4,2) is the
measured AMD layout. Select it explicitly for a fresh run, not a live migration.
"""
function device_solve!(levels;tol=1e-6,maxcycles=200,workgroup=(8,8,4),norm_workspace=nothing)
 r=Inf
 for it in 1:maxcycles
  device_vcycle!(levels;workgroup)
  if norm_workspace===nothing
   device_residual!(first(levels);workgroup)
   r=Float64(maximum(abs,first(levels).residual))/6
  else
   r=device_residual_norm!(norm_workspace,first(levels))
  end
  r<=tol && return (cycles=it,residual=r)
 end
 error("device multigrid residual $r exceeds $tol")
end
function device_mark!(levels,I)
 backend=KernelAbstractions.get_backend(first(levels).u)
 for l in eachindex(levels)
  level=levels[l]
  if l>1
   n=size(level.u,1);I=CartesianIndex(ntuple(d->clamp(round(Int,(I[d]+1)/2),2,n-1),3))
  end
  branch_mark_kernel!(backend,1)(level.u,level.occupied,Tuple(I)...,l==1 ? 1f0 : 0f0;ndrange=1)
 end
end
