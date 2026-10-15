# Single-workgroup patch prototype, kept separate from production MG.
include(joinpath(@__DIR__,"electrical_local_correction.jl"))
using Interior: KernelAbstractions
using .KernelAbstractions: @kernel, @index, @localmem, @synchronize

# Radius8 plus a one-cell read-only halo: 19^3*(4+1)=34,295 bytes LDS.
# Each color is independent, so CPU sweep order and parallel color order agree.
@kernel function electrical_patch8_kernel!(u,occupied,norm,ci,cj,ck,n,sweeps)
 lane=@index(Local,Linear)
 values=@localmem Float32 (6859,)
 fixed=@localmem UInt8 (6859,)
 maxima=@localmem Float32 (256,)
 @inbounds begin
 for q in lane:256:6859
  x=mod(q-1,19)+1;y=mod((q-1)÷19,19)+1;z=(q-1)÷361+1
  i=ci+x-10;j=cj+y-10;k=ck+z-10
  inside=1<=i<=n && 1<=j<=n && 1<=k<=n
  values[q]=inside ? u[i,j,k] : 0f0
  fixed[q]=!inside || i==1 || j==1 || k==1 || i==n || j==n || k==n || occupied[i,j,k]!=0 ? 0x01 : 0x00
 end
 @synchronize
 for sweep in 1:sweeps,color in 0:1
  for q in lane:256:6859
   x=mod(q-1,19)+1;y=mod((q-1)÷19,19)+1;z=(q-1)÷361+1
   if 1<x<19 && 1<y<19 && 1<z<19 && fixed[q]==0 && mod(ci+cj+ck+x+y+z-30,2)==color
    values[q]=(values[q-1]+values[q+1]+values[q-19]+values[q+19]+values[q-361]+values[q+361])*(1f0/6f0)
   end
  end
  @synchronize
 end
 # Check every changed stencil (patch plus halo). A halo neighbor outside
 # shared memory is outside the written patch, so its global value is unchanged.
 largest=0f0
 for q in lane:256:6859
  x=mod(q-1,19)+1;y=mod((q-1)÷19,19)+1;z=(q-1)÷361+1
  i=ci+x-10;j=cj+y-10;k=ck+z-10
  if 1<i<n && 1<j<n && 1<k<n && fixed[q]==0
   xm=x>1 ? values[q-1] : u[i-1,j,k]
   xp=x<19 ? values[q+1] : u[i+1,j,k]
   ym=y>1 ? values[q-19] : u[i,j-1,k]
   yp=y<19 ? values[q+19] : u[i,j+1,k]
   zm=z>1 ? values[q-361] : u[i,j,k-1]
   zp=z<19 ? values[q+361] : u[i,j,k+1]
   largest=max(largest,abs(xm+xp+ym+yp+zm+zp-6f0*values[q]))
  end
 end
 maxima[lane]=largest
 @synchronize
 for offset in (128,64,32,16,8,4,2,1)
  if lane<=offset;maxima[lane]=max(maxima[lane],maxima[lane+offset]);end
  @synchronize
 end
 if lane==1;norm[1]=maxima[1];end
 for q in lane:256:6859
  x=mod(q-1,19)+1;y=mod((q-1)÷19,19)+1;z=(q-1)÷361+1
  i=ci+x-10;j=cj+y-10;k=ck+z-10
  if 1<x<19 && 1<y<19 && 1<z<19 && 1<i<n && 1<j<n && 1<k<n
   u[i,j,k]=values[q]
  end
 end
 end # inbounds: guarded physical/shared indices
end
function electrical_device_patch8!(u,mask,I;sweeps=16,norm=similar(u,Float32,1))
 size(u)==size(mask) && all(==(size(u,1)),size(u)) || throw(DimensionMismatch("patch requires matching cubic fields"))
 sweeps isa Integer && sweeps>=0 || throw(ArgumentError("nonnegative integer sweeps required"))
 n=size(u,1);all(d->1<I[d]<n,1:3) || throw(ArgumentError("patch center must be interior"))
 backend=KernelAbstractions.get_backend(u)
 KernelAbstractions.get_backend(mask)==backend || throw(ArgumentError("patch arrays must share backend"))
 length(norm)==1 && KernelAbstractions.get_backend(norm)==backend || throw(ArgumentError("patch norm workspace mismatch"))
 electrical_patch8_kernel!(backend,256)(u,mask,norm,Tuple(I)...,n,sweeps;ndrange=256)
 (residual=Float64(maximum(norm))/6,)
end
