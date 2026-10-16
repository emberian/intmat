# Node-centered geometric multigrid for ∇·(σ∇V)=0, x Dirichlet and y/z no-flux.
mutable struct LensLevel
 v::Array{Float64,3};sigma::Array{Float64,3};rhs::Array{Float64,3};r::Array{Float64,3};h::Float64
 gx::Array{Float64,3};gy::Array{Float64,3};gz::Array{Float64,3};diagonal::Array{Float64,3}
end
function lens_levels(n)
 levels=LensLevel[]
 while true
  push!(levels,LensLevel(zeros(n,n,n),ones(n,n,n),zeros(n,n,n),zeros(n,n,n),2/(n-1),zeros(n,n,n),zeros(n,n,n),zeros(n,n,n),zeros(n,n,n)))
  n<=5 && break
  isodd(n) || error("Nested nodal grids require2^k+1");n=(n+1)÷2
 end
 levels
end
@inline lens_face(a,b)=2a*b/(a+b)
function lens_coefficients!(l)
 n=size(l.v,1);s=l.sigma
 Threads.@threads :static for k in 1:n
  @inbounds for j in 1:n,i in 1:n
   l.gx[i,j,k]=i<n ? lens_face(s[i,j,k],s[i+1,j,k]) : 0.
   l.gy[i,j,k]=j<n ? lens_face(s[i,j,k],s[i,j+1,k]) : 0.
   l.gz[i,j,k]=k<n ? lens_face(s[i,j,k],s[i,j,k+1]) : 0.
  end
 end
 Threads.@threads :static for k in 1:n
  @inbounds for j in 1:n,i in 1:n
   l.diagonal[i,j,k]=l.gx[i,j,k]+l.gy[i,j,k]+l.gz[i,j,k]+(i>1 ? l.gx[i-1,j,k] : 0.)+(j>1 ? l.gy[i,j-1,k] : 0.)+(k>1 ? l.gz[i,j,k-1] : 0.)
  end
 end
end
@inline function lens_sum(l,i,j,k)
 n=size(l.v,1);v=l.v
 @inbounds num=l.gx[i-1,j,k]*v[i-1,j,k]+l.gx[i,j,k]*v[i+1,j,k]+l.gy[i,j,k]*v[i,min(n,j+1),k]+l.gz[i,j,k]*v[i,j,min(n,k+1)]+(j>1 ? l.gy[i,j-1,k]*v[i,j-1,k] : 0.)+(k>1 ? l.gz[i,j,k-1]*v[i,j,k-1] : 0.)
 num,l.diagonal[i,j,k]
end
function lens_smooth!(l,sweeps)
 n=size(l.v,1);h2=l.h^2
 for sweep in 1:sweeps,parity in 0:1
  Threads.@threads :static for k in 1:n
   @inbounds for j in 1:n
    first=2+mod(parity-j-k,2)
    @simd for i in first:2:n-1
     num,den=lens_sum(l,i,j,k);l.v[i,j,k]+=.95*((num+l.rhs[i,j,k]*h2)/den-l.v[i,j,k])
    end
   end
  end
 end
end
function lens_residual!(l)
 n=size(l.v,1);fill!(l.r,0);h2=l.h^2
 Threads.@threads :static for k in 1:n
  @inbounds for j in 1:n,i in 2:n-1
   num,den=lens_sum(l,i,j,k);l.r[i,j,k]=l.rhs[i,j,k]-(den*l.v[i,j,k]-num)/h2
  end
 end
 maximum(abs,l.r)
end
function lens_restrict!(coarse,fine)
 n=size(coarse,1);nf=size(fine,1)
 Threads.@threads :static for k in 1:n
  @inbounds for j in 1:n,i in 1:n
   x=2i-1;y=2j-1;z=2k-1;total=0.
   for dz in -1:1,dy in -1:1,dx in -1:1
    weight=(dx==0 ? 2 : 1)*(dy==0 ? 2 : 1)*(dz==0 ? 2 : 1)
    total+=weight*fine[clamp(x+dx,1,nf),clamp(y+dy,1,nf),clamp(z+dz,1,nf)]
   end
   coarse[i,j,k]=total/64
  end
 end
end
function lens_prolong_add!(fine,coarse)
 n=size(fine,1)
 Threads.@threads :static for k in 1:n
  @inbounds for j in 1:n,i in 2:n-1
   a=(i+1)÷2;b=(j+1)÷2;c=(k+1)÷2;tx=iseven(i) ? .5 : 0.;ty=iseven(j) ? .5 : 0.;tz=iseven(k) ? .5 : 0.;v=0.
   for dz in 0:1,dy in 0:1,dx in 0:1
    weight=(dx==0 ? 1-tx : tx)*(dy==0 ? 1-ty : ty)*(dz==0 ? 1-tz : tz)
    weight==0 && continue
    v+=weight*coarse[a+dx,b+dy,c+dz]
   end
   fine[i,j,k]+=v
  end
 end
end
function lens_cycle!(levels,q=1)
 l=levels[q]
 if q==length(levels);lens_smooth!(l,60);return;end
 lens_smooth!(l,3);lens_residual!(l);c=levels[q+1];lens_restrict!(c.rhs,l.r);c.rhs[1,:,:].=0;c.rhs[end,:,:].=0;fill!(c.v,0)
 lens_cycle!(levels,q+1);lens_prolong_add!(l.v,c.v);lens_smooth!(l,3)
end
function lens_solve!(levels,sigma;tol=1e-4,maxcycles=40)
 l=first(levels);l.sigma.=sigma
 for q in 2:length(levels);lens_restrict!(levels[q].sigma,levels[q-1].sigma);end
 for level in levels;lens_coefficients!(level);end
 residual=Inf;scale=max(maximum(sigma),eps())
 for cycle in 1:maxcycles
  lens_cycle!(levels);residual=lens_residual!(l)/scale
  residual<tol && return (residual=residual,cycles=cycle)
 end
 (residual=residual,cycles=maxcycles)
end
