# CPU-only research prototype; no production electrical methods are overwritten.
using Interior, LinearAlgebra
using Interior: JSON3
BLAS.set_num_threads(1)
function branch_grid_size(n)
 n isa Integer && n>=5 && ispow2(n-1) || throw(ArgumentError("multigrid requires n=2^k+1, n>=5 (for example65,129,257)"))
 Int(n)
end
struct MGLevel
 u::Array{Float32,3}; rhs::Array{Float32,3}; residual::Array{Float32,3}; occupied::BitArray{3}
end
function rb!(u,b,mask,sweeps;omega=1f0)
 n=size(u,1)
 for sweep in 1:sweeps, color in 0:1
  @inbounds for k in 2:n-1,j in 2:n-1
   first=2+mod(color-j-k,2)
   @simd for i in first:2:n-1
    if !mask[i,j,k]
     avg=(u[i-1,j,k]+u[i+1,j,k]+u[i,j-1,k]+u[i,j+1,k]+u[i,j,k-1]+u[i,j,k+1]+b[i,j,k])*(1f0/6f0)
     u[i,j,k]=muladd(omega,avg-u[i,j,k],u[i,j,k])
    end
   end
  end
 end
end
function residual!(r,u,b,mask)
 n=size(u,1);m=0f0
 @inbounds for k in 2:n-1,j in 2:n-1,i in 2:n-1
  val=mask[i,j,k] ? 0f0 : b[i,j,k]+u[i-1,j,k]+u[i+1,j,k]+u[i,j-1,k]+u[i,j+1,k]+u[i,j,k-1]+u[i,j,k+1]-6u[i,j,k]
  r[i,j,k]=val;m=max(m,abs(val)/6)
 end
 m
end
function hierarchy(u,mask)
 branch_grid_size(size(u,1))
 size(u)==size(mask) && all(==(size(u,1)),size(u)) || throw(DimensionMismatch("multigrid requires matching cubic fields"))
 levels=[MGLevel(u,zeros(Float32,size(u)),zeros(Float32,size(u)),mask)]
 while size(last(levels).u,1)>5
  f=last(levels);n=(size(f.u,1)+1)÷2;cm=falses(n,n,n)
  for I in findall(f.occupied)
   J=CartesianIndex(ntuple(d->clamp(round(Int,(I[d]+1)/2),2,n-1),3));cm[J]=true
  end
  push!(levels,MGLevel(zeros(Float32,n,n,n),zeros(Float32,n,n,n),zeros(Float32,n,n,n),cm))
 end
 levels
end
function vcycle!(levels,l=1)
 f=levels[l]
 if l==length(levels);rb!(f.u,f.rhs,f.occupied,40);return;end
 rb!(f.u,f.rhs,f.occupied,3)
 residual!(f.residual,f.u,f.rhs,f.occupied)
 c=levels[l+1];fill!(c.u,0);fill!(c.rhs,0);nc=size(c.u,1)
 @inbounds for k in 2:nc-1,j in 2:nc-1,i in 2:nc-1
  c.occupied[i,j,k] && continue
  s=0f0
  for dz in -1:1,dy in -1:1,dx in -1:1
   weight=(dx==0 ? 2 : 1)*(dy==0 ? 2 : 1)*(dz==0 ? 2 : 1)
   s+=weight*f.residual[2i-1+dx,2j-1+dy,2k-1+dz]
  end
  c.rhs[i,j,k]=s/16 # restriction average /64; coarse h²=4 fine h²
 end
 vcycle!(levels,l+1)
 n=size(f.u,1)
 @inbounds for k in 2:n-1,j in 2:n-1,i in 2:n-1
  f.occupied[i,j,k] && continue
  ci=(i+1)÷2;cj=(j+1)÷2;ck=(k+1)÷2
  tx=iseven(i) ? .5f0 : 0f0;ty=iseven(j) ? .5f0 : 0f0;tz=iseven(k) ? .5f0 : 0f0
  s=0f0
  for dz in 0:1,dy in 0:1,dx in 0:1
   s+=(dx==0 ? 1-tx : tx)*(dy==0 ? 1-ty : ty)*(dz==0 ? 1-tz : tz)*c.u[ci+dx,cj+dy,ck+dz]
  end
  f.u[i,j,k]+=s
 end
 rb!(f.u,f.rhs,f.occupied,3)
end
function mgsolve!(levels;tol=2e-5,maxcycles=200)
 r=Inf
 for it in 1:maxcycles
  vcycle!(levels);f=first(levels);r=residual!(f.residual,f.u,f.rhs,f.occupied)
  r<=tol && return (cycles=it,residual=r)
 end
 error("MG did not converge: $r")
end
function process(n,mask)
 p=Interior.DielectricGrowth(zeros(Float32,n,n,n),copy(mask),falses(n,n,n),UInt64(1),0,0.,.05,3.,.8f0,2e-5,5000,1.9,Inf,0,false)
 p.potential[mask].=1
 for I in findall(mask);Interior.electrical_mark!(p,I);end
 p
end
function specimen(n)
 m=falses(n,n,n);c=(n+1)÷2
 # A nontrivial thin conductor: trunk and three forked arms; fresh fine-grid geometry.
 for t in 0:round(Int,n*.25)
  m[c,c,c+t]=true
  if t>n÷12
   m[c+t-n÷12,c,c+n÷12]=true
   m[c-t+n÷12,c,c+n÷12]=true
   m[c,c+t-n÷12,c+n÷12]=true
  end
 end
 m
end
function benchmark(ns)
 results=[]
 # Compile on a tiny domain first.
 mgsolve!(hierarchy(process(17,specimen(17)).potential,specimen(17)))
 for n in ns
  mask=specimen(n);p=process(n,mask);u=copy(p.potential);levels=hierarchy(u,mask)
  stats=nothing;t=@elapsed stats=mgsolve!(levels)
  exactres=Interior.potential_residual(process_with_u(n,mask,u))
  row=Dict("n"=>n,"mg_seconds"=>t,"mg_cycles"=>stats.cycles,"residual"=>exactres,"occupied"=>count(mask))
  if n<=65
   ts=@elapsed Interior.solve_potential!(p)
   row["sor_seconds"]=ts;row["sor_iterations"]=p.last_iterations;row["maximum_potential_difference"]=maximum(abs.(p.potential.-u))
   # Tight solve separates solver accuracy from the normal stopping tolerance.
   tight=mgsolve!(levels;tol=3e-7,maxcycles=500);p.tolerance=3e-7;Interior.solve_potential!(p)
   row["tight_max_difference"]=maximum(abs.(p.potential.-u));row["tight_residual"]=tight.residual
  end
  # Exact fine-grid addition with warm hierarchy rebuilt for its new Dirichlet constraint.
  c=(n+1)÷2;J=CartesianIndex(c+1,c,c);mask[J]=true;u[J]=1
  levels=hierarchy(u,mask);tw=@elapsed warm=mgsolve!(levels)
  row["warm_add_seconds"]=tw;row["warm_add_cycles"]=warm.cycles;row["warm_add_residual"]=warm.residual
  push!(results,row);println(JSON3.write(row));flush(stdout)
 end
 mkpath(joinpath(@__DIR__,"..","outputs","refinement"))
 open(joinpath(@__DIR__,"..","outputs","refinement","branch-solver-benchmark.json"),"w") do io;JSON3.pretty(io,results);end
end
function process_with_u(n,mask,u)
 p=process(n,mask);p.potential.=u;p
end
if abspath(PROGRAM_FILE)==@__FILE__
 benchmark(isempty(ARGS) ? [33,65,129] : parse.(Int,ARGS))
end

# Stateless selection is identical to electrical_event!: same candidate order,
# same per-face eta power, same hash; only the harmonic solve is replaceable.
function select_site(p,spacing)
 candidates=CartesianIndex{3}[];weights=Float64[]
 for I in CartesianIndices(p.frontier)
  p.frontier[I] || continue
  weight=0.
  for (axis,h) in enumerate(spacing)
   offset=CartesianIndex(ntuple(d->d==axis ? 1 : 0,3))
   for J in (I-offset,I+offset)
    p.occupied[J] && (weight+=(max(0.,1-Float64(p.potential[I]))/h)^p.eta)
   end
  end
  push!(candidates,I);push!(weights,weight)
 end
 isempty(candidates) && error("no frontier")
 total=sum(weights);target=Interior.electrical_uniform(p.seed,p.events+1)*total;acc=0.;selected=last(candidates)
 for (I,weight) in zip(candidates,weights)
  acc+=weight
  if acc>target;selected=I;break;end
 end
 selected
end
function mark_hierarchy!(levels,I)
 for l in 2:length(levels)
  n=size(levels[l].u,1)
  I=CartesianIndex(ntuple(d->clamp(round(Int,(I[d]+1)/2),2,n-1),3))
  levels[l].occupied[I]=true
 end
end
function growth_parity(;n=33,events=100,tol=3e-7)
 mask=falses(n,n,n);root=CartesianIndex(ntuple(_->(n+1)÷2,3));mask[root]=true
 p=process(n,mask);q=process(n,mask);p.tolerance=tol;q.tolerance=tol;p.eta=q.eta=4.
 levels=hierarchy(q.potential,q.occupied)
 Interior.solve_potential!(p);mgsolve!(levels;tol=tol,maxcycles=500)
 matches=0;maxdiff=0.;maxres=0.
 elapsed=@elapsed for event in 1:events
  I=select_site(p,(1.,1.,1.));J=select_site(q,(1.,1.,1.))
  I==J || error("growth diverged at event$event: $I != $J")
  matches+=1
  Interior.electrical_mark!(p,I);p.events+=1;Interior.solve_potential!(p)
  Interior.electrical_mark!(q,J);q.events+=1;mark_hierarchy!(levels,J);mgsolve!(levels;tol=tol,maxcycles=500)
  maxdiff=max(maxdiff,maximum(abs.(p.potential.-q.potential)));maxres=max(maxres,Interior.potential_residual(q))
  p.reached_ground && break
 end
 row=(n=n,matched_events=matches,max_potential_difference=maxdiff,max_residual=maxres,seconds=elapsed,grounded=p.reached_ground)
 println(JSON3.write(row));row
end
