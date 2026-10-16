include(joinpath(@__DIR__,"weather_lens_solver.jl"))
using JSON3
function check(n,heterogeneous)
 levels=lens_levels(n);l=first(levels);exact=similar(l.v)
 for k in 1:n,j in 1:n,i in 1:n
  x,y,z=(i-1)/(n-1),(j-1)/(n-1),(k-1)/(n-1)
  l.sigma[i,j,k]=heterogeneous ? exp(2sin(3x+2y)*cos(4z)) : 1.
  exact[i,j,k]=1-x+.05sin(pi*x)*cos(2pi*y)*cos(2pi*z)
 end
 l.v.=exact;lens_coefficients!(l);lens_residual!(l);l.rhs.=-l.r;l.v.=0;l.v[1,:,:].=exact[1,:,:];l.v[end,:,:].=exact[end,:,:]
 solved=lens_solve!(levels,copy(l.sigma);tol=1e-6,maxcycles=60);error=maximum(abs,l.v.-exact)
 println(JSON3.write((n=n,heterogeneous=heterogeneous,solved=solved,error=error)))
 @assert solved.residual<1e-6 && error<1e-5
end
for n in (17,33),heterogeneous in (false,true);check(n,heterogeneous);end
include(joinpath(@__DIR__,"weather_lens_checks.jl"))
let n=17
 l=first(lens_levels(n));sigma=ones(n,n,n)
 for i in 1:n;l.v[i,:,:].=1-(i-1)/(n-1);end
 c=lens_edge_checks(l.v,sigma,l.h,1.)
 println(JSON3.write((independent_constant_field=c,)))
 @assert c.independent_residual<1e-12 && abs(c.independent_power_defect)<1e-12
end
