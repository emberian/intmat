# Distributed variable-conductivity material, with closed precursor/mineral transfer.
using Interior, JSON3, TOML, SHA, Random, LinearAlgebra, Statistics
BLAS.set_num_threads(1)
include(joinpath(@__DIR__,"weather_lens_solver.jl"))
include(joinpath(@__DIR__,"weather_lens_checks.jl"))
@inline face(a,b)=2a*b/(a+b)
function potential!(v,sigma,voltage,tol,maxsweeps)
 n=size(v,1);residual=Inf
 for sweep in 1:maxsweeps
  for parity in 0:1
   Threads.@threads :static for k in 1:n
    @inbounds for j in 1:n,i in 2:n-1
     (i+j+k)%2==parity || continue
     c=sigma[i,j,k];num=0.;den=0.
     for (x,y,z) in ((i-1,j,k),(i+1,j,k),(i,max(1,j-1),k),(i,min(n,j+1),k),(i,j,max(1,k-1)),(i,j,min(n,k+1)))
      g=face(c,sigma[x,y,z]);num+=g*v[x,y,z];den+=g
     end
     v[i,j,k]+=1.55*(num/den-v[i,j,k])
    end
   end
  end
  if sweep%10==0
   residual=0.
   @inbounds for k in 1:n,j in 1:n,i in 2:n-1
    c=sigma[i,j,k];num=0.;den=0.
    for (x,y,z) in ((i-1,j,k),(i+1,j,k),(i,max(1,j-1),k),(i,min(n,j+1),k),(i,j,max(1,k-1)),(i,j,min(n,k+1)))
     g=face(c,sigma[x,y,z]);num+=g*v[x,y,z];den+=g
    end
    residual=max(residual,abs(num/den-v[i,j,k]))
   end
   residual<tol && return (residual=residual,sweeps=sweep)
  end
 end
 (residual=residual,sweeps=maxsweeps)
end
function lens_joule!(out,v,sigma,h)
 n=size(v,1)
 Threads.@threads :static for k in 1:n
  @inbounds for j in 1:n,i in 1:n
   q=0.;a=sigma[i,j,k]
   for (x,y,z) in ((max(1,i-1),j,k),(min(n,i+1),j,k),(i,max(1,j-1),k),(i,min(n,j+1),k),(i,j,max(1,k-1)),(i,j,min(n,k+1)))
    q+=.5face(a,sigma[x,y,z])*((v[i,j,k]-v[x,y,z])/h)^2
   end
   out[i,j,k]=q
  end
 end
end
function lens_material_step!(temp,scratch,mineral,precursor,damage,joule,spacing,dt,p)
 n=size(temp,1)
 Threads.@threads :static for k in 1:n
  @inbounds for j in 1:n,i in 1:n
   scratch[i,j,k]=max(0.,temp[i,j,k]+dt*(p.D*Interior.laplacian(temp,i,j,k,spacing)-p.cooling*temp[i,j,k]+p.jgain*joule[i,j,k]))
  end
 end
 deposited=0.;eroded=0.
 @inbounds for q in eachindex(temp)
  T=scratch[q];f=precursor[q];m=mineral[q]
  deposit=f*(-expm1(-dt*p.dgain*max(0.,T-p.dthreshold)));erode=m*(-expm1(-dt*p.egain*max(0.,T-p.ethreshold)))
  precursor[q]=f-deposit+erode;mineral[q]=m+deposit-erode;deposited+=deposit;eroded+=erode
  damage[q]=clamp(damage[q]+dt*(p.damagegain*max(0.,T-p.damagethreshold)*(1-damage[q])-p.recovery*damage[q]),0,1)
 end
 deposited,eroded
end
function main(recipe,out)
 c=TOML.parsefile(recipe);n=Int(c["n"]);dt=Float64(c["dt"]);nodal=get(c,"nodal_grid",false);h=nodal ? 2/(n-1) : 2/n;spacing=(h,h,h);rng=MersenneTwister(c["seed"])
 steps=Int(c["steps"]);cadence=Int(c["capture_every"]);voltage=Float64(c["voltage"]);tol=Float64(c["potential_tolerance"]);maxsweeps=Int(c["max_sweeps"])
 D=Float64(c["thermal_diffusion"]);cooling=Float64(c["cooling"]);jgain=Float64(c["joule_gain"]);cgain=Float64(c["conductivity_gain"])
 dgain=Float64(c["deposition_gain"]);egain=Float64(c["erosion_gain"]);dthreshold=Float64(c["deposition_threshold"]);ethreshold=Float64(c["erosion_threshold"])
 damagegain=Float64(c["damage_gain"]);damagethreshold=Float64(c["damage_threshold"]);recovery=Float64(c["damage_recovery"])
 thermal=(;D,cooling,jgain,dgain,egain,dthreshold,ethreshold,damagegain,damagethreshold,recovery)
 frozen=get(c,"immutable_conductivity",false);continuous=get(c,"continuous_texture",false)
 dt*(6D/h^2+cooling)<.9 || error("Thermal explicit step violates positivity bound")
 ispath(out) && error("Fresh output required");mkpath(out);cp(recipe,joinpath(out,"recipe.toml"));write(joinpath(out,"weather_lens.jl"),read(@__FILE__));write(joinpath(out,"weather_lens_solver.jl"),read(joinpath(@__DIR__,"weather_lens_solver.jl")));write(joinpath(out,"weather_lens_checks.jl"),read(joinpath(@__DIR__,"weather_lens_checks.jl")))
 v=zeros(n,n,n);mineral=zeros(n,n,n);precursor=zeros(n,n,n);temp=zeros(n,n,n);scratch=similar(temp);damage=zeros(n,n,n);sigma=similar(temp);joule=zeros(n,n,n);host=zeros(n,n,n)
 @inbounds for k in 1:n,j in 1:n,i in 1:n
  x,y,z=nodal ? ((i-1)*h-1,(j-1)*h-1,(k-1)*h-1) : (2(i-.5)/n-1,2(j-.5)/n-1,2(k-.5)/n-1)
  fold=.25sin(3.8x+1.1sin(3z))+.08sin(6z)
  lamella=sum(exp(-((y-fold-b)/.065)^2) for b in (-.42,0.,.42))
  host[i,j,k]=min(1.,lamella)*exp(-max(0.,abs(z)-.78)^2/.02)
  precursor[i,j,k]=host[i,j,k];mineral[i,j,k]=.025host[i,j,k]*(continuous ? 1+.25sin(17x+11z+7y)*sin(13z-9x) : .75+.5rand(rng));precursor[i,j,k]-=mineral[i,j,k]
  v[i,j,k]=voltage*(n-i)/(n-1)
 end
 levels=nodal ? lens_levels(n) : nothing
 nodal && (first(levels).v.=v)
 immutable_sigma=similar(sigma)
 initial=sum(precursor)+sum(mineral);history=Any[];frames=Any[];deposited=0.;eroded=0.;dissipation=0.;misses=0
 function capture(step,solved)
  path=joinpath(out,"frames",lpad(step,5,'0'));mkpath(path)
  export_volume((density=copy(host),structure=copy(mineral),activity=copy(joule),deposit=copy(mineral),temperature=copy(temp),damage=copy(damage),potential=copy(v),mask=ones(Float32,n,n,n)),path)
  if nodal
   manifest=JSON3.read(read(joinpath(path,"volume.json"),String),Dict{String,Any});manifest["sample_convention"]="nodes including physical endpoints";manifest["spacing"]=[h,h,h];open(io->JSON3.write(io,manifest),joinpath(path,"volume.json"),"w")
  end
  input_power=sum(face(sigma[1,j,k],sigma[2,j,k])*(v[1,j,k]-v[2,j,k])/h*voltage*h^2 for k in 1:n,j in 1:n)
  bulk_power=sum(joule)*h^3
  checks=lens_edge_checks(v,sigma,h,voltage)
  row=(step=step,time=step*dt,independent=checks,input_power=input_power,power_relative_defect=(bulk_power-input_power)/max(input_power,eps()),temperature_max=maximum(temp),mineral=sum(mineral),damage_mean=mean(damage),joule_power=sum(joule)*h^3,active_voxels=count(>(.02),joule),residual=solved.residual,sweeps=solved.sweeps,mass_residual=sum(precursor)+sum(mineral)-initial)
  push!(history,row);push!(frames,(step=step,time=step*dt,directory=relpath(path,out)));println(JSON3.write(row));flush(stdout)
 end
 cached_solve=(residual=Inf,sweeps=0)
 elapsed=@elapsed for step in 0:steps
  if !frozen || step==0
   @. sigma=.008+host*(.03+cgain*mineral^2)*exp(-6damage)
   step==0 && (immutable_sigma.=sigma)
  else
   sigma.=immutable_sigma
  end
  if frozen && step>0
   solved=cached_solve
  elseif nodal
   mg=lens_solve!(levels,sigma;tol=tol,maxcycles=maxsweeps);v.=first(levels).v;solved=(residual=mg.residual,sweeps=mg.cycles)
  else
   solved=potential!(v,sigma,voltage,tol,maxsweeps)
  end
  solved.residual>=tol && (misses+=1)
  cached_solve=solved
  (!frozen || step==0) && lens_joule!(joule,v,sigma,h)
  step%cadence==0 && capture(step,solved)
  step==steps && break
  dissipation+=sum(joule)*h^3*dt
  # Existing Interior periodic thermal Laplacian; Joule power replaces a toolpath.
  newdeposit,newerosion=lens_material_step!(temp,scratch,mineral,precursor,damage,joule,spacing,dt,thermal)
  deposited+=newdeposit;eroded+=newerosion
  temp,scratch=scratch,temp
 end
 all(a->all(isfinite,a),(v,mineral,precursor,temp,damage,sigma,joule)) || error("Nonfinite state")
 minimum(precursor)>=0 && minimum(mineral)>=0 || error("Negative material")
 minimum(v)>=-tol && maximum(v)<=voltage+tol || error("Potential maximum principle failed")
 misses==0 || error("Unconverged electrical states: $misses")
 report=(title=c["title"],recipe=c,sourcehash=bytes2hex(sha256(read(joinpath(out,"weather_lens.jl")))),solverhash=bytes2hex(sha256(read(joinpath(out,"weather_lens_solver.jl")))),frames=frames,history=history,elapsed=elapsed,thermal_positivity_factor=dt*(6D/h^2+cooling),maximum_local_material_error=maximum(abs,precursor.+mineral.-host),initial_material=initial,final_material=sum(precursor)+sum(mineral),total_deposition=deposited,total_erosion=eroded,electrical_dissipation=dissipation,unconverged_solves=misses,description="Designed folded precursor lamellae and fixed opposing electrodes. Variable-conductivity potential uses nodal multigrid or the pilot weighted red-black SOR. Thermal diffusion is periodic; electrical boundaries are fixed x electrodes and no-flux y/z. Joule heating deposits conductive mineral at moderate temperature; high temperature erodes it and creates recovering insulating damage. Precursor and mineral exchange conservatively. Electrical supply/cooling and damage law are open and phenomenological. Temperature is nondimensional; conversion does not debit latent heat or chemical free energy.")
 open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w")
end
main(ARGS[1],ARGS[2])
