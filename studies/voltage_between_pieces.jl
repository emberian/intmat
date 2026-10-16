# Optical contours of a genuinely rupturing 3D conductive spring material.
using Interior, LinearAlgebra, SHA
using Interior: JSON3
BLAS.set_num_threads(1)
isdefined(Interior,:BondNetwork)||Base.include(Interior,joinpath(@__DIR__,"../src/mechanics.jl"))
function graph_apply!(out,x,p,g,den,pinned)
 out.=den.*x
 @inbounds for e in eachindex(g)
  g[e]==0&&continue;i=p.bonds[1,e];j=p.bonds[2,e]
  !pinned[i]&&!pinned[j]&&(out[i]-=g[e]*x[j];out[j]-=g[e]*x[i])
 end
 out[pinned].=0
end
function voltage!(v,p,g,den,r,z,d,q,b,pinned;tol=1e-7,maxiter=1000)
 fill!(den,0);fill!(b,0);pinned.=p.grip.!=0
 @inbounds for e in eachindex(g)
  i=p.bonds[1,e];j=p.bonds[2,e];den[i]+=g[e];den[j]+=g[e]
 end
 # Components without any electrode receive a constant arbitrary gauge; no Joule heat.
 adjacency=[Int[] for _ in eachindex(v)]
 @inbounds for e in eachindex(g)
  g[e]>0||continue;i=p.bonds[1,e];j=p.bonds[2,e];push!(adjacency[i],j);push!(adjacency[j],i)
 end
 seen=falses(length(v));queue=Int[];component=Int[];floating=0
 for root in eachindex(v)
  seen[root]&&continue;empty!(queue);empty!(component);push!(queue,root);seen[root]=true;touch=false;head=1
  while head<=length(queue)
   i=queue[head];head+=1;push!(component,i);touch|=p.grip[i]!=0
   for j in adjacency[i];!seen[j]&&(seen[j]=true;push!(queue,j));end
  end
  if !touch;floating+=1;pinned[component].=true;v[component].=0;end
 end
 for i in eachindex(v);p.grip[i]!=0&&(v[i]=(p.grip[i]+1)/2);end
 @inbounds for e in eachindex(g)
  i=p.bonds[1,e];j=p.bonds[2,e]
  if pinned[j]&&!pinned[i];b[i]+=g[e]*v[j];elseif pinned[i]&&!pinned[j];b[j]+=g[e]*v[i];end
 end
 # Solve only free potentials; fixed values enter the right-hand side.
 saved=copy(v);v[pinned].=0;graph_apply!(q,v,p,g,den,pinned);r.=b.-q;r[pinned].=0
 z.=r./max.(den,1e-30);d.=z;rz=dot(r,z);normb=max(norm(b),1e-30);iterations=0
 for k=1:maxiter
  norm(r)<=tol*normb&&break
  graph_apply!(q,d,p,g,den,pinned);dq=dot(d,q);dq>0||error("CG breakdown $dq");alpha=rz/dq;v.+=alpha.*d;r.-=alpha.*q;z.=r./max.(den,1e-30);next=dot(r,z);d.=z.+(next/rz).*d;rz=next;iterations=k
 end
 v[pinned].=saved[pinned]
 # Independent edge-current residual and source/Joule identity.
 fill!(q,0);joule=0.
 @inbounds for e in eachindex(g)
  i=p.bonds[1,e];j=p.bonds[2,e];current=g[e]*(v[i]-v[j]);q[i]+=current;q[j]-=current;joule+=current*(v[i]-v[j])
 end
 # All-fixed/all-floating graphs have no free equations and zero free residual.
 residual=maximum((abs(q[i])/max(den[i],1e-30) for i in eachindex(v) if !pinned[i]);init=0.);power=sum((v[i]*q[i] for i in eachindex(v) if p.grip[i]!=0);init=0.)
 residual<2e-6||error("independent graph residual $residual")
 (iterations=iterations,residual=residual,joule_power=joule,electrode_power=power,relative_power_defect=abs(joule-power)/max(joule,1e-30),floating_components=floating)
end
function produce(out;dims=(24,18,10),steps=1300,weakening=80.,scaled=false,snapshot=nothing,electrical_cadence=nothing)
 ispath(out)&&error("fresh output required");mkpath(out)
 ratio=scaled ? prod(((1.7/(dims[1]-1))/(1.7/23),(1.5/(dims[2]-1))/(1.5/17),(.95/(dims[3]-1))/(.95/9)))^(1/3) : 1.
 duration=steps*.0007;steps=ceil(Int,steps/ratio);dt=duration/steps;cadence=electrical_cadence===nothing ? max(1,round(Int,25/ratio)) : Int(electrical_cadence);cadence>0||error("electrical cadence must be positive")
 w=World(Grid(4,4,4);seed=1944);p=Interior.BondNetwork(w;size=dims,extent=(1.7,1.5,.95),mass=.008*ratio^3,stiffness=6*ratio,damping=5,failure_strain=.18,weak_seam=.42,stretch_rate=.26,shear_rate=.19,jitter=.14)
 n=length(p.node_ids);m=length(p.bond_ids);v=zeros(n);dose=zeros(n);temperature=zeros(n);base=copy(p.threshold);g=zeros(m);g0=zeros(m);den=zeros(n);r=zeros(n);z=zeros(n);d=zeros(n);q=zeros(n);b=zeros(n);pinned=falses(n);heat=zeros(n)
 # Quenched conductivity grains and two small weak notches are designed inputs.
 for e in 1:m
  i,j=p.bonds[:,e];x,y,zz=(p.reference[:,i]+p.reference[:,j])/2
  g0[e]=ratio^2*(.24+.76*(.5+.5cos(9y+4zz+2sin(6x)))^2)/p.rest_length[e]
  base[e]*=.72+.38*(.5+.5sin(11y+8zz-3x))
  if abs(x-.10sin(5y+4zz))<.04 && (y<-.32 || y>.38&&zz>.1);base[e]*=.45;end
  p.threshold[e]=base[e]
 end
 for i=1:n;x,y,zz=p.reference[:,i];v[i]=x/1.7+.5;p.load[3,i]=ratio^3*.012sin(5y+7zz)*exp(-4x*x);end
 cp(@__FILE__,joinpath(out,"source-study.jl"));history=Any[];lastsolve=1;initialenergy=Interior.mechanical_energy(p).total;solved=nothing
 elapsed=@elapsed for step=1:steps
  if step==1||step%cadence==0
   g.=g0.*(1 .-p.damage);solved=voltage!(v,p,g,den,r,z,d,q,b,pinned);fill!(heat,0)
   for e=1:m;i=p.bonds[1,e];j=p.bonds[2,e];power=g[e]*(v[i]-v[j])^2;heat[i]+=.5power/ratio^3;heat[j]+=.5power/ratio^3;end
   interval=step==1 ? 0. : (step-lastsolve)*dt;lastsolve=step
   for i=1:n;temperature[i]=exp(-interval*2)*temperature[i]+(-expm1(-interval*2))*heat[i];dose[i]+=interval*temperature[i];end
   for e=1:m;i=p.bonds[1,e];j=p.bonds[2,e];p.threshold[e]=base[e]/(1+weakening*(dose[i]+dose[j]));end
   snapshot===nothing||snapshot((stage=:electrical,step=step-1,time=w.time,process=p,potential=v,dose=dose,temperature=temperature,solver=solved))
  end
  step!(w,p,dt)
  if step%100==0||step==steps
   push!(history,(step=step,broken=count(==(1),p.damage),maxdose=maximum(dose),solver=solved));println(last(history));flush(stdout)
  end
 end
 g.=g0.*(1 .-p.damage);solved=voltage!(v,p,g,den,r,z,d,q,b,pinned)
 snapshot===nothing||snapshot((stage=:final,step=steps,time=w.time,process=p,potential=v,dose=dose,temperature=temperature,solver=solved))
 energy=Interior.mechanical_energy(p);energy_residual=energy.total-initialenergy-energy.external_work
 open(io->JSON3.write(io,(dims=dims,position=p.position,reference=p.reference,bonds=p.bonds,damage=p.damage,strain=p.peak_strain,potential=v,dose=dose,temperature=temperature)),joinpath(out,"state.json"),"w")
 report=(title="The voltage between the pieces",family="branches",shape=dims,nodes=n,bonds=m,steps=steps,dt=dt,weakening=weakening,pitch_ratio=ratio,physical_duration=duration,electrical_cadence=cadence,scaling=scaled,history=history,broken=count(==(1),p.damage),solver=solved,elapsed_seconds=elapsed,energy=energy,mechanical_energy_residual=energy_residual,source_sha256=Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,"..",f)))) for f in ("studies/voltage_between_pieces.jl","src/mechanics.jl")),model="Shared BondNetwork; graph conductivity vanishes on failed bonds. Independently checked PCG potential; local Joule dose weakens tensile thresholds. External electrodes/grips power exploratory law; not calibrated electromechanics.",authored="Quenched geometry/conductivity grains, two small weak notches, distributed load and grip motion. Optical contours are computed potential observations, not added matter.")
 open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w")
end
if abspath(PROGRAM_FILE)==@__FILE__
 produce(ARGS[1];dims=Tuple(parse.(Int,split(get(ENV,"PIECES_DIMS","24,18,10"),','))),steps=parse(Int,get(ENV,"PIECES_STEPS","1300")),weakening=parse(Float64,get(ENV,"PIECES_WEAKENING","80")),scaled=get(ENV,"PIECES_SCALED","0")=="1")
end
