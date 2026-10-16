using Interior, JSON3, LinearAlgebra, SHA
BLAS.set_num_threads(1)
isdefined(Interior,:BondNetwork) || Base.include(Interior,joinpath(@__DIR__,"..","src","mechanics.jl"))
function main()
 out=get(ENV,"INTMAT_FAULT_OUTPUT",joinpath(@__DIR__,"..","outputs","fault-remembers"));mkpath(out)
 capture=get(ENV,"INTMAT_FAULT_CAPTURE","0")=="1"
 frame_steps=Set(round.(Int,range(0,1450,length=40)))
 capture && mkpath(joinpath(out,"states"))
 function snapshot(k,dims,p,V,T,dose)
  open(io->JSON3.write(io,(step=k,time=k*.0007,dims=dims,position=p.position,reference=p.reference,bonds=p.bonds,damage=p.damage,strain=p.peak_strain,dose=dose,temperature=T,potential=V)),joinpath(out,"states",lpad(k,4,'0')*".json"),"w")
 end
 dims=(30,27,7); w=World(Grid(8,8,8);seed=7126)
 p=Interior.BondNetwork(w;size=dims,extent=(1.65,1.7,.65),mass=.008,stiffness=6.,damping=5.,failure_strain=.20,weak_seam=.7,stretch_rate=.28,shear_rate=.21,jitter=.17)
 n=length(p.node_ids); m=length(p.bond_ids); V=zeros(n); T=zeros(n); dose=zeros(n); num=zeros(n); den=zeros(n); heat=zeros(n); cond=zeros(m)
 # Lamination is mechanical and electrical anisotropy, present before loading.
 for e in 1:m
  i,j=p.bonds[:,e]; a=p.reference[:,i]; b=p.reference[:,j]; z=(a[3]+b[3])/2
  cross=abs(a[3]-b[3])>.04
  p.stiffness[e]*=cross ? .20 : 1.
  p.threshold[e]*=cross ? .50 : .72+.30*(.5+.5*sin(13(a[2]+b[2])+19z))
  # Three offset flaws seed competing fronts rather than prescribe open geometry.
  x=(a[1]+b[1])/2;y=(a[2]+b[2])/2
  if abs(x-.16*sin(5y+8z))<.045 && (y<-.12 || z>.1 && y>.25);p.threshold[e]*=.46;end
 end
 for i in 1:n
  x,y,z=p.reference[:,i];p.load[3,i]=.014*sin(4y+9z)*exp(-5x*x)
 end
 base=copy(p.threshold);for i in 1:n; V[i]=(p.reference[1,i]/1.65+.5);end
 capture && snapshot(0,dims,p,V,T,dose)
 dt=.0007; steps=1450
 elapsed=@elapsed for k in 1:steps
  if k%10==0
   fill!(den,0)
   for e in 1:m
    i=p.bonds[1,e];j=p.bonds[2,e];cond[e]=p.damage[e]==1 ? 0. : (abs(p.reference[3,i]-p.reference[3,j])>.04 ? .045 : 1.)/p.rest_length[e]
    den[i]+=cond[e];den[j]+=cond[e]
   end
   for sweep in 1:24
    fill!(num,0)
    for e in 1:m;i=p.bonds[1,e];j=p.bonds[2,e];g=cond[e];num[i]+=g*V[j];num[j]+=g*V[i];end
    for i in 1:n; V[i]=p.grip[i]!=0 ? (p.grip[i]+1)/2 : den[i]>0 ? .72V[i]+.28num[i]/den[i] : V[i];end
   end
   fill!(heat,0)
   for e in 1:m;i=p.bonds[1,e];j=p.bonds[2,e];q=cond[e]*(V[i]-V[j])^2;heat[i]+=q/2;heat[j]+=q/2;end
   for i in 1:n;T[i]=.985T[i]+.015heat[i];dose[i]+=.007T[i];end
   for e in 1:m;i=p.bonds[1,e];j=p.bonds[2,e];p.threshold[e]=base[e]/(1+8(dose[i]+dose[j]));end
  end
  step!(w,p,dt)
  capture && k in frame_steps && snapshot(k,dims,p,V,T,dose)
  k%250==0 && println(k," broken ",count(==(1),p.damage)," dose ",maximum(dose))
 end
 Interior.validate(w,p)
 open(io->JSON3.write(io,(dims=dims,position=p.position,reference=p.reference,bonds=p.bonds,damage=p.damage,strain=p.peak_strain,dose=dose,temperature=T,potential=V)),joinpath(out,"state.json"),"w")
 open(io->JSON3.pretty(io,(title="The fault remembers",steps=steps,dt=dt,seed=7126,nodes=n,bonds=m,broken=count(==(1),p.damage),maxdose=maximum(dose),elapsed=elapsed,energy=Interior.mechanical_energy(p),sourcehash=bytes2hex(sha256(read(@__FILE__))),description="A seven-lamina 3D spring material driven in tension and shear. Surviving bonds conduct a relaxed voltage field; Joule dissipation accumulates dose and lowers tensile thresholds. Broken bonds stop conducting. Membrane faces exist only where their perimeter bonds survive. This is an exploratory art law, not calibrated electromechanics.")),joinpath(out,"run.json"),"w")
end
main()
