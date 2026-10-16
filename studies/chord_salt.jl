# Wave-written finite material: authored packets/support, computed interference memory.
using Interior, JSON3, SHA, TOML, Statistics
function chord_run(out;n=64,steps=900,dt=.006,feedback=12.,amplitude=1.)
 ispath(out) && error("fresh output required");mkpath(out)
 w=World(Grid(n,n,n);seed=144);p=AcousticField(w;speed=.44,material_feedback=feedback,damping=.065,source_amplitude=0,activity_gain=0,envelope_time=.18)
 h=2/n;dt<Interior.stability_limit(w,p)||error("CFL")
 host=zeros(Float32,n,n,n);initial=zeros(Float32,n,n,n);dose=zeros(Float32,n,n,n)
 for k=1:n,j=1:n,i=1:n
  x,y,z=((i-.5)*h-1,(j-.5)*h-1,(k-.5)*h-1)
  # Ellipsoidal lens and open diagonal cut are intentional composition.
  r=sqrt((x/.86)^2+(y/.63)^2+(z/.83)^2)
  H=clamp((1-r)/.065,0,1)*clamp((y+.43+.33x-.13z)/.05,0,1)
  host[i,j,k]=H;initial[i,j,k]=.015H
  q=0.
  for (u,v,t,c,phase) in ((1.,.32,.22,-.34,0.),(-.27,1.,.41,.29,1.1),(.35,-.24,1.,-.18,2.2))
   d=(u*x+v*y+t*z)/sqrt(u*u+v*v+t*t)
   q+=amplitude*.7*cos(34d+phase)*exp(-((d-c)/.19)^2)*exp(-((x/.86)^8+(y/.86)^8+(z/.86)^8))
  end
  p.pressure[i,j,k]=q;p.previous[i,j,k]=q
 end
 # Zero velocity initial data must use leapfrog half acceleration, same symmetric operator.
 for k=1:n,j=1:n,i=1:n
  lap=0.;u=p.pressure[i,j,k]
  for (ii,jj,kk) in ((mod1(i-1,n),j,k),(mod1(i+1,n),j,k),(i,mod1(j-1,n),k),(i,mod1(j+1,n),k),(i,j,mod1(k-1,n)),(i,j,mod1(k+1,n)))
   c2=.44^2/(1+feedback*(initial[i,j,k]+initial[ii,jj,kk])/2);lap+=c2*(p.pressure[ii,jj,kk]-u)/h^2
  end
  p.previous[i,j,k]+=dt^2*.5lap
 end
 w.fields.deposit.=initial;w.fields.nutrient.=host.-initial
 history=Any[];elapsed=@elapsed for step=1:steps
  Interior.advance!(w,p,dt);w.time+=dt
  # Exact exponential transfer closes each voxel's finite material inventory.
  @inbounds for q in eachindex(host)
   intensity=p.envelope[q];dose[q]+=dt*intensity
   transfer=w.fields.nutrient[q]*(-expm1(-dt*5*max(0.,intensity-.018)))
   w.fields.deposit[q]+=transfer;w.fields.nutrient[q]=host[q]-w.fields.deposit[q]
  end
  if step%100==0||step==steps
   push!(history,(step=step,time=w.time,peak_pressure=maximum(abs,p.pressure),converted=sum(w.fields.deposit.-initial),peak_envelope=maximum(p.envelope)))
   println(last(history));flush(stdout)
  end
 end
 localerr=maximum(abs.(w.fields.deposit.+w.fields.nutrient.-host));all(isfinite,p.pressure)||error("nonfinite")
 export_volume((density=host,deposit=Array(w.fields.deposit),structure=Array(w.fields.deposit),activity=dose,pressure=Array(p.pressure),envelope=Array(p.envelope),mask=host),joinpath(out,"volume"))
 meta=(id="chord-salt-$n",title="A chord held in salt",shape=(n,n,n),steps=steps,dt=dt,time=w.time,feedback=feedback,amplitude=amplitude,elapsed_seconds=elapsed,local_inventory_error=localerr,history=history,source_sha256=bytes2hex(sha256(read(@__FILE__))),acoustic_sha256=bytes2hex(sha256(read(joinpath(@__DIR__,"../src/acoustics.jl")))),model="Existing periodic symmetric-face AcousticField; finite local exponential precursor-to-mineral conversion. Wave speed decreases with mineral; externally prepared packets, no ongoing forcing; conversion has no latent-energy debit.",authored="Three angled Gaussian-enveloped coherent packets; ellipsoidal support and open diagonal cut. Internal converted fringes computed.")
 open(io->JSON3.pretty(io,meta),joinpath(out,"run.json"),"w")
 meta
end
if abspath(PROGRAM_FILE)==@__FILE__
 chord_run(ARGS[1];n=parse(Int,get(ENV,"CHORD_N","64")),steps=parse(Int,get(ENV,"CHORD_STEPS","900")),dt=parse(Float64,get(ENV,"CHORD_DT",".006")),feedback=parse(Float64,get(ENV,"CHORD_FEEDBACK","12")),amplitude=parse(Float64,get(ENV,"CHORD_AMPLITUDE","1")))
end
