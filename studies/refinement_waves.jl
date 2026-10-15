# Fixed-extent crossed scalar-wave refinement; uses the validated core propagator.
using Interior, TOML, SHA
using Interior: JSON3
function initial_arrays(n::Int,c1::NTuple{3,Float32},r1::Float32,a1::Float32,c2::NTuple{3,Float32},r2::Float32,a2::Float32,medium::Float32,speed::Float32,feedback::Float32,dt::Float32)
 q=zeros(Float32,n,n,n); mineral=similar(q);h=2f0/n
 @inbounds for k=1:n,j=1:n,i=1:n
  x=(i-.5f0)*h-1; y=(j-.5f0)*h-1;z=(k-.5f0)*h-1
  q[i,j,k]=a1*exp(-((x-c1[1])^2+(y-c1[2])^2+(z-c1[3])^2)/(2r1^2))+a2*exp(-((x-c2[1])^2+(y-c2[2])^2+(z-c2[3])^2)/(2r2^2))
  mineral[i,j,k]=medium*exp(-((z-.12f0*sin(5x+2y))/.055f0)^2)*exp(-((x/.78f0)^8+(y/.78f0)^8))
 end
 previous=copy(q)
 # Leapfrog half-acceleration gives zero initial velocity in the heterogeneous medium.
 @inbounds for k=1:n,j=1:n,i=1:n
  divergence=0f0
  for axis=1:3,side in (-1,1)
   v=mod1((i,j,k)[axis]+side,n); neighbor=axis==1 ? (v,j,k) : axis==2 ? (i,v,k) : (i,j,v)
   coeff=speed^2/(1+feedback*(mineral[i,j,k]+mineral[neighbor...])/2)
   divergence+=coeff*(q[neighbor...]-q[i,j,k])/h^2
  end
  previous[i,j,k]+=.5f0*dt^2*divergence
 end
 q,mineral,previous
end
function crossed_wave(config; output=nothing)
 c=TOML.parsefile(config); n=c["resolution"]; dt=Float32(c["dt"]); steps=c["steps"]
 output===nothing || !ispath(output) || error("refusing existing output")
 backend=get(c,"backend","cpu")=="metal" ? Metal.MetalBackend() : CPU()
 w=World(Grid(n,n,n);seed=114,backend)
 p=AcousticField(w;speed=c["speed"],frequency=c["frequency"],source_amplitude=c["source_amplitude"],source_radius=c["source_radius"],source_center=Tuple(c["source_center"]),material_feedback=c["material_feedback"],damping=c["damping"],activity_gain=.45,intensity_scale=.0007,envelope_time=.12)
 q,mineral,previous=initial_arrays(n,Tuple(Float32.(c["pulse1_center"])),Float32(c["pulse1_radius"]),Float32(c["pulse1_amplitude"]),Tuple(Float32.(c["pulse2_center"])),Float32(c["pulse2_radius"]),Float32(c["pulse2_amplitude"]),Float32(c["medium_amplitude"]),p.speed,p.material_feedback,dt)
 copyto!(w.fields.deposit,mineral);copyto!(w.fields.nutrient,1f0.-mineral)
 copyto!(p.pressure,q);copyto!(w.fields.pressure,q)
 copyto!(p.previous,previous);p.timestep=dt
 q=nothing;previous=nothing;mineral=nothing;GC.gc()
 procs=(p,Deposition(rate=.45,threshold=.18,dissolution=.0002))
 limit=minimum(stability_limit(w,x) for x in procs);dt<=limit || error("CFL $dt > $limit")
 println(JSON3.write((stage="initialized",shape=w.grid.size,limit=limit)));flush(stdout)
 elapsed=@elapsed for step=1:steps
  step!(w,procs,dt)
  if step%50==0;println(JSON3.write((stage="advance",step=step,time=w.time)));flush(stdout);end
 end
 finite=all(isfinite,p.pressure)&&all(isfinite,w.fields.deposit)
 finite || error("nonfinite")
 output===nothing && return (finite=finite,time=w.time,pressure=Array(p.pressure))
 mkpath(output);mkpath(joinpath(output,"native"));mkpath(joinpath(output,"volume"))
 for (name,a) in (("pressure",p.pressure),("previous",p.previous),("intensity",p.envelope),("activity",w.fields.activity),("deposit",w.fields.deposit),("nutrient",w.fields.nutrient))
  open(io->write(io,Array(a)),joinpath(output,"native",name*".f32"),"w")
 end
 pressure=Array(p.pressure); deposit=Array(w.fields.deposit)
 sample=sort(abs.(vec(pressure)[1:997:end])); scale=sample[ceil(Int,.995length(sample))]
 crest=clamp.(abs.(pressure)./max(scale,1f-8),0,1);mask=similar(crest)
 for I in CartesianIndices(mask)
  x,y,z=ntuple(d->2f0*(I[d]-.5f0)/n-1,3)
  mask[I]=clamp((1-max(abs(x)/.94,abs(y)/.94,abs(z)/.94))*40,0,1)
 end
 export_volume((density=.035f0.+.9f0.*deposit.+.38f0.*crest.^2,activity=crest.^2,structure=deposit,pressure=pressure,mask=mask),joinpath(output,"volume"))
 cp(config,joinpath(output,"recipe.toml"))
 report=(id=c["id"],title=c["title"],description=c["description"],shape=w.grid.size,extent=w.grid.extent,time=w.time,steps=steps,dt=dt,params=c,elapsed_seconds=elapsed,finite=finite,sourcehash=bytes2hex(sha256(read(@__FILE__))),source_sha256=Interior.arena_sources(),optical_scale=scale,optical_transfer="99.5th-percentile-clamped absolute pressure squared; artistic crest surrogate, no physical light emission",native_layout="Float32 x-fastest column-major; same grid and extent as volume",metrics=(pressure_extrema=extrema(p.pressure),deposit_extrema=extrema(w.fields.deposit),activity_extrema=extrema(w.fields.activity)))
 open(io->JSON3.pretty(io,report),joinpath(output,"run.json"),"w")
 job=(id=c["id"],family="waves",title=c["title"],description=c["description"],volume="volume",renderer="volume",size=2400,samples=128,args=["--palette","acoustic-pearl","--emission","0","--absorption",".6","--scattering",".12","--density-scale",".8","--polished-carrier","--carrier-shape","box","--orbit",".65","--elevation",".34"])
 open(io->JSON3.pretty(io,job),joinpath(output,"render-job.json"),"w")
 println("READY ",output);report
end
if abspath(PROGRAM_FILE)==@__FILE__
 get(TOML.parsefile(ARGS[1]),"backend","cpu")=="metal" && (@eval using Metal)
 Base.invokelatest(crossed_wave,ARGS[1];output=length(ARGS)>1 ? abspath(ARGS[2]) : nothing)
end
