using Interior, SHA, TOML
using Interior: JSON3
Interior.record_process_fields(p::AcousticField) = NamedTuple{fieldnames(typeof(p))}(Tuple(getfield(p,n) for n in fieldnames(typeof(p))))
const ROOT=normpath(joinpath(@__DIR__,".."))
const SPECS=[
("waves-crossed-choir","Crossed choir","Unequal, opposite-sign pressure pulses cross a continuously driven off-axis transducer.","cross",240,"glass-ice"),
("waves-pearl-cavity","Pearl cavity","A pulsed scalar wave refracts around an off-center mineral shell; pressure crests are an emissive optical surrogate.","cavity",300,"acoustic-pearl"),
("waves-memory-rings","Memory rings","Concentric pressure bands write persistent mineral rings while their instantaneous crests move onward.","rings",260,"amber-gold"),
("waves-oblique-fronts","Oblique fronts","Three tilted planar pressure packets collide through a stratified mineral medium.","fronts",210,"glass-ice"),
("waves-neural-scroll","Neural scroll","A fixed recurrent two-channel tanh field unfolds a twisted phase initialization into mineral-writing excitation.","scroll",300,"acoustic-pearl"),
("waves-neural-lattice","Neural lantern lattice","Structured recurrent neural excitation competes with a mineral lattice, leaving a porous activity memory.","lattice",280,"amber-gold")]
function runone(spec;n=96,steps_override=0)
 slug,title,description,mode,steps,palette=spec; steps_override>0 && (steps=steps_override)
 out=joinpath(ROOT,"outputs","collection-20261015",slug); mkpath(out)
 w=World(Grid(n,n,n);seed=114); dt=.01f0
 p=AcousticField(w;speed=.44,frequency=mode=="rings" ? 2.2 : 1.5,source_amplitude=mode=="cross" ? 55 : 0,source_radius=.055,source_center=(-.38,.13,-.12),material_feedback=12,damping=.035,activity_gain=.45,intensity_scale=.0007,envelope_time=.12)
 neural=mode in ("scroll","lattice")
 nf=neural ? NeuralField(w;initial_noise=0,diffusion=.00045,response_time=.4,nutrient_feedback=.02,material_feedback=1.2,recurrence=(1.8,-1.6,1.4,1.7),gain=1) : nothing
 for k=1:n,j=1:n,i=1:n
  x,y,z=ntuple(d->2f0*((i,j,k)[d]-.5f0)/n-1,3); r=sqrt(x*x+y*y+z*z)
  g(a,b,c,s)=exp(-((x-a)^2+(y-b)^2+(z-c)^2)/(2s^2))
  q=mode=="cross" ? 1.8g(-.43,-.10,.08,.075)-1.3g(.35,.17,-.17,.10) : mode=="cavity" ? 2g(-.18,.08,.0,.07) : mode=="rings" ? .8cos(29r)*exp(-(r/.50)^6) : mode=="fronts" ? .85exp(-((x+.38+.3y)/.06)^2)*exp(-((y/.75)^6+(z/.7)^6))-.65exp(-((x-.3-.35z)/.07)^2)*exp(-((y/.7)^6+(z/.65)^6))+.6exp(-((y+.27-.3z)/.07)^2)*exp(-((x/.7)^6+(z/.65)^6)) : 0
  p.pressure[i,j,k]=q; p.previous[i,j,k]=q
  mineral=mode=="cavity" ? .8exp(-((sqrt((x-.13)^2+y^2+z^2)-.49)/.035)^2)*(1-.85exp(-((y-.28)^2+(z-.18)^2)/.035)) : mode=="fronts" ? .25*(.5+.5sin(13z+3y))^5 : mode=="lattice" ? .32*(cos(12x)*cos(12y)*cos(12z))^6 : 0
  w.fields.deposit[i,j,k]=mineral; w.fields.nutrient[i,j,k]=1-mineral
  if neural
   phase=mode=="scroll" ? 3atan(y,x)+12z+8r : 3sin(5x)*cos(5y)+8z+5r
   env=exp(-(r/.83)^8)
   nf.latent1[i,j,k]=.8sin(phase)*env; nf.latent2[i,j,k]=.8cos(phase)*env
  end
 end
 dep=Deposition(rate=.7,threshold=neural ? .28 : .18,dissolution=.0001)
 procs=neural ? (nf,dep) : (p,Excitable(diffusion=.00035,threshold=.2,recovery_rate=.15,gain=2.5,material_feedback=.7),dep)
 initialmass=material_mass(w); elapsed=@elapsed for step=1:steps
  step!(w,procs,dt)
  step%50==0 && println("$slug $step/$steps t=$(w.time)")
 end
 all(a->all(isfinite,a),values(w.fields)) || error("nonfinite")
 base=optical_fields(w); crest=neural ? max.(Array(w.fields.activity),0) : abs.(Array(p.pressure)); crestscale=sort(vec(copy(crest)))[ceil(Int,.995length(crest))]; crest=clamp.(crest./max(crestscale,1f-8),0,1)
 # Explicit artistic transfer retains the unmodified pressure and mineral channels.
 fields=(density=.035f0 .+ .9f0.*Array(w.fields.deposit).+.38f0.*crest.^2,activity=crest.^2,structure=Array(w.fields.deposit),pressure=Array(w.fields.pressure),excitation=Array(w.fields.activity),mask=base.mask)
 export_volume(fields,joinpath(out,"volume")); record_state(w,procs,joinpath(out,"state");step=steps)
 recipe=Dict("id"=>slug,"title"=>title,"mode"=>mode,"resolution"=>n,"steps"=>steps,"dt"=>dt,"study"=>"studies/collection_waves.jl","model"=>neural ? "fixed recurrent tanh; no training" : "damped scalar wave plus excitable conservative deposition; not calibrated acoustics, cavitation or sonoluminescence")
 for path in (joinpath(out,"recipe.toml"),joinpath(ROOT,"experiments","collection","waves",slug*".toml"))
  mkpath(dirname(path)); open(io->TOML.print(io,recipe),path,"w")
 end
 metrics=(pressure_extrema=extrema(w.fields.pressure),activity_extrema=extrema(w.fields.activity),deposit_extrema=extrema(w.fields.deposit),initial_mass=initialmass,final_mass=material_mass(w),crest_fraction=sum(crest.>.2)/length(crest))
 open(io->JSON3.pretty(io,(id=slug,family="waves",shape=w.grid.size,title=title,description=description,time=w.time,params=recipe,metrics=metrics,elapsed_seconds=elapsed,sourcehash=bytes2hex(sha256(read(@__FILE__))),optical_transfer="99.5th-percentile-clamped absolute pressure squared (or positive neural activity squared) drives emissive surrogate; density combines mineral and crest; no physical light emission claim")),joinpath(out,"run.json"),"w")
 args=["--palette",palette,"--emission","0.5","--absorption","0.7","--scattering","0.12","--density-scale","0.8","--polished-carrier","--carrier-shape","box","--ior","1.38","--orbit","0.65","--elevation","0.34"]
 open(io->JSON3.pretty(io,(id=slug,family="waves",title=title,description=description,volume="volume",renderer="volume",args=args,size=1600,samples=96)),joinpath(out,"render-job.json"),"w")
 println("READY $slug $metrics")
end
function refreshoptics()
 for spec in SPECS
  slug=spec[1]; out=joinpath(ROOT,"outputs","collection-20261015",slug); path=joinpath(out,"run.json")
  isfile(path) || continue
  meta=JSON3.read(read(path,String),Dict{String,Any}); dims=Tuple(Int.(meta["shape"])); vdir=joinpath(out,"volume")
  field(name)=reshape(open(io->read!(io,Vector{Float32}(undef,prod(dims))),joinpath(vdir,name*".f32")),dims)
  neural=spec[4] in ("scroll","lattice"); raw=neural ? max.(field("excitation"),0) : abs.(field("pressure"))
  scale=sort(vec(copy(raw)))[ceil(Int,.995length(raw))]; crest=clamp.(raw./max(scale,1f-8),0,1)
  deposit=field("structure")
  export_volume((density=.035f0 .+.9f0.*deposit.+.38f0.*crest.^2,activity=crest.^2,structure=deposit,pressure=field("pressure"),excitation=field("excitation"),mask=field("mask")),vdir)
  meta["optical_transfer"]="99.5th-percentile-clamped absolute pressure squared (or positive neural activity squared); emissive surrogate, not physical light emission"
  meta["optical_scale"]=scale; meta["metrics"]["crest_fraction"]=sum(crest.>.2)/length(crest)
  meta["sourcehash"]=bytes2hex(sha256(read(@__FILE__)))
  open(io->JSON3.pretty(io,meta),path,"w")
  println("REFRESHED $slug scale=$scale")
 end
end
if "--refresh" in ARGS
 refreshoptics(); exit()
end
for (index,spec) in enumerate(SPECS)
 !isempty(ARGS) && !(string(index) in split(ARGS[1],',')) && continue
 runone(spec;n=parse(Int,get(ENV,"WAVES_N","96")))
 GC.gc()
end
