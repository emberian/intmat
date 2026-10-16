# Finite precipitation writes a chemical-potential memory into a conserved phase.
using Interior, TOML, SHA, Statistics, FFTW, LinearAlgebra
using Interior: JSON3
using Interior.KernelAbstractions: @kernel, @index
BLAS.set_num_threads(1);FFTW.set_num_threads(min(2,Threads.nthreads()))
@kernel function brine_react!(a,b,mineral,interface,dt,rate,basal,nucleation_gradient,saturation_gradient)
 i=@index(Global,Linear)
 gate=basal+(1f0-basal)*clamp((interface[i]-nucleation_gradient)/(saturation_gradient-nucleation_gradient),0f0,1f0)
 extent=min(a[i],b[i],dt*rate*gate*a[i]*b[i])
 a[i]-=extent;b[i]-=extent;mineral[i]+=2f0*extent
end
function brine_initialize!(w,p,a,b,c)
 L=Float32(c["extent"]);h=L/size(a,1)
 @inbounds for I in CartesianIndices(a)
  x,y,z=ntuple(d->(I[d]-.5f0)*h-L/2,3);u,v,s=2f0*Float32(pi).*(x,y,z)./L
  if get(c,"phase_initialization","lamellar")=="spinodal"
   modes=((3,2,1),(2,3,-1),(1,-2,3),(-3,1,2),(2,-2,2),(3,1,-2),(1,3,2),(-2,2,3),(4,1,0),(0,3,3),(3,0,-3),(2,2,-3))
   perturbation=0f0
   for (q,mode) in enumerate(modes)
    perturbation+=sin(mode[1]*u+mode[2]*v+mode[3]*s+Float32(.37q+.113q*q))
   end
   p.phi[I]=Float32(c["phase_perturbation_rms"])*perturbation/sqrt(6f0)
  else
   nodal=.8f0*sin(2s+.45f0*sin(2u)+.5f0*cos(v))+.5f0*cos(2u)*cos(2v)+.15f0*sin(3u)
   p.phi[I]=tanh(2.2f0*nodal)
  end
  exponent=get(c,"reservoir_shape","sheet")=="lens" ? 6 : 8
  sheet=exp(-((z-.9f0*sin(.35f0*x)-.6f0*cos(.4f0*y))/Float32(c["reservoir_thickness"]))^4-(x/8f0)^exponent-(y/7f0)^exponent)
  a[I]=Float32(c["a_inventory"])*sheet;b[I]=Float32(c["b_inventory"])*(1f0-sheet)
 end
 Interior.phase_observe!(w,p)
end
function brine_arm(c,out,arm;n=c["resolution"],steps=c["steps"])
 ispath(out)&&error("Fresh output required: $out");mkpath(out)
 0<=get(c,"nucleation_gradient",0)<get(c,"saturation_gradient",1/c["interface_reaction_gain"]) || error("Ordered physical gradient thresholds required")
 estimate=128n^3;estimate<=c["memory_budget_mb"]*1024^2||error("Memory plan exceeds recipe")
 sources=("studies/brine_memory.jl","src/phase.jl","src/implicit_phase.jl","src/compartments.jl","src/fields.jl")
 for f in sources;cp(f,joinpath(out,basename(f)));end
 w=World(Grid(n,n,n;extent=ntuple(_->c["extent"],3));fields=(:nutrient,:deposit,:activity),scratch_fields=(:nutrient,))
 p=ImplicitPhaseSeparation(w;noise=0,kappa=c["kappa"],mobility=c["mobility"],mineral_pinning=0,nutrient_diffusion=0,interface_activity=0,activity_relaxation=0,stabilization=c["stabilization"])
 a=zeros(Float32,n,n,n);b=similar(a);scratch=similar(a);permeability=similar(a)
 brine_initialize!(w,p,a,b,c);dt=Float32(c["dt"]);Da=Float32(c["diffusion_a"]);Db=Float32(c["diffusion_b"]);h=w.grid.spacing[1]
 dt*6max(Da,Db)/h^2<.9||error("Diffusion positivity bound")
 feedback=arm in ("closed","phase-only") ? Float32(c["mineral_phase_feedback"]) : 0f0
 barrier=arm in ("closed","transport-only") ? Float32(c["mineral_transport_barrier"]) : 0f0
 volume=Float64(h)^3;initialmass=(sum(Float64,a)+sum(Float64,b))*volume;initialphase=sum(Float64,p.phi)*volume
 history=Any[];frames=Any[];times=Float64[]
 function capture(step)
  directory=joinpath(out,"frames",lpad(step,5,'0'));mkpath(directory)
  for (name,field) in (("a",a),("b",b),("mineral",w.fields.deposit),("phi",p.phi),("interface",p.interface))
   write(joinpath(directory,name*".f32"),field)
  end
  total=(sum(Float64,a)+sum(Float64,b)+sum(Float64,w.fields.deposit))*volume
  row=(step=step,time=w.time,mineral_mass=sum(Float64,w.fields.deposit)*volume,mineral_max=maximum(w.fields.deposit),phase_mass=sum(Float64,p.phi)*volume,phase_mass_error=(sum(Float64,p.phi)*volume-initialphase)/c["extent"]^3,initial_chemical_mass=initialmass,chemical_mass=total,chemical_mass_absolute_error=total-initialmass,chemical_mass_relative_error=(total-initialmass)/initialmass,phase_range=extrema(p.phi),interface_integral=sum(Float64,p.interface)*volume,min_a=minimum(a),min_b=minimum(b))
  all(isfinite,p.phi)&&row.min_a>=0&&row.min_b>=0&&abs(row.chemical_mass_relative_error)<2e-6&&abs(row.phase_mass_error)<2e-6||error("Conservation/finite-state failure")
  push!(history,row);push!(frames,(step=step,time=w.time,directory=relpath(directory,out)));println(JSON3.write((;arm,row...)));flush(stdout)
 end
 capture(0)
 for step in 1:steps
  elapsed=@elapsed begin
   @. p.chemical_forcing=feedback*w.fields.deposit*p.phi
   step!(w,p,dt)
   @. permeability=1f0/(1f0+Float32(c["phase_transport_barrier"])*(1f0+tanh(2f0*p.phi))/2f0+barrier*w.fields.deposit)
   Interior.medium_diffuse!(scratch,a,permeability,w.grid,Da,dt)
   Interior.medium_diffuse!(scratch,b,permeability,w.grid,Db,dt)
   Interior.launch!(w,brine_react!,a,b,w.fields.deposit,p.interface,dt,Float32(c["reaction_rate"]),Float32(c["basal_reaction"]),Float32(get(c,"nucleation_gradient",0)),Float32(get(c,"saturation_gradient",1/c["interface_reaction_gain"])))
  end
  step>1&&push!(times,elapsed)
  (step%max(1,steps÷3)==0||step==steps)&&capture(step)
 end
 report=(execution_base_commit=get(ENV,"BRINE_EXECUTION_COMMIT","unrecorded"),id=c["id"],title=c["title"],family="interfaces",arm=arm,shape=(n,n,n),extent=ntuple(_->c["extent"],3),spacing=Tuple(w.grid.spacing),dt=dt,steps=steps,time=w.time,frames=frames,history=history,median_step_seconds=median(times),evolution_seconds=sum(times),estimated_working_bytes=estimate,recipe=c,source_sha256=Dict(f=>bytes2hex(sha256(read(joinpath(out,basename(f))))) for f in sources),simulation="Existing stabilized semi-implicit periodic conserved Cahn–Hilliard phase. Existing conservative harmonic-face solute diffusion. Study-local bounded A+B→2 mineral reaction, gated by actual phase gradient and recipe nucleation threshold. New mineral lowers later permeability and supplies chemical-potential derivative feedback*mineral*phi, pinning mixing interfaces. Finite chemical inventory and signed phase composition are separately conserved; no thermodynamic energy closure or calibrated chemistry is claimed.",initialization="Recipe-selected physical Fourier perturbation or corrugated lamellar nodal field; A in one finite bent lens/sheet reservoir and B outside it. No pre-existing mineral. Optical boundaries are observations, not added simulation walls.")
 open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w");report
end
function brine_study(recipe,out;pilot=false)
 c=TOML.parsefile(recipe);ispath(out)&&error("Fresh study output required");mkpath(out);cp(recipe,joinpath(out,"recipe.toml"))
 closed=brine_arm(c,joinpath(out,"closed"),"closed")
 if pilot
  control=brine_arm(c,joinpath(out,"no-mineral-feedback"),"no-mineral-feedback")
  readfield(arm,name)=reinterpret(Float32,read(joinpath(out,arm,last(closed.frames).directory,name*".f32")))
  phierr=maximum(abs.(readfield("closed","phi").-readfield("no-mineral-feedback","phi")));merr=maximum(abs.(readfield("closed","mineral").-readfield("no-mineral-feedback","mineral")))
  comparison=(schema="brine-mineral-feedback-ablation",shape=closed.shape,time=closed.time,phase_max_difference=phierr,mineral_max_difference=merr,closed_mineral_mass=last(closed.history).mineral_mass,no_feedback_mineral_mass=last(control.history).mineral_mass,interpretation="Same initialization and reaction; remove only mineral-to-phase forcing and mineral-to-transport barrier. Demonstrates a material history affects later dynamics, not emergent biological behavior.")
  open(io->JSON3.pretty(io,comparison),joinpath(out,"comparison.json"),"w");println(JSON3.write(comparison))
 end
end
if abspath(PROGRAM_FILE)==@__FILE__;brine_study(ARGS[1],ARGS[2];pilot="--pilot" in ARGS);end
