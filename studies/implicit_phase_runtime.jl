# Whole coupled brine timing; native equations remain in their original modules.
using Interior, FFTW, LinearAlgebra, Statistics, TOML
using Interior: JSON3
isdefined(@__MODULE__,:brine_initialize!)||include("brine_memory.jl")
struct ReusableImplicitPhase{P,R}
    phase::P
    runtime::R
end
Interior.process_ports(p::ReusableImplicitPhase)=Interior.process_ports(p.phase)
Interior.check_dt(w,p::ReusableImplicitPhase,dt)=Interior.check_dt(w,p.phase,dt)
Interior.advance!(w,p::ReusableImplicitPhase,dt::Real)=Interior.advance!(w,p.phase,dt,p.runtime)
struct UncachedImplicitPhase{P}
    phase::P
end
Interior.process_ports(p::UncachedImplicitPhase)=Interior.process_ports(p.phase)
Interior.check_dt(w,p::UncachedImplicitPhase,dt)=Interior.check_dt(w,p.phase,dt)
Interior.advance!(w,p::UncachedImplicitPhase,dt::Real)=Interior.advance!(w,p.phase,dt,Interior.ImplicitPhaseRuntime(p.phase))
function phase_runtime_fixture(c,n)
 w=World(Grid(n,n,n;extent=ntuple(_->c["extent"],3));fields=(:nutrient,:deposit,:activity),scratch_fields=(:nutrient,))
 p=ImplicitPhaseSeparation(w;noise=0,kappa=c["kappa"],mobility=c["mobility"],mineral_pinning=0,nutrient_diffusion=0,interface_activity=0,activity_relaxation=0,stabilization=c["stabilization"])
 a=zeros(Float32,n,n,n);b=similar(a);scratch=similar(a);permeability=similar(a)
 brine_initialize!(w,p,a,b,c)
 (;w,p,a,b,scratch,permeability,runtime=Interior.ImplicitPhaseRuntime(p))
end
function coupled_runtime_step!(x,c,dt,cached;reference=false)
 (;w,p,a,b,scratch,permeability)=x
 feedback=Float32(c["mineral_phase_feedback"]);phase_barrier=Float32(c["phase_transport_barrier"]);mineral_barrier=Float32(c["mineral_transport_barrier"])
 @. p.chemical_forcing=feedback*w.fields.deposit*p.phi
 step!(w,reference ? ReferenceImplicitPhase(p,x.runtime) : cached ? ReusableImplicitPhase(p,x.runtime) : p,dt)
 @. permeability=1f0/(1f0+phase_barrier*(1f0+tanh(2f0*p.phi))/2f0+mineral_barrier*w.fields.deposit)
 Interior.medium_diffuse!(scratch,a,permeability,w.grid,Float32(c["diffusion_a"]),dt)
 Interior.medium_diffuse!(scratch,b,permeability,w.grid,Float32(c["diffusion_b"]),dt)
 Interior.launch!(w,brine_react!,a,b,w.fields.deposit,p.interface,dt,Float32(c["reaction_rate"]),Float32(c["basal_reaction"]),Float32(c["nucleation_gradient"]),Float32(c["saturation_gradient"]))
 x
end
function benchmark_phase_runtime(output;resolutions=(96,192),steps=5)
 FFTW.set_num_threads(2);BLAS.set_num_threads(1)
 c=TOML.parsefile(joinpath(@__DIR__,"..","experiments/art/brine-lace-96.toml"));results=[]
 for n in resolutions
  dt=Float32(c["dt"]*(96/n)^2);modes=[]
  for cached in (false,true)
   GC.gc();setup=@elapsed x=phase_runtime_fixture(c,n)
   initialphase=sum(Float64,x.p.phi);initialmass=sum(Float64,x.a)+sum(Float64,x.b)
   cold=@elapsed coupled_runtime_step!(x,c,dt,cached)
   times=[@elapsed coupled_runtime_step!(x,c,dt,cached) for _ in 1:steps]
   phase_error=(sum(Float64,x.p.phi)-initialphase)/n^3
   mass_error=(sum(Float64,x.a)+sum(Float64,x.b)+sum(Float64,x.w.fields.deposit)-initialmass)/initialmass
   push!(modes,(;cached,setup_seconds=setup,cold_seconds=cold,step_seconds=times,median_seconds=median(times),phase_error,mass_error,owned_bytes=Base.summarysize(x)))
   x=nothing;GC.gc()
  end
  push!(results,(;resolution=n,dt,steps,modes,speedup=modes[1].median_seconds/modes[2].median_seconds))
  println(JSON3.write(last(results)));flush(stdout)
 end
 open(io->JSON3.pretty(io,(;threads=Threads.nthreads(),fftw_threads=FFTW.get_num_threads(),results)),output,"w")
 results
end


# A source-derived oracle removes only the zero-D guard. All loop equations and
# the complete phase schedule stay sourced from their production definitions.
function install_phase_coupling_reference!()
 source=Meta.parseall(read(joinpath(@__DIR__,"..","src","phase.jl"),String))
 kernel=only(filter(ex->ex isa Expr&&ex.head==:macrocall&&ex.args[1]==Symbol("@kernel")&&ex.args[end].args[1].args[1]==:phase_couple_kernel!,source.args))
 function unguard(ex)
  ex isa Expr||return ex
  ex.head==:if&&ex.args[1]==:(D!=0)&&return unguard(ex.args[2])
  Expr(ex.head,map(unguard,ex.args)...)
 end
 kernel=unguard(deepcopy(kernel));kernel.args[end].args[1].args[1]=:phase_couple_reference_kernel!
 Core.eval(Interior,kernel)
 source=Meta.parseall(read(joinpath(@__DIR__,"..","src","implicit_phase.jl"),String))
 fn=only(filter(ex->ex isa Expr&&ex.head==:function&&ex.args[1] isa Expr&&ex.args[1].head==:call&&ex.args[1].args[1]==:implicit_phase_advance!&&length(ex.args[1].args)==5,source.args))
 replace_symbol(ex)=ex isa Expr ? Expr(ex.head,map(replace_symbol,ex.args)...) : ex===:phase_couple_kernel! ? :phase_couple_reference_kernel! : ex
 fn=replace_symbol(deepcopy(fn));fn.args[1].args[1]=:phase_runtime_reference_advance!
 insert!(fn.args[2].args,1,:(check_dt(w,p,dt)))
 Core.eval(Interior,fn)
 nothing
end
struct ReferenceImplicitPhase{P,R};phase::P;runtime::R;end
Interior.process_ports(p::ReferenceImplicitPhase)=Interior.process_ports(p.phase)
Interior.check_dt(w,p::ReferenceImplicitPhase,dt)=Interior.check_dt(w,p.phase,dt)
Interior.advance!(w,p::ReferenceImplicitPhase,dt::Real)=Base.invokelatest(getfield(Interior,:phase_runtime_reference_advance!),w,p.phase,dt,p.runtime)
function benchmark_zero_diffusion_phase(output;resolutions=(96,192),steps=10)
 install_phase_coupling_reference!()
 c=TOML.parsefile(joinpath(@__DIR__,"..","experiments/art/brine-lace-96.toml"));rows=[]
 for n in resolutions
  x=phase_runtime_fixture(c,n);dt=Float32(.04*(96/n)^2)
  coupled_runtime_step!(x,c,dt,true;reference=true);coupled_runtime_step!(x,c,dt,true)
  times=[]
  for i in 1:steps
   reference=isodd(i);seconds=@elapsed coupled_runtime_step!(x,c,dt,true;reference)
   push!(times,(;reference,seconds))
  end
  push!(rows,(;n,dt,times));println(JSON3.write(last(rows)));flush(stdout)
  x=nothing;GC.gc()
 end
 open(io->JSON3.pretty(io,rows),output,"w")
 rows
end

if abspath(PROGRAM_FILE)==@__FILE__;benchmark_phase_runtime(ARGS[1]);end
