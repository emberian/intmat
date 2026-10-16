# Complete CPU/native phase and coupled brine gate, with identical observation
# schedules. No alternative physical evolution or study reaction is implemented.
started=time_ns()
include(joinpath(@__DIR__,"device_cellular_setup.jl"))
include(joinpath(@__DIR__,"..","studies","implicit_phase_runtime.jl"))
using SHA, Serialization
module CPUImplicitPhaseTests
using Interior, Test, Random
if get(ENV,"PHASE_SKIP_CPU_GATE","0")!="1"
include(joinpath(@__DIR__,"..","test","implicit_phase.jl"))
include(joinpath(@__DIR__,"..","test","implicit_phase_runtime.jl"))
include(joinpath(@__DIR__,"..","test","compartments.jl"))
end
end
const RECIPE=TOML.parsefile(joinpath(@__DIR__,"..","experiments/art/brine-lace-96.toml"))
function device_fixture(host)
 w=World(host.w.grid;backend,seed=host.w.seed,fields=keys(host.w.fields),scratch_fields=keys(host.w.scratch))
 for k in keys(w.fields);copyto!(w.fields[k],host.w.fields[k]);end
 for k in keys(w.scratch);copyto!(w.scratch[k],host.w.scratch[k]);end
 p=ImplicitPhaseSeparation(w;noise=0,kappa=host.p.kappa,mobility=host.p.mobility,mineral_pinning=0,
  nutrient_diffusion=host.p.nutrient_diffusion,interface_activity=host.p.interface_activity,activity_relaxation=host.p.activity_relaxation,stabilization=host.p.stabilization)
 copyto!(p.phi,host.p.phi);copyto!(p.chemical_forcing,host.p.chemical_forcing);Interior.phase_observe!(w,p)
 upload(a)=begin b=Interior.KernelAbstractions.zeros(backend,eltype(a),size(a)...);copyto!(b,a);b end
 (;w,p,a=upload(host.a),b=upload(host.b),scratch=upload(host.scratch),permeability=upload(host.permeability),runtime=Interior.ImplicitPhaseRuntime(p))
end
function host_observation(x)
 (;phi=Array(x.p.phi),interface=Array(x.p.interface),chemical_potential=Array(x.p.chemical_potential),
  a=Array(x.a),b=Array(x.b),mineral=Array(x.w.fields.deposit),nutrient=Array(x.w.fields.nutrient),activity=Array(x.w.fields.activity),time=x.w.time)
end
function check_observation(a,b;spacing,kappa,bulk,tolerance=2e-5)
 a.time==b.time||error("world clocks differ")
 names=keys(Base.structdiff(a,(time=a.time,)))
 errors=NamedTuple{names}(Tuple(maximum(abs,Float64.(a[k]).-Float64.(b[k])) for k in names))
 all(k->k==:chemical_potential||errors[k]<tolerance,names)||error("native primary physical parity failed: $errors")
 # Chemical potential includes a discrete Laplacian: high-resolution derivatives
 # amplify one composition ULP. Bound that derived error by its actual operator.
 maxphi=max(maximum(abs,a.phi),maximum(abs,b.phi))
 maxmu=max(1.,maximum(abs,a.chemical_potential),maximum(abs,b.chemical_potential))
 mu_bound=(Float64(bulk)*(3Float64(maxphi)^2+1)+4Float64(kappa)*sum(inv(Float64(h)^2) for h in spacing))*errors.phi+16eps(Float32)*maxmu
 errors.chemical_potential<=mu_bound||error("chemical-potential error exceeds composition/rounding bound: $errors > $mu_bound")
 all(isfinite,b.phi)&&minimum(b.a)>=0&&minimum(b.b)>=0||error("native finite/positivity failed")
 (;max_error=maximum(errors),errors,chemical_potential_bound=mu_bound)
end
function checkpoint_restart(x)
 # One alias-preserving host physical checkpoint, with no native resources.
 cp=Interior.implicit_phase_checkpoint(x.w,x.p)
 io=IOBuffer();serialize(io,(cp,cp.phase,Array(x.a),Array(x.b),Array(x.scratch),Array(x.permeability)));seekstart(io)
 cp,alias,ha,hb,hs,hperm=deserialize(io)
 cp.phase===alias||error("host checkpoint process alias was lost")
 w,p=Interior.implicit_phase_restore(cp;backend)
 p.real_buffer!==x.p.real_buffer||error("restart did not reconstruct arrays")
 upload(a)=begin b=Interior.KernelAbstractions.zeros(backend,eltype(a),size(a)...);copyto!(b,a);b end
 (;w,p,a=upload(ha),b=upload(hb),scratch=upload(hs),permeability=upload(hperm),runtime=Interior.ImplicitPhaseRuntime(p))
end
function phase_boundaries()
 checks=0
 for grid in (Grid(9,7,5),Grid(1,9,7),Grid(1,1,1))
  cw=World(grid);gw=World(grid;backend)
  fill!(cw.fields.deposit,2);fill!(gw.fields.deposit,2)
  kw=(noise=0,mean=.3,mineral_pinning=.7,nutrient_diffusion=.01,activity_relaxation=.2)
  cp=ImplicitPhaseSeparation(cw;kw...);gp=ImplicitPhaseSeparation(gw;kw...)
  fill!(cp.chemical_forcing,4);fill!(gp.chemical_forcing,4)
  cr=Interior.ImplicitPhaseRuntime(cp);gr=Interior.ImplicitPhaseRuntime(gp)
  initial=copy(cp.phi)
  for _ in 1:3
   step!(cw,ReusableImplicitPhase(cp,cr),.001)
   step!(gw,ReusableImplicitPhase(gp,gr),.001)
  end
  maximum(abs,cp.phi.-Array(gp.phi))<2e-7||error("odd/singleton native phase differs")
  maximum(abs,cp.phi.-initial)<2e-7||error("constant chemical forcing changed zero mode")
  cw.time==gw.time||error("native boundary clock differs")
  caught=false
  try Interior.advance!(gw,gp,.001,cr);catch e;caught=e isa ArgumentError;end
  caught||error("wrong backend runtime accepted")
  caught=false
  try Interior.advance!(gw,gp,-1);catch e;caught=e isa ArgumentError;end
  caught||error("invalid dt accepted")
  # Preserve a real physical alias through device-host-device transfer.
  aliased=ImplicitPhaseSeparation(gp.reference,gp.phi,gp.stabilization,gp.eigenvalues,gp.real_buffer,gp.spectrum,gp.nonlinear_spectrum)
  checkpoint=Interior.implicit_phase_checkpoint(gw,aliased)
  aw,ap=Interior.implicit_phase_restore(checkpoint;backend)
  ap.chemical_forcing===ap.phi||error("device checkpoint field alias lost")
  ap.chemical_forcing .= .25f0
  maximum(abs,Array(ap.phi).-.25f0)==0||error("restored alias did not update physical buffer")
  checks+=7
 end
 checks
end
function memory_snapshot()
 free,total=AMDGPU.info()
 scope=split(only(readlines("/proc/self/cgroup")),":";limit=3)[3]
 peak=try parse(Int,strip(read("/sys/fs/cgroup"*scope*"/memory.peak",String))) catch;nothing end
 (;device_free=free,device_total=total,device_used=total-free,host_cgroup_peak=peak)
end
function runtime_memory(r)
 arrays=(r.forward.buffer,r.forward.workarea,r.inverse_real.buffer,r.inverse_real.workarea)
 (;forward_buffer=length(arrays[1])*sizeof(eltype(arrays[1])),forward_workspace=length(arrays[2])*sizeof(eltype(arrays[2])),
  inverse_buffer=length(arrays[3])*sizeof(eltype(arrays[3])),inverse_workspace=length(arrays[4])*sizeof(eltype(arrays[4])))
end
function transport_boundaries()
 checks=0
 upload(a)=begin b=Interior.KernelAbstractions.zeros(backend,eltype(a),size(a)...);copyto!(b,a);b end
 for dims in ((1,1,1),(1,7,3),(9,1,4),(8,9,1),(8,9,5),(2,2,2)),mode in (:unequal,:closed,:constant,:zero_diffusion)
  grid=Grid(dims...;extent=(3,4,5))
  a=reshape(Float32.(abs.(sin.(1:prod(dims)))),dims)
  perm=reshape(Float32.(abs.(cos.(1:prod(dims)))),dims)
  mode==:closed&&fill!(perm,0);mode==:constant&&fill!(a,.3f0)
  d=mode==:zero_diffusion ? 0f0 : .03f0
  cpu=copy(a);gpu=upload(a);gp=upload(perm);scratch=similar(gpu);out=similar(cpu)
  for _ in 1:5
   Interior.medium_diffuse!(out,cpu,perm,grid,d,.001f0)
   Interior.medium_diffuse!(scratch,gpu,gp,grid,d,.001f0)
  end
  observed=Array(gpu)
  maximum(abs,cpu.-observed)<2e-7||error("native periodic/anisotropic transport differs")
  abs(sum(Float64,observed)-sum(Float64,a))<2e-6||error("native transport mass drift")
  minimum(observed)>=0||error("native transport lost positivity")
  caught=false
  try Interior.medium_diffuse!(gpu,gpu,gp,grid,d,.001f0);catch e;caught=e isa ArgumentError;end
  caught||error("device diffusion accepted aliased scratch")
  caught=false
  try Interior.medium_diffuse!(gp,gpu,gp,grid,d,.001f0);catch e;caught=e isa ArgumentError;end
  caught||error("device diffusion accepted permeability output alias")
  checks+=5
 end
 checks
end
function main(output)
baseline_memory=memory_snapshot();boundary_checks=phase_boundaries();transport_checks=transport_boundaries();reports=[]
for n in (96,192)
 dt=Float32(.04*(96/n)^2)
 setup_cpu=@elapsed cpu=phase_runtime_fixture(RECIPE,n)
 setup_gpu=@elapsed gpu=device_fixture(cpu)
 memory_after_setup=memory_snapshot();runtime_bytes=runtime_memory(gpu.runtime)
 initial=host_observation(cpu);initialphase=sum(Float64,initial.phi);initialmass=sum(Float64,initial.a)+sum(Float64,initial.b)
 trials=[];restarts=[];maxerror=0.
 for mode in (:phase,:brine)
  for trial in 1:4
   count=trial==1 ? 1 : 5
   ct=@elapsed begin
    for _ in 1:count
     if mode==:phase;step!(cpu.w,ReusableImplicitPhase(cpu.p,cpu.runtime),dt)
     else coupled_runtime_step!(cpu,RECIPE,dt,true);end
    end
    co=host_observation(cpu)
   end
   gt=@elapsed begin
    for _ in 1:count
     if mode==:phase;step!(gpu.w,ReusableImplicitPhase(gpu.p,gpu.runtime),dt)
     else coupled_runtime_step!(gpu,RECIPE,dt,true);end
    end
    go=host_observation(gpu)
   end
   parity=check_observation(co,go;spacing=cpu.w.grid.spacing,kappa=cpu.p.kappa,bulk=cpu.p.bulk)
   maxerror=max(maxerror,parity.max_error)
   push!(trials,(;mode,trial,count,cpu_seconds=ct,gpu_seconds=gt,speedup=ct/gt,errors=parity.errors,chemical_potential_bound=parity.chemical_potential_bound))
   println(JSON3.write((;n,last=last(trials),maxerror)));flush(stdout)
  end
  # A restored state must rebuild resources and continue from the same phase.
  before=host_observation(gpu);rt=@elapsed gpu=checkpoint_restart(gpu)
  push!(restarts,(;mode,seconds=rt))
  before==host_observation(gpu)||error("restart changed physical payload")
 end
 final=host_observation(gpu)
 phase_error=(sum(Float64,final.phi)-initialphase)/n^3
 mass_error=(sum(Float64,final.a)+sum(Float64,final.b)+sum(Float64,final.mineral)-initialmass)/initialmass
 abs(phase_error)<2e-6&&abs(mass_error)<2e-6||error("native conservation failed")
 push!(reports,(;n,dt,setup_cpu,setup_gpu,memory_after_setup,memory_after_run=memory_snapshot(),runtime_bytes,trials,restarts,maxerror,phase_error,mass_error,time=gpu.w.time,
  cpu_owned_bytes=Base.summarysize(cpu),physical_bytes=sum(length(a)*sizeof(eltype(a)) for a in (gpu.p.phi,gpu.p.chemical_potential,gpu.p.interface,gpu.p.scratch,gpu.p.chemical_forcing,gpu.p.eigenvalues,gpu.p.real_buffer,gpu.p.spectrum,gpu.p.nonlinear_spectrum,gpu.a,gpu.b,gpu.scratch,gpu.permeability,values(gpu.w.fields)...,values(gpu.w.scratch)...))))
 cpu=nothing;gpu=nothing;GC.gc();Interior.KernelAbstractions.synchronize(backend)
end
open(io->JSON3.pretty(io,(;elapsed_seconds=(time_ns()-started)/1e9,baseline_memory,boundary_checks,transport_checks,reports,precision="physical Float32; FFT Float64/ComplexF64",source=Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,"..",f)))) for f in ("src/implicit_phase.jl","src/phase.jl","src/compartments.jl","studies/brine_memory.jl")))),output,"w")
println("native complete phase and brine gate passed");flush(stdout)

end
main(ARGS[1])
