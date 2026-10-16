# New resolution producer: existing brine equations, worker-owned FFT resources,
# bounded host/device memory and streamed artifacts. No historical run is loaded.
const NATIVE_BRINE_STARTED=time_ns()
include(joinpath(@__DIR__,"..","scripts","device_cellular_setup.jl"))
include("implicit_phase_runtime.jl")
using SHA
const NATIVE_BRINE_SOURCES=("src/core.jl","src/phase.jl","src/implicit_phase.jl","src/compartments.jl",
 "studies/brine_memory.jl","studies/implicit_phase_runtime.jl","studies/native_brine_producer.jl","scripts/device_cellular_setup.jl")
@kernel function brine_sample_kernel!(out,field,n,m)
 i,j,k=@index(Global,NTuple)
 q=ntuple(d->((2*(i,j,k)[d]-1)*n)÷(2m)+1,3)
 out[i,j,k]=field[q...]
end
function producer_memory()
 free,total=AMDGPU.info()
 scope=split(only(readlines("/proc/self/cgroup")),":";limit=3)[3]
 readvalue(name)=parse(Int,strip(read("/sys/fs/cgroup"*scope*"/"*name,String)))
 (;device_free=free,device_total=total,host_current=readvalue("memory.current"),host_peak=readvalue("memory.peak"))
end
function producer_upload(host)
 w,p=Interior.implicit_phase_transfer(host.w,host.p,backend)
 upload(a)=begin out=Interior.KernelAbstractions.zeros(backend,eltype(a),size(a)...);copyto!(out,a);out end
 (;w,p,a=upload(host.a),b=upload(host.b),scratch=upload(host.scratch),permeability=upload(host.permeability),runtime=Interior.ImplicitPhaseRuntime(p))
end
function producer_scalars(x,initial;strict=true)
 s64(a)=sum(v->Float64(v),a)
 phase=s64(x.p.phi);chemical=s64(x.a)+s64(x.b)+s64(x.w.fields.deposit)
 phase_error=(phase-initial.phase)/length(x.p.phi)
 chemical_error=(chemical-initial.chemical)/initial.chemical
 arange=extrema(x.a);brange=extrema(x.b);prange=extrema(x.p.phi);mrange=extrema(x.w.fields.deposit)
 valid=all(isfinite,prange)&&all(isfinite,arange)&&all(isfinite,brange)&&all(isfinite,mrange)&&mrange[1]>=0&&arange[1]>=0&&brange[1]>=0&&abs(phase_error)<2e-6&&abs(chemical_error)<2e-6
 strict&&!valid&&error("brine finite/positivity/conservation gate failed: phase=$phase_error chemical=$chemical_error")
 (;valid,time=x.w.time,mineral_sum=s64(x.w.fields.deposit),mineral_max=maximum(x.w.fields.deposit),phase_mean=phase/length(x.p.phi),
  phase_mean_error=phase_error,chemical_relative_error=chemical_error,phase_range=prange,a_range=arange,b_range=brange,mineral_range=mrange)
end
producer_chemical_type(c)=get(c,"chemical_storage","Float32")=="Float64" ? Float64 : get(c,"chemical_storage","Float32")=="Float32" ? Float32 : error("unsupported chemical_storage")
function producer_preflight(c,out)
 n=Int(c["resolution"]);N=n^3;H=(n÷2+1)*n*n;memory=producer_memory()
 chemical_type=producer_chemical_type(c);extra=4N*(sizeof(chemical_type)-4)
 physical=60N+36H+extra;native_active=60N+68H+extra
 # Upper-bound initialization and transfer transients beyond the owned CPU
 # fixture: two Float32 fields, 256MiB staging, plus 512MiB safety margin.
 host_estimate=memory.host_current+physical+8N+768*1024^2
 device_margin=2*1024^3
 m=Int(c["replay_resolution"]);frames=1+cld(Int(c["steps"]),Int(c["replay_stride"]))
 artifact_estimate=N*(4+3sizeof(chemical_type))+frames*4m^3*4+1024^2
 h=Float64(c["extent"])/n;dt=Float32(c["dt"])
 positivity=Float64(dt)*6max(c["diffusion_a"],c["diffusion_b"])/h^2
 positivity<.9||error("explicit solute positivity bound exceeded")
 host_estimate<=8*1024^3||error("host initialization preflight exceeds8GiB")
 native_active+device_margin<memory.device_free||error("native active state plus2GiB margin exceeds available device memory")
 artifact_estimate<=Int(c["artifact_budget_bytes"])||error("artifact preflight exceeds budget")
 (;memory,physical_bytes=physical,native_active_bytes=native_active,host_peak_estimate=host_estimate,device_margin,artifact_estimate,
  dt,positivity,planned_steps=Int(c["steps"]),planned_physical_time=Float64(dt)*Int(c["steps"]))
end
function producer_capture(x,out,step,m,frames;full=false)
 directory=joinpath(out,full ? "final" : "replay",lpad(step,6,'0'));mkpath(directory)
 fields=(phi=x.p.phi,a=x.a,b=x.b,mineral=x.w.fields.deposit)
 sampler=full ? nothing : Interior.KernelAbstractions.zeros(backend,Float32,m,m,m)
 for (name,field) in pairs(fields)
  if full
   # One physical field at a time; never a full host observation/checkpoint.
   let host=Array(field);write(joinpath(directory,string(name)*(eltype(field)===Float64 ? ".f64" : ".f32")),host);end
  else
   brine_sample_kernel!(backend,256)(sampler,field,size(field,1),m;ndrange=(m,m,m))
   Interior.KernelAbstractions.synchronize(backend)
   write(joinpath(directory,string(name)*".f32"),Array(sampler))
  end
  GC.gc(false)
 end
 push!(frames,(;step,time=x.w.time,kind=full ? "final" : "replay",shape=full ? size(x.p.phi) : (m,m,m),dtypes=full ? (phi="Float32",a=string(eltype(x.a)),b=string(eltype(x.b)),mineral=string(eltype(x.w.fields.deposit))) : (phi="Float32",a="Float32",b="Float32",mineral="Float32"),directory=relpath(directory,out)))
end
function native_brine_producer(recipe,out;preflight_only=false,capture=producer_capture)
 ispath(out)&&error("fresh output directory required");mkpath(out)
 c=TOML.parsefile(recipe);cp(recipe,joinpath(out,"recipe.toml"))
 sources=Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,"..",f)))) for f in NATIVE_BRINE_SOURCES)
 for (role,file) in ((:implicit,"src/implicit_phase.jl"),(:phase,"src/phase.jl"),(:core,"src/core.jl"))
  Interior.IMPLICIT_PHASE_CHECKPOINT_SOURCE[role]==sources[file]||error("loaded phase source differs from frozen producer")
 end
 # Compile the actual shared execution path, scalar reductions and replay
 # sampler at small resolution before considering the large allocation.
 small=phase_runtime_fixture(c,16;chemical_type=producer_chemical_type(c));warm=producer_upload(small)
 initial=(phase=sum(Float64,small.p.phi),chemical=sum(Float64,small.a)+sum(Float64,small.b))
 warm_seconds=@elapsed begin
  coupled_runtime_step!(warm,c,Float32(c["dt"]),true)
  producer_scalars(warm,initial)
  capture(warm,joinpath(out,"warmup"),0,8,Any[])
 end
 small=nothing;warm=nothing;GC.gc();Interior.KernelAbstractions.synchronize(backend)
 preflight=producer_preflight(c,out)
 open(io->JSON3.pretty(io,(;preflight,warm_seconds,sources,recipe_sha256=bytes2hex(sha256(read(recipe))))),joinpath(out,"preflight.json"),"w")
 println(JSON3.write((;event="preflight_passed",preflight,warm_seconds)));flush(stdout)
 preflight_only&&return preflight
 n=Int(c["resolution"]);dt=Float32(c["dt"]);steps=Int(c["steps"])
 initial_setup=@elapsed host=phase_runtime_fixture(c,n;chemical_type=producer_chemical_type(c))
 initial=(phase=sum(Float64,host.p.phi),chemical=sum(Float64,host.a)+sum(Float64,host.b))
 transfer_setup=@elapsed x=producer_upload(host)
 host=nothing;GC.gc();Interior.KernelAbstractions.synchronize(backend)
 after_upload=producer_memory()
 after_upload.device_free>1024^3||error("less than1GiB device headroom after initialization")
 after_upload.host_current<7*1024^3||error("less than1GiB host headroom after initialization")
 frames=Any[];history=Any[];times=Float64[];completed=0;stopped="planned_steps"
 m=Int(c["replay_resolution"]);history_stride=Int(c["history_stride"]);replay_stride=Int(c["replay_stride"])
 evolution_started=time_ns();failure=nothing;export_failure=nothing
 try
 push!(history,(step=0,producer_scalars(x,initial;strict=false)...))
 capture(x,out,0,m,frames)
 for step in 1:steps
  # Preserve a useful final artifact before the external20-minute deadline.
  if (time_ns()-NATIVE_BRINE_STARTED)/1e9>Float64(c["wall_seconds"])
   stopped="wall_budget";break
  end
  elapsed=@elapsed coupled_runtime_step!(x,c,dt,true)
  step>5&&push!(times,elapsed);completed=step
  if step%history_stride==0||step==steps
   row=(step=step,producer_scalars(x,initial;strict=false)...);push!(history,row)
   println(JSON3.write((;event="history",row...,elapsed_seconds=(time_ns()-NATIVE_BRINE_STARTED)/1e9)));flush(stdout)
   if !row.valid;stopped="numerical_guard";break;end
  end
  (step%replay_stride==0||step==steps)&&capture(x,out,step,m,frames)
 end
 catch err
  stopped="execution_error";failure=sprint(showerror,err,catch_backtrace())
 end
 try
 isempty(history)||last(history).step==completed||push!(history,(step=completed,producer_scalars(x,initial;strict=false)...))
 catch err
  failure=something(failure, "final scalar observation failed: "*sprint(showerror,err,catch_backtrace()))
 end
 try
 capture(x,out,completed,m,frames;full=true)
 catch err
  export_failure=sprint(showerror,err,catch_backtrace())
 end
 # Write failure metadata even when diagnostic export itself fails. These are
 # observations of stopped state, never resumable checkpoints.
 success=stopped=="planned_steps"&&failure===nothing&&export_failure===nothing&&!isempty(history)&&last(history).valid
 open(io->JSON3.pretty(io,(;success,stopped,completed_steps=completed,physical_time=x.w.time,failure,export_failure,resumable=false)),joinpath(out,"status.json"),"w")
 final_memory=producer_memory()
 bytes=sum(filesize(joinpath(root,f)) for (root,_,files) in walkdir(out) for f in files)
 if bytes>=Int(c["artifact_budget_bytes"]);success=false;export_failure="artifact budget exceeded";end
 open(io->JSON3.pretty(io,(;success,stopped,completed_steps=completed,physical_time=x.w.time,failure,export_failure,resumable=false)),joinpath(out,"status.json"),"w")
 report=(schema="Interior.NativeBrineProducer",base_commit=get(ENV,"BRINE_BASE_COMMIT","unknown"),recipe=c,sources,
  recipe_sha256=bytes2hex(sha256(read(recipe))),runtime=(julia=string(VERSION),rocfft=string(AMDGPU.rocFFT.version()),library=AMDGPU.librocfft,optlevel=ENV["DEVICE_AMD_OPTLEVEL"]),preflight,warm_seconds,initial_setup,transfer_setup,after_upload,final_memory,
  success,failure,export_failure,resumable=false,completed_steps=completed,stopped,physical_time=x.w.time,precision="phase/forcing/permeability Float32; chemical $(producer_chemical_type(c)); FFT Float64/ComplexF64",frames,history,
  median_warm_step_seconds=isempty(times) ? nothing : median(times),evolution_seconds=(time_ns()-evolution_started)/1e9,
  wall_seconds=(time_ns()-NATIVE_BRINE_STARTED)/1e9,artifact_bytes=bytes,
  replay="Point samples at recorded replay shapes; scalar histories use full native fields. Final files retain full resolution. No migrated state or duplicated checkpoint.")
 open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w")
 println(JSON3.write((;event=success ? "completed" : "stopped",success,completed,physical_time=x.w.time,stopped,artifact_bytes=bytes,wall_seconds=report.wall_seconds)));flush(stdout)
 report
end
if abspath(PROGRAM_FILE)==@__FILE__
 fresh_output=!ispath(ARGS[2])
 try
  native_brine_producer(ARGS[1],ARGS[2];preflight_only="--preflight-only" in ARGS)
 catch err
  # A rejected reused output path must never modify its historical artifacts.
  if fresh_output
   mkpath(ARGS[2]);open(io->JSON3.pretty(io,(success=false,resumable=false,failure=sprint(showerror,err,catch_backtrace()))),joinpath(ARGS[2],"failure.json"),"w")
  end
  rethrow()
 end
end
