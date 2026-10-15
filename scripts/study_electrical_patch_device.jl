#!/usr/bin/env julia
# Fresh, bounded fixture study. Never opens a production run for writing.
using Interior, SHA, Serialization, Statistics
using Interior: JSON3, KernelAbstractions
include(joinpath(@__DIR__,"..","studies","electrical_patch_kernel.jl"))
include(joinpath(@__DIR__,"..","studies","branch_metal_multigrid.jl"))
function electrical_patch_fixture(geometry_path,prefix,events)
 g=JSON3.read(read(geometry_path,String));n=Int(g.shape[1])
 all(==(n),g.shape) || error("noncubic geometry")
 mask=falses(n,n,n);sites=Dict{Int,CartesianIndex{3}}()
 for (xyz,birth) in zip(g.positions,g.birth)
  I=CartesianIndex(Tuple(Int.(xyz)));b=Int(birth)
  if b<=prefix;mask[I]=true;else;sites[b]=I;end
 end
 next=[sites[b] for b in prefix+1:prefix+events]
 (n=n,mask=mask,sites=next,prefix=prefix,geometry_sha256=bytes2hex(sha256(read(geometry_path))))
end
function electrical_patch_device_study(backend_name,geometry_path,output;prefix=4096,events=9)
 ispath(output)&&error("refusing existing output")
 mkpath(output)
 hardware=Dict{String,Any}()
 upload=if backend_name=="amd"
  @eval using AMDGPU
  device=AMDGPU.device();props=AMDGPU.HIP.properties(device)
  hardware["device"]=AMDGPU.HIP.name(device)
  hardware["shared_memory_per_block"]=Int(props.sharedMemPerBlock)
  hardware["max_threads_per_block"]=Int(props.maxThreadsPerBlock)
  hardware["shared_memory_per_block"]>=36*1024 || error("device LDS limit insufficient for patch")
  hardware["max_threads_per_block"]>=256 || error("device thread limit insufficient for patch")
  println(JSON3.write(hardware));flush(stdout)
  AMDGPU.ROCArray
 elseif backend_name=="cpu";identity
 else;error("backend must be cpu or amd")
 end
 fixture=electrical_patch_fixture(geometry_path,prefix,events)
 n=fixture.n;tol=1.6e-5;group=(32,4,2)
 sources=Dict{String,String}()
 for f in ("studies/electrical_local_correction.jl","studies/electrical_patch_kernel.jl","studies/branch_solver_prototype.jl","studies/branch_metal_multigrid.jl","scripts/study_electrical_patch_device.jl")
  bytes=read(joinpath(@__DIR__,"..",f));sources[f]=bytes2hex(sha256(bytes));write(joinpath(output,basename(f)),bytes)
 end
 # Compile both paths without mutating the fixture. No native timing includes compilation.
 tiny=process(17,specimen(17));warm=device_levels(hierarchy(tiny.potential,tiny.occupied),upload)
 device_solve!(warm;tol,workgroup=group,norm_workspace=DeviceMGNorm(first(warm)))
 expected=Array(first(warm).u);electrical_patch!(expected,tiny.occupied,CartesianIndex(8,8,8);radius=8,sweeps=16)
 electrical_device_patch8!(first(warm).u,first(warm).occupied,CartesianIndex(8,8,8))
 reinterpret(UInt32,vec(Array(first(warm).u)))==reinterpret(UInt32,vec(expected)) || error("device patch differs from scalar CPU patch")
 warm=nothing;GC.gc()
 p=process(n,fixture.mask)
 a=device_levels(hierarchy(copy(p.potential),copy(p.occupied)),upload)
 norm_a=DeviceMGNorm(first(a));initial=device_solve!(a;tol,workgroup=group,norm_workspace=norm_a)
 initial_u=Array(first(a).u)
 # This cold reconstruction matches geometry/operator, not production's last
 # finite-tolerance iterate. Preserve that distinction and the exact generated field.
 serialize(joinpath(output,"fixture.jls"),merge(fixture,(potential=initial_u,tolerance=tol,initial_residual=initial.residual)))
 b=device_levels(hierarchy(copy(initial_u),copy(p.occupied)),upload)
 norm_b=DeviceMGNorm(first(b));patch_norm=similar(first(b).u,Float32,1)
 prior=initial.residual;rows=Any[];backend=KernelAbstractions.get_backend(first(a).u)
 timed(f)=begin
  KernelAbstractions.synchronize(backend);start=time_ns();stats=f();KernelAbstractions.synchronize(backend)
  merge(stats,(seconds=(time_ns()-start)*1e-9,))
 end
 for (offset,I) in enumerate(fixture.sites)
  Interior.electrical_mark!(p,I);p.events+=1;device_mark!(a,I);device_mark!(b,I)
  baseline=()->device_solve!(a;tol,workgroup=group,norm_workspace=norm_a)
  candidate=()->begin
   patch=electrical_device_patch8!(first(b).u,first(b).occupied,I;norm=patch_norm)
   bound=max(prior,patch.residual)
   if bound<=tol;(cycles=0,residual=bound,patch_residual=patch.residual)
   else
    stats=device_solve!(b;tol,workgroup=group,norm_workspace=norm_b)
    merge(stats,(patch_residual=patch.residual,))
   end
  end
  ta,tb=if isodd(offset);(timed(baseline),timed(candidate))
   else;cb=timed(candidate);ba=timed(baseline);(ba,cb)
  end
  prior=tb.residual
  au=Array(first(a).u);bu=Array(first(b).u)
  p.potential.=au;ra=Interior.potential_residual(p)
  p.potential.=bu;rb=Interior.potential_residual(p)
  max(ra,rb)<=2e-5 || error("independent global residual failed")
  all(==(1f0),bu[p.occupied]) || error("conductor changed")
  boundary=all(iszero,bu[1,:,:])&&all(iszero,bu[end,:,:])&&all(iszero,bu[:,1,:])&&all(iszero,bu[:,end,:])&&all(iszero,bu[:,:,1])&&all(iszero,bu[:,:,end])
  boundary||error("ground boundary changed")
  push!(rows,(event=prefix+offset,site=Tuple(I),baseline=ta,candidate=tb,
   independent_baseline_residual=ra,independent_candidate_residual=rb,
   maximum_potential_difference=maximum(abs.(au.-bu)),speedup=ta.seconds/tb.seconds))
  println(JSON3.write(last(rows)));flush(stdout)
 end
 report=(scope="complete harmonic solve including patch dispatch and synchronization; recorded geometry sequence, not coupled growth parity",
  backend=backend_name,hardware=hardware,julia=string(VERSION),threads=Threads.nthreads(),n=n,prefix=prefix,events=events,
  geometry_sha256=fixture.geometry_sha256,source_sha256=sources,tolerance=tol,
  initial_state="cold reconstructed exact occupied prefix, not saved production potential",lds_bytes=35319,
  median_paired_speedup=median(r.speedup for r in rows),rows=rows)
 open(io->JSON3.pretty(io,report),joinpath(output,"results.json"),"w")
 report
end
if abspath(PROGRAM_FILE)==@__FILE__
 length(ARGS)>=3||error("usage: cpu|amd GEOMETRY_JSON FRESH_OUTPUT [PREFIX=4096] [EVENTS=9]")
 electrical_patch_device_study(ARGS[1],ARGS[2],ARGS[3];prefix=length(ARGS)>3 ? parse(Int,ARGS[4]) : 4096,events=length(ARGS)>4 ? parse(Int,ARGS[5]) : 9)
end
