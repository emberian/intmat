using Interior,Serialization,Statistics,SHA
using Interior:JSON3
include(joinpath(@__DIR__,"..","studies","lattice_split_bank.jl"))
function split_convert(w,p,mod,kind)
 weights=kind==:split ? mod.lattice_splitbank() : permutedims(p.weights,(2,3,1))
 q=mod.LatticeEcology((name==:weights ? weights : deepcopy(getfield(p,name)) for name in fieldnames(typeof(p)))...)
 v=deepcopy(w);Interior.validate(v,q);v,q
end
function split_window(w,p,n)
 @timed for _ in 1:n;step!(w,p,.1);end
end
function split_main(out)
 ispath(out)&&error("fresh output required");mkpath(out)
 mods=install_lattice_split_variants(joinpath(out,"variants"))
 Base.invokelatest(split_runs,out,mods)
end
function split_runs(out,mods)
 coefficient=lattice_split_coefficient_audit(mods.split);controller=lattice_split_controller_audit(mods.split)
 rows=Any[]
 for dims in ((32,32,1),(64,64,1),(32,32,32))
  w=LatticeEcologyStudy.lattice_world(Grid(dims...;extent=dims));p=LatticeEcologyStudy.LatticeEcology(w;occupancy=.5)
  for _ in 1:64;step!(w,p,.1);end
  initial=LatticeEcologyStudy.lattice_statistics(w,p)
  states=Dict(:baseline=>(w,p),:packed=>split_convert(w,p,mods.packed,:packed),:split=>split_convert(w,p,mods.split,:split))
  n=dims[3]==1 ? 100 : 20
  for arm in (:baseline,:packed,:split);x,y=deepcopy(states[arm]);split_window(x,y,2);end
  samples=Any[];ends=Dict{Symbol,Any}()
  for repeat in 1:3
   for arm in (isodd(repeat) ? (:baseline,:packed,:split) : (:split,:packed,:baseline))
    x,y=deepcopy(states[arm]);GC.gc();t=split_window(x,y,n)
    push!(samples,(repeat=repeat,arm=arm,seconds=t.time,bytes=t.bytes,gc_seconds=t.gctime));ends[arm]=(x,y)
   end
  end
  errors=Dict{String,Any}()
  a,b=ends[:baseline]
  for arm in (:packed,:split)
   x,y=ends[arm];differences=Dict{String,Any}()
   for name in fieldnames(typeof(b))
    name==:weights&&continue
    v=getfield(b,name);z=getfield(y,name)
    if v isa AbstractArray
     differences[string(name)]=eltype(v)<:AbstractFloat ? (exact=isequal(v,z),maxabs=maximum(abs,Float64.(v).-Float64.(z);init=0.)) : (exact=isequal(v,z),changed=count(!iszero,v.!=z))
    end
   end
   mod=arm==:split ? mods.split : mods.packed
   stat=mod.lattice_statistics(x,y)
   errors[string(arm)]=(arrays=differences,statistics=stat,mass_residual=stat.material-initial.material,energy_residual=stat.energy_accounted-initial.energy_accounted)
   serialize(joinpath(out,"state-$(join(dims,'x'))-$arm.jls"),(world=x,process=y))
  end
  row=(dims=dims,warmup=64,steps=n,initial=initial,samples=samples,errors=errors,
   baseline_final=LatticeEcologyStudy.lattice_statistics(a,b),speed_vs_baseline=median(s.seconds for s in samples if s.arm==:baseline)/median(s.seconds for s in samples if s.arm==:split),
   speed_vs_packed=median(s.seconds for s in samples if s.arm==:packed)/median(s.seconds for s in samples if s.arm==:split))
  push!(rows,row);println(JSON3.write(row));flush(stdout)
  open(io->JSON3.pretty(io,(coefficient=coefficient,controller=controller,rows=rows,core=Interior.arena_sources(),julia=string(VERSION),threads=Threads.nthreads())),joinpath(out,"results.json"),"w")
 end
end
if abspath(PROGRAM_FILE)==@__FILE__;split_main(only(ARGS));end
