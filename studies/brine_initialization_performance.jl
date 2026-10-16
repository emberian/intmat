# Source-derived old initializer oracle; never a second simulator.
include("implicit_phase_runtime.jl")
using SHA
const BRINE_INITIALIZER_REFERENCE_COMMIT="9c9f1a5"
function install_brine_initializer_reference!()
 root=normpath(joinpath(@__DIR__,".."))
 reference_file=get(ENV,"BRINE_INITIALIZER_REFERENCE_FILE",nothing)
 source=reference_file===nothing ? read(`git -C $root show $(BRINE_INITIALIZER_REFERENCE_COMMIT*":studies/brine_memory.jl")`,String) : read(reference_file,String)
 parsed=Meta.parseall(source)
 fn=only(filter(ex->ex isa Expr&&ex.head===:function&&ex.args[1] isa Expr&&ex.args[1].head===:call&&ex.args[1].args[1]===:brine_initialize!,parsed.args))
 fn=deepcopy(fn);fn.args[1].args[1]=:brine_initialize_reference!
 Core.eval(@__MODULE__,fn)
 bytes2hex(sha256(source))
end
function brine_initialization_equal(a,b)
 a.w.fields==b.w.fields&&a.a==b.a&&a.b==b.b&&a.p.phi==b.p.phi&&a.p.interface==b.p.interface&&a.p.chemical_potential==b.p.chemical_potential&&a.p.chemical_forcing==b.p.chemical_forcing
end
function benchmark_brine_initialization(output;resolutions=(96,192),repeats=3)
 reference_sha=install_brine_initializer_reference!()
 c=TOML.parsefile(joinpath(@__DIR__,"..","experiments/art/brine-porous-384-f64.toml"));rows=Any[]
 warm=phase_runtime_fixture(c,8;chemical_type=Float64)
 Base.invokelatest(brine_initialize_reference!,warm.w,warm.p,warm.a,warm.b,c)
 for n in resolutions
  allocation_and_specialized_seconds=@elapsed a=phase_runtime_fixture(c,n;chemical_type=Float64)
  b=phase_runtime_fixture(c,n;chemical_type=Float64)
  times=Any[]
  for repeat in 1:repeats
   # Alternate ordering; both timings include all fields and phase_observe!.
   old=0.;new=0.
   for reference in (isodd(repeat),!isodd(repeat))
    seconds=@elapsed if reference
     Base.invokelatest(brine_initialize_reference!,b.w,b.p,b.a,b.b,c)
    else
     brine_initialize!(a.w,a.p,a.a,a.b,c)
    end
    if reference;old=seconds;else;new=seconds;end
   end
   brine_initialization_equal(a,b)||error("initializer changed physical fields at$n")
   push!(times,(;reference_seconds=old,specialized_seconds=new))
  end
  push!(rows,(;n,chemical_storage="Float64",allocation_and_specialized_seconds,times,
   speedup=median(getproperty.(times,:reference_seconds))/median(getproperty.(times,:specialized_seconds)),bitwise_physical_parity=true))
  println(JSON3.write(last(rows)));flush(stdout);a=nothing;b=nothing;GC.gc()
 end
 report=(;reference_commit=BRINE_INITIALIZER_REFERENCE_COMMIT,reference_file=get(ENV,"BRINE_INITIALIZER_REFERENCE_FILE",nothing),reference_sha,current_sha=bytes2hex(sha256(read(joinpath(@__DIR__,"brine_memory.jl")))),threads=Threads.nthreads(),rows,
  timing="Complete brine_initialize! including all chemical/phase fields and phase_observe!; planning/allocation excluded from paired medians and reported separately for specialized factory.")
 open(io->JSON3.pretty(io,report),output,"w");report
end
if abspath(PROGRAM_FILE)==@__FILE__;benchmark_brine_initialization(ARGS[1]);end
