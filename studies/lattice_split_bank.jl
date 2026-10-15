# Source-derived arithmetic experiment; never replaces the frozen lattice module.
using SHA
isdefined(@__MODULE__,:LatticeEcologyStudy)||include("lattice_ecology.jl")
const LATTICE_SPLIT_BASE_SOURCE=read(joinpath(@__DIR__,"lattice_ecology.jl"),String)
bytes2hex(sha256(LATTICE_SPLIT_BASE_SOURCE))==LatticeEcologyStudy.MODEL_SHA||error("loaded lattice source differs from disk")
function lattice_split_source()
 source=LATTICE_SPLIT_BASE_SOURCE
 lo=findfirst("function lattice_genotypes(",source);hi=findnext("    bank=Array{Float32}",source,last(lo))
 lo!==nothing&&hi!==nothing||error("genotype component markers changed")
 components=replace(source[first(lo):first(hi)-1],"function lattice_genotypes("=>"function lattice_split_components(")*"    (;base,basis)\nend\n"
 helpers=raw"""
struct LatticeSplitBank{A} <: AbstractArray{Float32,3}
 low::A
 high::A
 loci::Int
 split::Int
end
Base.size(b::LatticeSplitBank)=(1<<b.loci,INPUTS,OUTPUTS)
Base.IndexStyle(::Type{<:LatticeSplitBank})=IndexCartesian()
Base.dataids(b::LatticeSplitBank)=(Base.dataids(b.low)...,Base.dataids(b.high)...)
Base.sizeof(b::LatticeSplitBank)=sizeof(b.low)+sizeof(b.high)
KernelAbstractions.get_backend(::LatticeSplitBank)=CPU()
@inline function Base.getindex(b::LatticeSplitBank,code::Integer,col::Integer,output::Integer)
 @boundscheck checkbounds(b,code,col,output)
 c=Int(code)-1;mask=(1<<b.split)-1
 @inbounds b.low[col,output,(c&mask)+1]+b.high[col,output,(c>>b.split)+1]
end
function lattice_splitbank(;loci=12,seed=911,strength=.065)
 base,basis=lattice_split_components(;loci,seed,strength)
 split=div(loci,2);low=Array{Float32}(undef,INPUTS,OUTPUTS,1<<split)
 high=Array{Float32}(undef,INPUTS,OUTPUTS,1<<(loci-split))
 for output in 1:OUTPUTS,input in 1:INPUTS
  for code in 0:(1<<split)-1
   value=base[input,output]
   for bit in 0:split-1;value+=(((code>>bit)&1)==1 ? 1f0 : -1f0)*basis[bit+1,input,output];end
   low[input,output,code+1]=value
  end
  for code in 0:(1<<(loci-split))-1
   value=0f0
   for bit in 0:loci-split-1;value+=(((code>>bit)&1)==1 ? 1f0 : -1f0)*basis[split+bit+1,input,output];end
   high[input,output,code+1]=value
  end
 end
 all(isfinite,low)&&all(isfinite,high)||throw(ArgumentError("split allele table overflows Float32"))
 LatticeSplitBank(low,high,loci,split)
end
"""
 constructor="    weights=lattice_genotypes(;loci,seed=genotype_seed);bank=z(Float32,size(weights)...);copyto!(bank,weights)"
 length(findall(constructor,source))==1||error("lattice constructor marker changed")
 variant=replace(source,"module LatticeEcologyStudy"=>"module LatticeSplitBankStudy",
  constructor=>"    w.backend isa CPU||throw(ArgumentError(\"split-bank study is CPU-only\"));bank=lattice_splitbank(;loci,seed=genotype_seed)",
  "mutable struct LatticeEcology{"=>components*helpers*"\nmutable struct LatticeEcology{",
  "function lattice_checkpoint(w,p,path)\n"=>"function lattice_checkpoint(w,p,path)\n    throw(ArgumentError(\"split-bank checkpoint adapter is not implemented; preserve study observations/source\"))\n",
  "function restore_lattice(path;backend=CPU())\n"=>"function restore_lattice(path;backend=CPU())\n    throw(ArgumentError(\"split-bank checkpoint adapter is not implemented\"))\n")
 variant
end
function lattice_packed_control_source()
 replace(LATTICE_SPLIT_BASE_SOURCE,"module LatticeEcologyStudy"=>"module LatticePackedBankControl",
  "weights=lattice_genotypes(;loci,seed=genotype_seed)"=>"weights=permutedims(lattice_genotypes(;loci,seed=genotype_seed),(2,3,1))",
  "size(p.weights)==(1<<p.loci,INPUTS,OUTPUTS)"=>"size(p.weights)==(INPUTS,OUTPUTS,1<<p.loci)","bank[code,col,output]"=>"bank[col,output,code]")
end
function install_lattice_split_variants(output)
 isdir(output)||mkpath(output)
 for (file,source) in (("split.jl",lattice_split_source()),("packed.jl",lattice_packed_control_source()))
  path=joinpath(output,file);ispath(path)&&error("fresh variant files required")
  write(path,source);Base.include(@__MODULE__,path)
 end
 (split=Base.invokelatest(getproperty,@__MODULE__,:LatticeSplitBankStudy),packed=Base.invokelatest(getproperty,@__MODULE__,:LatticePackedBankControl))
end
# Independent Float64 target and forward error bound for the quantized additive
# coefficients. Both are separate from either Float32 accumulation path.
function lattice_split_coefficient_audit(mod;loci=12,seed=911,strength=.065)
 base,basis=mod.lattice_split_components(;loci,seed,strength)
 original=LatticeEcologyStudy.lattice_genotypes(;loci,seed,strength)
 split=mod.lattice_splitbank(;loci,seed,strength)
 u=Float64(eps(Float32))/2;gamma(n)=n*u/(1-n*u)
 maximum_delta=0.;old_error=0.;new_error=0.;bound=0.;violations=0
 for output in axes(original,3),input in axes(original,2),code in axes(original,1)
  target=Float64(base[input,output]);magnitude=abs(target)
  for bit in 0:loci-1
   term=Float64(basis[bit+1,input,output]);target+=(((code-1)>>bit)&1)==1 ? term : -term;magnitude+=abs(term)
  end
  local_bound=(gamma(loci)+gamma(max(split.split,loci-split.split)+1))*magnitude
  a=Float64(original[code,input,output]);b=Float64(split[code,input,output])
  maximum_delta=max(maximum_delta,abs(a-b));old_error=max(old_error,abs(a-target));new_error=max(new_error,abs(b-target));bound=max(bound,local_bound)
  violations+=abs(a-b)>local_bound
 end
 (;loci,seed,strength,coefficients=length(original),maximum_delta,old_error,new_error,conservative_delta_bound=bound,violations,
  full_bytes=sizeof(original),split_bytes=sizeof(split))
end
function lattice_split_controller_audit(mod;samples=256,seed=721)
 old=LatticeEcologyStudy.lattice_genotypes();bank=mod.lattice_splitbank();rng=LatticeEcologyStudy.MersenneTwister(seed)
 gamma14=(14Float64(eps(Float32))/2)/(1-14Float64(eps(Float32))/2)
 maximum_dot_difference=0.;maximum_bound=0.;hidden_error=0.;action_error=0.;violations=0
 for sample in 1:samples
  code=sample<=2 ? (sample==1 ? 1 : 4096) : rand(rng,1:4096)
  inputs=sample==1 ? ntuple(_->1f0,14) : sample==2 ? ntuple(i->isodd(i) ? -1f0 : 1f0,14) : Tuple(2f0.*rand(rng,Float32,14).-1f0)
  for output in 1:11
   a=LatticeEcologyStudy.network_output(old,code,inputs,output)
   b=mod.network_output(bank,code,inputs,output)
   difference=0.;magnitude=0.
   for col in 1:14
    x=Float64(old[code,col,output]);y=Float64(bank[code,col,output])
    difference+=abs(x-y);magnitude+=abs(x)+abs(y)
   end
   bound=difference+gamma14*magnitude
   maximum_dot_difference=max(maximum_dot_difference,abs(Float64(a)-b));maximum_bound=max(maximum_bound,bound)
   violations+=abs(Float64(a)-b)>bound
   if output<=4;hidden_error=max(hidden_error,abs(Float64(tanh(a))-tanh(b)))
   else;action_error=max(action_error,abs(Float64(LatticeEcologyStudy.gate(a))-mod.gate(b)));end
  end
 end
 (;samples,neurons=samples*11,maximum_dot_difference,conservative_dot_bound=maximum_bound,violations,
  hidden_error,action_error,real_tanh_lipschitz_bound=maximum_bound,real_logistic_lipschitz_bound=maximum_bound/4,
  scope="bounded |input|<=1; FMA forward rounding bound; real-map Lipschitz bounds exclude transcendental implementation rounding")
end
