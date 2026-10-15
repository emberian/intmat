# Isolated order-preserving contact-search factorization. No core method changes.
# Cache only a subsequence of the original bin walk. Any changed bin/list order,
# identity, radius or exceeded skin invalidates reuse for the rest of this call.
function dc_contact_definition(x,name)
 x isa Expr||return nothing
 if x.head===:function
  signature=x.args[1]
  signature isa Expr&&signature.head===:call&&signature.args[1]===name&&return deepcopy(x)
 end
 for arg in x.args
  result=dc_contact_definition(arg,name);result===nothing||return result
 end
 nothing
end
function dc_contact_source_install!()
 path=joinpath(@__DIR__,"..","src","device_cellular_kernels.jl");bytes=read(path)
 bytes2hex(SHA.sha256(bytes))==DEVICE_CELLULAR_SOURCE_PIN.kernels||error("loaded contact source differs from disk")
 tree=Meta.parseall(String(copy(bytes)))
 force=dc_contact_definition(tree,:dc_forces_at!);force===nothing&&error("force source missing")
 force.args[1].args[1]=:dc_cached_forces_at!;push!(force.args[1].args,:cache)
 # Preserve the exact arithmetic of the original neighbour contribution body.
 replaced=Ref(0)
 function rewrite!(x;allpairs=false)
  x isa Expr||return x
  if x.head===:for&&x.args[1] isa Expr&&x.args[1].head===:block&&
    occursin("dc_search_count",string(x.args[1]))
   contribution=nothing
   function find_contribution(y)
    y isa Expr||return
    if y.head===:if&&y.args[1]==:(j!=i);contribution=deepcopy(y);return;end
    for arg in y.args;find_contribution(arg);end
   end
   find_contribution(x);contribution===nothing&&error("force contribution source missing")
   replaced[]+=1
   if allpairs
    return :(for j in 1:length(s.active)
      if s.active[j]==1;$contribution;end
     end)
   end
   return :(for slot in 1:cache.count[i]
      j=cache.indices[i,slot]
      $contribution
    end)
  end
  for i in eachindex(x.args);x.args[i]=rewrite!(x.args[i];allpairs);end
  x
 end
 allpairs=dc_contact_definition(tree,:dc_forces_at!)
 allpairs.args[1].args[1]=:dc_allpairs_forces_at!
 rewrite!(allpairs;allpairs=true);replaced[]==1||error("expected one all-pairs force bin walk")
 Core.eval(@__MODULE__,allpairs)
 replaced[]=0
 rewrite!(force);replaced[]==1||error("expected one cached force bin walk")
 insert!(force.args[2].args,1,:(if cache.invalid[1]!=0||cache.overflow[i]!=0
   s.active[i]==1&&(cache.fallbacks[i]+=Int32(1))
   return cache.strict_order ? dc_forces_at!(s,e,q,g,dt,i) : dc_allpairs_forces_at!(s,e,q,g,dt,i)
  end))
 insert!(force.args[2].args,2,:(s.active[i]==1&&(cache.hits[i]+=Int32(1))))
 Core.eval(@__MODULE__,force)
 geometry=dc_contact_definition(tree,:dc_geometry!);geometry===nothing&&error("geometry source missing")
 geometry.args[1].args[1]=:dc_cached_geometry!;push!(geometry.args[1].args,:cache)
 push!(geometry.args[2].args,:(dc_contact_validate_at!(s,e,g,cache,i)))
 Core.eval(@__MODULE__,Expr(:macrocall,Symbol("@kernel"),LineNumberNode(0),geometry))
 commit=dc_contact_definition(tree,:dc_commit_positions!)
 commit.args[1].args[1]=:dc_cached_commit_positions!
 append!(commit.args[1].args,[:g,:cache])
 push!(commit.args[2].args,:(dc_contact_validate_at!(s,e,g,cache,i)))
 Core.eval(@__MODULE__,Expr(:macrocall,Symbol("@kernel"),LineNumberNode(0),commit))
 (;source=String(copy(bytes)),force=string(force),allpairs=string(allpairs),geometry=string(geometry),commit=string(commit))
end
# Wrapped coordinates and only adjacent-period subtraction: padding dominates
# Float32 coordinate/subtraction rounding; rejecting extra candidates is optional.
@inline dc_contact_padding(g)=64f0*eps(Float32)*max(g.extent[1],g.extent[2],g.extent[3])
@inline function dc_contact_validate_at!(s,e,g,cache,i)
 bad=s.active[i]!=cache.active[i]||s.ids[i]!=cache.ids[i]||s.radii[i]!=cache.radii[i]
 if cache.strict_order;bad|=s.cell_bin[i]!=cache.bin[i]||s.next[i]!=cache.next[i];end
 if s.active[i]==1
  padding=dc_contact_padding(g)
  for d in 1:Int(g.nd)
   delta=s.position[i,d]-cache.origin[i,d];L=g.extent[d];delta-=floor(delta/L+.5f0)*L
   bad|=abs(delta)>cache.skin/2f0-padding
  end
 end
 bad&&(KernelAbstractions.@atomic cache.invalid[1]|=Int32(1))
 nothing
end
@inline function dc_contact_search_range(s,e,g,cache,i,d)
 d>g.nd&&return Int32(0)
 width=g.extent[d]/Float32(g.bins[d])
 radius=s.radii[i]+reinterpret(Float32,e.metadata[7])+cache.skin+dc_contact_padding(g)
 padding=8f0*eps(Float32)*Float32(g.bins[d])
 Int32(min(Float32(g.bins[d]),ceil(radius/width+padding)))
end
@kernel function dc_build_contact_candidates!(s,e,g,cache)
 i=@index(Global,Linear)
 cache.count[i]=0;cache.overflow[i]=0;cache.hits[i]=0;cache.fallbacks[i]=0
 cache.active[i]=s.active[i];cache.ids[i]=s.ids[i];cache.radii[i]=s.radii[i]
 cache.bin[i]=s.cell_bin[i];cache.next[i]=s.next[i]
 for d in 1:3;cache.origin[i,d]=s.position[i,d];end
 if e.status[1]==0&&s.active[i]==1
  key=dc_bin(s,g,i);range=ntuple(d->cache.strict_order ? dc_search_range(s,e,g,i,d) : dc_contact_search_range(s,e,g,cache,i,d),3)
  count=Int32(0)
  for dz in 0:dc_search_count(g.bins[3],range[3])-1,dy in 0:dc_search_count(g.bins[2],range[2])-1,dx in 0:dc_search_count(g.bins[1],range[1])-1
   bin=(dc_neighbor_bin(key[1],dx,g.bins[1],range[1]),dc_neighbor_bin(key[2],dy,g.bins[2],range[2]),dc_neighbor_bin(key[3],dz,g.bins[3],range[3]))
   j=e.heads[dc_index(bin...,g.bins)]
   while j>0
    if j!=i
     included=true;bound=s.radii[i]+s.radii[j]+cache.skin+dc_contact_padding(g)
     for d in 1:Int(g.nd);included&=abs(dc_delta(s,g,i,j,d))<=bound;end
     if included
      count+=Int32(1)
      if count<=size(cache.indices,2);cache.indices[i,count]=j
      else;cache.overflow[i]=1;end
     end
    end
    j=s.next[j]
   end
  end
  cache.count[i]=min(count,Int32(size(cache.indices,2)))
 end
end
const DEVICE_CONTACT_CANDIDATE_SOURCE=dc_contact_source_install!()
@kernel function dc_cached_forces!(s,e,q,g,dt,cache)
 i=@index(Global,Linear)
 dc_cached_forces_at!(s,e,q,g,dt,i,cache)
end
function device_contact_candidate_cache(p;capacity=64,skin=.025f0,strict_order=true)
 C=Int(p.geometry.capacity);g=p.geometry
 capacity>0&&isfinite(skin)&&4dc_contact_padding(g)<skin&&(!strict_order||skin<=.15f0)||throw(ArgumentError("positive capacity and resolvable physical skin required"))
 backend=KernelAbstractions.get_backend(p.state.active)
 z(T,dims...)=KernelAbstractions.zeros(backend,T,dims...)
 (indices=z(Int32,C,capacity),count=z(Int32,C),overflow=z(Int32,C),invalid=z(Int32,1),
  origin=z(Float32,C,3),radii=z(Float32,C),active=z(Int32,C),ids=z(Int32,C),bin=z(Int32,C),next=z(Int32,C),
  hits=z(Int32,C),fallbacks=z(Int32,C),skin=Float32(skin),strict_order=Bool(strict_order))
end
function device_candidate_mechanics!(p,delta,cache)
 # Strict mode retains bin rebuild/order; coverage mode retains complete pairs
 # and frozen contribution order, avoiding30 geometry launches. Invalid/overflow
 # rows use complete deterministic all-pairs forces, never stale bins.
 stage=p.ownership.stages.mechanics;s=stage.state;e=stage.environment;q=stage.parameters;g=stage.geometry
 C=Int(p.geometry.capacity);substeps=Int(p.parameters.mechanical_substeps);dt=Float32(delta)/Float32(substeps)
 fill!(cache.invalid,Int32(0))
 device_cellular_launch!(p,dc_clear_heads!,e;ndrange=length(e.heads))
 device_cellular_launch!(p,dc_geometry!,s,e,q,g;ndrange=C)
 device_cellular_launch!(p,dc_build_contact_candidates!,s,e,g,cache;ndrange=C)
 for step in 1:substeps
  if step>1&&cache.strict_order
   device_cellular_launch!(p,dc_clear_heads!,e;ndrange=length(e.heads))
   device_cellular_launch!(p,dc_cached_geometry!,s,e,q,g,cache;ndrange=C)
  end
  device_cellular_launch!(p,dc_cached_forces!,s,e,q,g,dt,cache;ndrange=C)
  device_cellular_launch!(p,dc_move!,s,e,q,g,dt;ndrange=C)
  if cache.strict_order
   device_cellular_launch!(p,dc_commit_positions!,s,e;ndrange=C)
  else
   device_cellular_launch!(p,dc_cached_commit_positions!,s,e,g,cache;ndrange=C)
  end
 end
 device_cellular_launch!(p,dc_clear_heads!,e;ndrange=length(e.heads))
 device_cellular_launch!(p,dc_geometry!,s,e,q,g;ndrange=C)
 nothing
end

struct DeviceContactCandidateCarrier{P,C}
 process::P
 cache::C
end
@inline Base.getproperty(p::DeviceContactCandidateCarrier,key::Symbol)=key in (:process,:cache) ? getfield(p,key) : getproperty(getfield(p,:process),key)
function Base.setproperty!(p::DeviceContactCandidateCarrier,key::Symbol,value)
 key===:attempted_steps||throw(ArgumentError("only attempt bookkeeping can be forwarded"))
 setproperty!(p.process,key,value)
end
device_cellular_validate(w,p::DeviceContactCandidateCarrier)=device_cellular_validate(w,p.process)
device_cellular_timestep_limit(p::DeviceContactCandidateCarrier)=device_cellular_timestep_limit(p.process)
device_cellular_status!(p::DeviceContactCandidateCarrier)=device_cellular_status!(p.process)
device_support_descriptor(p::DeviceContactCandidateCarrier)=device_support_descriptor(p.process)
device_cellular_mechanics!(p::DeviceContactCandidateCarrier,dt,::Val{mode}) where mode=device_candidate_mechanics!(p.process,dt,p.cache)
function device_candidate_step!(w,p::DeviceContactCandidateCarrier,dt;anabolism=nothing)
 device_cellular_validate(w,p);delta=Float32(dt)
 isfinite(delta)&&0<delta<=device_cellular_timestep_limit(p)||throw(ArgumentError("invalid candidate timestep"))
 p.attempted_steps+=1
 device_cellular_enqueue!(w,p,delta,:projected,anabolism,UInt32(p.attempted_steps))
 device_cellular_status!(p);p
end

arena_process_supported(::DeviceContactCandidateCarrier)=false
