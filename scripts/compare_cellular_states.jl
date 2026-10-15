#!/usr/bin/env julia
# Compare trusted, bare Julia-serialized Arenas produced by the checkpoint probe.
# Both inputs are deserialized copies. Only proposal scratch values are zeroed;
# sizes, types, all other values, and the complete mutable alias graph must agree.
# Opt-in --materialize-genomes validates owned indexed banks and converts only
# CellularGenomeView weights/marks to their logical dense values on these copies.
using Interior, Serialization, JSON3
const compared_nodes=Ref(0)
const compared_float_bits=Ref(0)
const compared_arrays=Ref(0)
function differences(a,b,path="arena",seen=Set{Tuple{UInt,UInt}}(),forward=Dict{UInt,UInt}(),backward=Dict{UInt,UInt}())
 compared_nodes[]+=1
 typeof(a)==typeof(b)||return (path,string(typeof(a)),string(typeof(b)))
 if a isa Union{Float32,Float64}
  compared_float_bits[]+=1
  T=sizeof(a)==4 ? UInt32 : UInt64
  return reinterpret(T,a)==reinterpret(T,b) ? nothing : (path,a,b)
 end
 if a isa Union{Number,AbstractString,Symbol,Nothing,Type,Function}
  return isequal(a,b) ? nothing : (path,a,b)
 end
 if ismutabletype(typeof(a))
  key=(objectid(a),objectid(b))
  haskey(forward,key[1]) && forward[key[1]]!=key[2] && return (path,"aliasgraph forward mismatch",key)
  haskey(backward,key[2]) && backward[key[2]]!=key[1] && return (path,"aliasgraph reverse mismatch",key)
  forward[key[1]]=key[2];backward[key[2]]=key[1]
  key in seen&&return nothing;push!(seen,key)
 end
 if a isa AbstractArray
  compared_arrays[]+=1
  axes(a)==axes(b)||return (path*".axes",axes(a),axes(b))
  for i in eachindex(a,b);r=differences(a[i],b[i],path*"[$i]",seen,forward,backward);r===nothing||return r;end
 elseif a isa AbstractDict
  Set(keys(a))==Set(keys(b))||return (path*".keys",keys(a),keys(b))
  for k in keys(a);r=differences(a[k],b[k],path*"[$k]",seen,forward,backward);r===nothing||return r;end
 else
  for f in fieldnames(typeof(a));r=differences(getfield(a,f),getfield(b,f),path*"."*string(f),seen,forward,backward);r===nothing||return r;end
 end
 nothing
end
paths=filter(!=("--materialize-genomes"),ARGS)
length(paths)>=2 && iseven(length(paths)) || error("usage: compare_cellular_states.jl BASELINE_ARENA CANDIDATE_ARENA [BASELINE_ARENA CANDIDATE_ARENA ...] [--materialize-genomes]")
function main()
failed=false
for i in 1:2:length(paths)
 x=open(deserialize,paths[i]);y=open(deserialize,paths[i+1])
 x isa Arena && y isa Arena || error("inputs must be bare serialized Arenas, not checkpoint envelopes")
 canonicalized_lengths=Int[]
 for arena in (x,y), p in arena.processes
  if p isa CellularEcology
   push!(canonicalized_lengths,length(p.stencil.proposals))
   fill!(p.stencil.proposals,0.0)
   if "--materialize-genomes" in ARGS
    if p.weights isa Interior.CellularGenomeView || p.plastic_marks isa Interior.CellularGenomeView
     Interior._cellular_validate_storage(p) # Validate owned shared bank/unique slots before normalization.
     p.weights isa Interior.CellularGenomeView && (p.weights=Array(p.weights))
     p.plastic_marks isa Interior.CellularGenomeView && (p.plastic_marks=Array(p.plastic_marks))
    end
   end
  end
 end
 compared_nodes[]=0;compared_float_bits[]=0;compared_arrays[]=0
 seen=Set{Tuple{UInt,UInt}}();forward=Dict{UInt,UInt}();backward=Dict{UInt,UInt}()
 mismatch=differences(x,y,"arena",seen,forward,backward)
 failed |= mismatch!==nothing
 println(JSON3.write((;baseline=abspath(paths[i]),candidate=abspath(paths[i+1]),
     exact=mismatch===nothing,mismatch,compared_nodes=compared_nodes[],
     compared_float_bits=compared_float_bits[],compared_arrays=compared_arrays[],
     mutable_alias_pairs=length(forward),canonicalized_field="CellularEcology.stencil.proposals",
     genome_views_materialized="--materialize-genomes" in ARGS,
     canonicalized_lengths,canonicalization="zero values on deserialized copies only; retain length/type/alias checks")))
end
return failed ? 1 : 0
end
exit(main())
