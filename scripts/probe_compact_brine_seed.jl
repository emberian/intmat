# Tiny correctness/resource gate for data-only initialization, never a producer.
include("device_cellular_setup.jl")
include(joinpath(@__DIR__,"..","studies","compact_brine_runtime.jl"))
using Serialization
c=TOML.parsefile(joinpath(@__DIR__,"..","experiments/art/brine-porous-384-f64.toml"))
fields(x)=(;phi=Array(x.p.phi),interface=Array(x.p.interface),forcing=Array(x.p.chemical_forcing),a=Array(x.a),b=Array(x.b),mineral=Array(x.w.fields.deposit),permeability=Array(x.permeability),time=x.w.time)
rows=[]
for n in (16,32)
 seed=compact_brine_seed(c,n)
 host_bytes=sum(sizeof,(seed.p.phi,seed.p.interface,seed.a,seed.b,seed.w.fields.deposit))
 host_bytes==compact_brine_seed_memory(n) || error("host seed byte accounting differs")
 x=compact_brine_seed_upload(seed,c,backend)
 ref=compact_brine_transfer(compact_brine_fixture(c,n),backend)
 fields(x)==fields(ref) || error("native seed initialization differs")
 x.scratch===x.p.real_buffer || error("native seed lost workspace alias")
 for _ in 1:5;compact_brine_step!(x,c,.0025f0);compact_brine_step!(ref,c,.0025f0);end
 fields(x)==fields(ref) || error("native seed future differs")
 cp=compact_brine_checkpoint(x,c);io=IOBuffer();serialize(io,cp);seekstart(io)
 y=compact_brine_restore(deserialize(io),c;backend)
 fields(x)==fields(y) && y.scratch===y.p.real_buffer || error("native seed checkpoint differs")
 for _ in 1:3;compact_brine_step!(x,c,.0025f0);compact_brine_step!(y,c,.0025f0);end
 fields(x)==fields(y) || error("native seed restart future differs")
 plans=(x.runtime.forward.buffer,x.runtime.forward.workarea,x.runtime.inverse_real.buffer,x.runtime.inverse_real.workarea)
 plan_bytes=map(a->length(a)*sizeof(eltype(a)),plans)
 push!(rows,(;n,host_seed_bytes=host_bytes,initialized_and_future_exact=true,checkpoint_and_future_exact=true,time=x.w.time,plan_bytes))
 println(JSON3.write(last(rows)));flush(stdout)
 x=nothing;ref=nothing;y=nothing;seed=nothing;cp=nothing;io=nothing;GC.gc()
end
open(io->JSON3.pretty(io,(;runtime=device_cellular_runtime(),source=COMPACT_BRINE_SOURCE,script_sha=bytes2hex(SHA.sha256(read(@__FILE__))),rows)),ARGS[1],"w")
