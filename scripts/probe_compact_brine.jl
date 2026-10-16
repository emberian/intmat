# Bounded whole-step parity/restart/timing gate. Never launches a producer.
include("device_cellular_setup.jl")
include(joinpath(@__DIR__,"..","studies","compact_brine_runtime.jl"))
using Serialization, Statistics
const COMPACT_GATE_RECIPE=TOML.parsefile(joinpath(@__DIR__,"..","experiments/art/brine-porous-384-f64.toml"))
function reference_upload(host)
 w,p=Interior.implicit_phase_transfer(host.w,host.p,backend)
 upload(a)=begin b=Interior.KernelAbstractions.zeros(backend,eltype(a),size(a)...);copyto!(b,a);b end
 (;w,p,a=upload(host.a),b=upload(host.b),scratch=upload(host.scratch),permeability=upload(host.permeability),runtime=Interior.ImplicitPhaseRuntime(p))
end
observation(x)=(;phi=Array(x.p.phi),interface=Array(x.p.interface),forcing=Array(x.p.chemical_forcing),a=Array(x.a),b=Array(x.b),mineral=Array(x.w.fields.deposit),permeability=Array(x.permeability),time=x.w.time)
function difference(x,y)
 a=observation(x);b=observation(y)
 a.time==b.time || error("clocks differ")
 errors=Dict(String(k)=>maximum(abs,Float64.(a[k]).-Float64.(b[k])) for k in keys(a) if k!=:time)
 maximum(values(errors))<=2e-5 || error("complete physical fields differ: $errors")
 errors
end
function runtime_memory(r)
 kind=="amd" || return nothing
 arrays=(r.forward.buffer,r.forward.workarea,r.inverse_real.buffer,r.inverse_real.workarea)
 (;forward_buffer=length(arrays[1])*sizeof(eltype(arrays[1])),forward_workspace=length(arrays[2])*sizeof(eltype(arrays[2])),inverse_buffer=length(arrays[3])*sizeof(eltype(arrays[3])),inverse_workspace=length(arrays[4])*sizeof(eltype(arrays[4])))
end
rows=[]
for (n,positive_feedback) in ((16,false),(32,false),(96,false),(32,true))
 local c=deepcopy(COMPACT_GATE_RECIPE)
 positive_feedback && (c["basal_reaction"]=.15)
 host_ref=phase_runtime_fixture(c,n;chemical_type=Float64);host_compact=compact_brine_fixture(c,n)
 if positive_feedback
  # Synthetic funded initial mineral activates feedback before the first step.
  # The ordinary recipe fixture remains independent and unchanged above.
  for host in (host_ref,host_compact)
   @. host.w.fields.deposit=.05*min(host.a,host.b)
   @. host.a-=host.w.fields.deposit/2
   @. host.b-=host.w.fields.deposit/2
  end
 end
 ref=reference_upload(host_ref);x=compact_brine_transfer(host_compact,backend)
 dt=.0025f0
 initial_mineral=sum(v->Float64(v),x.w.fields.deposit)
 initial_phase=sum(v->Float64(v),x.p.phi)
 initial=sum(v->Float64(v),x.a)+sum(v->Float64(v),x.b)+sum(v->Float64(v),x.w.fields.deposit)
 for _ in 1:20;coupled_runtime_step!(ref,c,dt,true);compact_brine_step!(x,c,dt);end
 errors=difference(ref,x)
 mass=(sum(v->Float64(v),x.a)+sum(v->Float64(v),x.b)+sum(v->Float64(v),x.w.fields.deposit)-initial)/initial
 phase_error=(sum(v->Float64(v),x.p.phi)-initial_phase)/length(x.p.phi)
 ranges=(phi=extrema(x.p.phi),interface=extrema(x.p.interface),permeability=extrema(x.permeability),forcing=extrema(x.p.chemical_forcing),a=extrema(x.a),b=extrema(x.b),mineral=extrema(x.w.fields.deposit))
 all(all(isfinite,r) for r in values(ranges)) && all(ranges[k][1]>=0 for k in (:a,:b,:mineral)) && abs(mass)<2e-6 && abs(phase_error)<2e-6 || error("unchanged brine finite/positivity/conservation guard failed")
 if positive_feedback
  sum(v->Float64(v),x.w.fields.deposit)>initial_mineral || error("positive fixture failed to exercise reaction")
  maximum(abs,x.p.chemical_forcing)>0 && ranges.interface[2]>0 && ranges.permeability[1]<ranges.permeability[2] || error("positive fixture failed to exercise mineral feedback, gradient or permeability")
 end
 payload=compact_brine_checkpoint(x,c)
 io=IOBuffer();serialize(io,payload);seekstart(io);y=compact_brine_restore(deserialize(io),c;backend)
 y.scratch===y.p.real_buffer || error("restart workspace alias lost")
 immediate_restart=difference(x,y)
 maximum(values(immediate_restart))==0 || error("restart physical state differs")
 for _ in 1:3;compact_brine_step!(x,c,dt);compact_brine_step!(y,c,dt);end
 restart=difference(x,y)
 maximum(values(restart))==0 || error("restart future differs")
 for _ in 1:3;coupled_runtime_step!(ref,c,dt,true);end
 timings=[]
 for block in 1:6
  compact=iseven(block)
  elapsed=@elapsed begin
   for _ in 1:5
    compact ? compact_brine_step!(x,c,dt) : coupled_runtime_step!(ref,c,dt,true)
   end
   Interior.KernelAbstractions.synchronize(backend)
  end
  push!(timings,(;compact,seconds_per_step=elapsed/5))
 end
 final_errors=difference(ref,x)
 push!(rows,(;n,fixture=positive_feedback ? "synthetic-funded-mineral-positive-reaction" : "original-recipe",errors,final_errors,restart,immediate_restart,initial_chemical_sum=initial,initial_mineral_sum=initial_mineral,final_mineral_sum=sum(v->Float64(v),x.w.fields.deposit),mass,phase_error,ranges,timings,memory=compact_brine_memory(n),runtime_memory=runtime_memory(x.runtime)))
 println(JSON3.write(last(rows)));flush(stdout)
 ref=nothing;x=nothing;y=nothing;payload=nothing;io=nothing;host_ref=nothing;host_compact=nothing;GC.gc()
end
open(io->JSON3.pretty(io,(;description="bounded complete physical brine step; Float64 chemicals unchanged; timings include validations and clocks, exclude field download/file IO",runtime=device_cellular_runtime(),source=COMPACT_BRINE_SOURCE,script_sha=bytes2hex(SHA.sha256(read(@__FILE__))),rows)),ARGS[1],"w")
