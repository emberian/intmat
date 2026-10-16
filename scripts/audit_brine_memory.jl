# Independent inventory audit of saved fields, including both stoichiometric invariants.
using JSON3, SHA
function audit_brine(run,out)
 r=JSON3.read(read(joinpath(run,"run.json"),String));shape=Tuple(Int.(r.shape));volume=prod(Float64.(r.extent)./shape);cells=prod(shape)
 readfield(frame,name)=begin
  bytes=read(joinpath(run,String(frame.directory),name*".f32"));length(bytes)==4cells||error("Field byte length");(reinterpret(Float32,bytes),bytes2hex(sha256(bytes)))
 end
 a0,_=readfield(first(r.frames),"a");b0,_=readfield(first(r.frames),"b");phi0,_=readfield(first(r.frames),"phi");initial_a=sum(Float64,a0)*volume;initial_b=sum(Float64,b0)*volume;initial_phi=sum(Float64,phi0);previous=zeros(Float32,cells);rows=Any[]
 for frame in r.frames
  a,asha=readfield(frame,"a");b,bsha=readfield(frame,"b");m,msha=readfield(frame,"mineral");phi,phisha=readfield(frame,"phi")
  ma=sum(Float64,a)*volume;mb=sum(Float64,b)*volume;mm=sum(Float64,m)*volume
  ae=ma+mm/2-initial_a;be=mb+mm/2-initial_b;increment=minimum(m.-previous);phase_mean=(sum(Float64,phi)-initial_phi)/cells
  minimum(a)>=0&&minimum(b)>=0&&minimum(m)>=0&&increment>=0||error("Positivity or irreversible-memory failure")
  abs(ae)/initial_a<2e-6&&abs(be)/initial_b<2e-6&&abs(phase_mean)<2e-6||error("Stoichiometric or signed-phase conservation failure")
  push!(rows,(step=frame.step,time=frame.time,a_mass=ma,b_mass=mb,mineral_mass=mm,a_plus_half_m_absolute_error=ae,b_plus_half_m_absolute_error=be,total_mass_absolute_error=ae+be,a_invariant_relative_error=ae/initial_a,b_invariant_relative_error=be/initial_b,phase_mean_error=phase_mean,minimum_mineral_increment=increment,phase_range=extrema(phi),source_sha256=(a=asha,b=bsha,mineral=msha,phi=phisha)))
  previous=copy(m)
 end
 result=(schema="independent-brine-field-audit",shape=shape,cell_volume=volume,initial_a_inventory=initial_a,initial_b_inventory=initial_b,source_run_sha256=bytes2hex(sha256(read(joinpath(run,"run.json")))),auditor_sha256=bytes2hex(sha256(read(@__FILE__))),history=rows,interpretation="Direct sums of independently saved fields verify A+M/2 and B+M/2 separately, signed phase mean and nondecreasing mineral. No hidden mean adjustment or budget correction is applied. Floating-point residuals are reported against each inventory scale; no chemical free-energy closure is asserted.")
 open(io->JSON3.pretty(io,result),out,"w");println(JSON3.write(result));result
end
if abspath(PROGRAM_FILE)==@__FILE__;audit_brine(ARGS[1],ARGS[2]);end
