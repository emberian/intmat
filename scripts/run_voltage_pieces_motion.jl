include(joinpath(@__DIR__,"../studies/voltage_between_pieces.jl"))
function capture_run(out)
 ispath(out)&&error("fresh capture required");mkpath(out);frames=Any[]
 function observe(s)
  s.step>=319||return
  p=s.process;ordinal=length(frames);folder=joinpath(out,"states");mkpath(folder)
  ordinal==0&&open(io->JSON3.write(io,(dims=(48,36,20),reference=p.reference,bonds=p.bonds)),joinpath(out,"base.json"),"w")
  file=joinpath(folder,lpad(ordinal,5,'0')*".json")
  open(io->JSON3.write(io,(step=s.step,time=s.time,position=p.position,damage=p.damage,strain=p.peak_strain,potential=s.potential,dose=s.dose,temperature=s.temperature)),file,"w")
  energy=Interior.mechanical_energy(p)
  push!(frames,(index=ordinal,step=s.step,time=s.time,file=relpath(file,out),sha256=bytes2hex(sha256(read(file))),broken=count(==(1),p.damage),maxdose=maximum(s.dose),solver=s.solver,mechanical_energy=energy,mechanical_relative_work_residual=(energy.total-energy.external_work)/max(abs(energy.external_work),1e-30)))
  println("CAPTURE ",ordinal," step ",s.step," broken ",last(frames).broken);flush(stdout)
 end
 produce(joinpath(out,"simulation");dims=(48,36,20),steps=750,scaled=true,snapshot=observe,electrical_cadence=32)
 length(frames)==40||error("Expected40actual states, got$(length(frames))")
 meta=(title="The field gives way",actual_stored_states=true,restart_checkpoints=false,temporal_interpolation=false,native_shape=(48,36,20),native_nodes=34560,native_bonds=419228,frames=frames,base_sha256=bytes2hex(sha256(read(joinpath(out,"base.json")))),global_dose_optical_max=maximum(f.maxdose for f in frames),source_sha256=Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,"..",f)))) for f in ("studies/voltage_between_pieces.jl","scripts/run_voltage_pieces_motion.jl","src/mechanics.jl")),qualification="Freshboundedtrajectory. Electrical cadence32native steps vs52in publishedstill, a temporaldiscretization refinement; samephysicalbody/forcing/rates and .525nominalduration. Eachcapture has consistent currentdamage, positions and freshlysolved potential. Observations, not fullrestartcheckpoints.")
 open(io->JSON3.pretty(io,meta),joinpath(out,"capture.json"),"w")
 cp(@__FILE__,joinpath(out,"source-runner.jl"))
end
capture_run(ARGS[1])
