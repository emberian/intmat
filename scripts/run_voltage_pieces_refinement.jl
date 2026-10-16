# Fresh high-resolution realization with compact actual observations, no law duplicate.
include(joinpath(@__DIR__,"../studies/voltage_between_pieces.jl"))
using TOML
function refine(recipe,out)
 c=TOML.parsefile(recipe);ispath(out)&&error("fresh output required");mkpath(out);cp(recipe,joinpath(out,"recipe.toml"));dims=Tuple(Int.(c["shape"]));frames=Any[]
 files=("studies/voltage_between_pieces.jl","scripts/run_voltage_pieces_refinement.jl","src/mechanics.jl")
 pins=Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,"..",f)))) for f in files)
 for f in files;cp(joinpath(@__DIR__,"..",f),joinpath(out,"source-"*basename(f)));end
 header=(title=c["title"],actual_stored_states=true,restart_checkpoints=false,temporal_interpolation=false,native_shape=dims,native_nodes=c["nodes"],native_bonds=c["bonds"],source_sha256=pins,recipe_sha256=bytes2hex(sha256(read(recipe))),qualification=c["interpretation"])
 function observe(s)
  s.time>=c["capture_start_time"]||return
  p=s.process;ordinal=length(frames);folder=joinpath(out,"states");mkpath(folder)
  ordinal==0&&open(io->JSON3.write(io,(dims=dims,reference=p.reference,bonds=p.bonds)),joinpath(out,"base.json"),"w")
  file=joinpath(folder,lpad(ordinal,5,'0')*".json")
  open(io->JSON3.write(io,(step=s.step,time=s.time,position=p.position,damage=p.damage,potential=s.potential,dose=s.dose,temperature=s.temperature)),file,"w")
  energy=Interior.mechanical_energy(p)
  push!(frames,(index=ordinal,step=s.step,time=s.time,file=relpath(file,out),sha256=bytes2hex(sha256(read(file))),broken=count(==(1),p.damage),maxdose=maximum(s.dose),solver=s.solver,mechanical_energy=energy,mechanical_relative_work_residual=(energy.total-energy.external_work)/max(abs(energy.external_work),1e-30)))
  open(io->JSON3.pretty(io,(;header...,status="running",frames=frames)),joinpath(out,"capture-progress.json"),"w")
  println("CAPTURE ",ordinal," step ",s.step," t ",s.time," broken ",last(frames).broken);flush(stdout)
 end
 produce(joinpath(out,"simulation");dims=dims,steps=c["reference_steps"],scaled=true,snapshot=observe,electrical_cadence=c["electrical_cadence_native_steps"])
 length(frames)==c["expected_captures"]||error("capture count $(length(frames))")
 meta=(;header...,status="complete",frames=frames,base_sha256=bytes2hex(sha256(read(joinpath(out,"base.json")))),global_dose_optical_max=maximum(f.maxdose for f in frames))
 open(io->JSON3.pretty(io,meta),joinpath(out,"capture.json"),"w")
end
refine(ARGS[1],ARGS[2])
