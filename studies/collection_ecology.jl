using Interior, Random, SHA, TOML, LinearAlgebra
using Interior: JSON3
BLAS.set_num_threads(1)
const ROOT=normpath(joinpath(@__DIR__,".."))
const SPECS=[
("ecology-gyroid-gardens","Gyroid gardens","gyroid",0.75,.045,1.3,.7,72.,.035,"glass-ice"),
("ecology-tidal-lamellae","Tidal lamellae","lamella",.15,.08,1.8,.35,55.,.025,"amber-gold"),
("ecology-pearl-niches","Pearl niches","pearls",.55,.025,1.2,.5,90.,.08,"acoustic-pearl"),
("ecology-root-canals","Root canals","canals",.85,.09,1.6,.6,65.,.035,"glass-ice"),
("ecology-finite-orchard","Finite orchard","orchard",.2,.015,1.9,.45,60.,0.,"amber-gold"),
("ecology-crossfeed-foam","Crossfeed foam","foam",.65,.065,1.5,.4,80.,.055,"acoustic-pearl")]
function runone(spec;steps=1000,n=32)
 slug,title,mode,coop,brate,uptake,assim,life,feed,palette=spec
 out=joinpath(ROOT,"outputs","collection-20261015",slug); mkpath(out)
 cfg=TOML.parsefile(joinpath(ROOT,"experiments","arena","metabolic-3d.toml"))
 cfg["seed"]=findfirst(==(spec),SPECS)*79+117; cfg["ecology"]["population"]=640; cfg["ecology"]["capacity"]=2800
 cfg["ecology"]["lifespan"]=life; cfg["ecology"]["speed"]=.65; cfg["ecology"]["secretion_rate"]=.045
 cfg["grid"]["shape"]=[n,n,n];cfg["grid"]["extent"]=[24.,24.,24.]
 cfg["metabolism"]["cooperation"]=coop;cfg["metabolism"]["boundary_rate"]=brate;cfg["metabolism"]["uptake_rate"]=uptake;cfg["metabolism"]["assimilation"]=assim
 cfg["metabolism"]["byproduct_fraction"]=mode=="foam" ? .8 : .2
 cfg["resources"]["source_rates"]=[feed,feed,0.];cfg["resources"]["spatial_niches"]=false
 cfg["medium"]["remodeling_gain"]=mode=="canals" ? 3. : .7; cfg["medium"]["membrane_pinning"]=mode=="pearls" ? 3. : 1.
 cfg["phase"]=Dict("mean"=>mode=="pearls" ? -.3 : 0.,"mobility"=>.25)
 a=build_arena(cfg); w=a.world; resources,eco,boundary,phase=a.processes; c=eco.controller
 for k=1:n,j=1:n,i=1:n
  x,y,z=ntuple(d->2pi*((i,j,k)[d]-.5)/n,3)
  f=mode=="gyroid" ? sin(x)*cos(y)+sin(y)*cos(z)+sin(z)*cos(x) : mode=="lamella" ? cos(2z+.7sin(x)+.4cos(y)) : mode=="pearls" ? cos(2x)*cos(2y)*cos(2z)-.45 : mode=="canals" ? cos(x+.6sin(z))+cos(y-.6cos(z)) : mode=="orchard" ? cos(2x)+cos(2y)+cos(2z)-.8 : sin(2x)*cos(y)+sin(2y)*cos(z)+sin(2z)*cos(x)
  phase.phi[i,j,k]=Float32(.52tanh(2f))
  niche=Float32(.5+.47tanh(2f));resources.source_masks[1][i,j,k]=niche;resources.source_masks[2][i,j,k]=1-niche
  resources.fields[1][i,j,k]=1.2niche;resources.fields[2][i,j,k]=1.2(1-niche)
 end
 Interior.phase_observe!(w,phase); Interior.medium_observe!(eco.controller.medium)
 history=Any[];elapsed=@elapsed for tick=1:steps
  arena_step!(a)
  if tick%100==0
   push!(history,merge((tick=tick,time=w.time,),ecology_statistics(eco)));println("$slug $tick/$steps $(ecology_statistics(eco))");flush(stdout)
  end
 end
 optical=optical_fields(w);interface=Array(phase.interface); iscale=maximum(interface); visible=clamp.(interface./max(iscale,1f-8),0,1)
 fields=(density=.015f0.+.38f0.*visible.^2 .+ 1.8f0.*boundary.film .+ .6f0.*w.fields.deposit,activity=.12f0.*visible.^3,structure=Array(w.fields.deposit).+boundary.film,phase=Array(phase.phi),interface=interface,film=Array(boundary.film),deposit=Array(w.fields.deposit),substrate_a=Array(resources.fields[1]),substrate_b=Array(resources.fields[2]),byproduct=Array(resources.fields[3]),mask=optical.mask)
 export_volume(fields,joinpath(out,"volume"));record_state(w,a.processes,joinpath(out,"state");step=a.steps)
 pop=(dimension=3,ids=c.ids,parents=c.parents,lineages=c.lineages,generations=c.generations,positions=[collect(c.positions[:,i]) for i in eachindex(c.ids)],traits=(allocation=eco.allocation,assimilation=eco.assimilation,efficiency=eco.efficiency,cooperation=eco.cooperation),statistics=ecology_statistics(eco))
 open(io->JSON3.pretty(io,pop),joinpath(out,"population.json"),"w")
 description="Three-dimensional metabolic population in $(mode) phase niches, with finite source budgets, inherited substrate allocation and paid boundary maintenance. Snapshot after $steps steps; no adaptation or multicellularity claim."
 cfg["collection"]=Dict("id"=>slug,"title"=>title,"mode"=>mode,"steps"=>steps,"study"=>"studies/collection_ecology.jl")
 for path in (joinpath(out,"recipe.toml"),joinpath(ROOT,"experiments","collection","ecology",slug*".toml"))
  mkpath(dirname(path));open(io->TOML.print(io,cfg),path,"w")
 end
 open(io->JSON3.pretty(io,(id=slug,family="ecology",title=title,description=description,shape=w.grid.size,steps=steps,time=w.time,metrics=ecology_statistics(eco),history=history,elapsed_seconds=elapsed,sourcehash=bytes2hex(sha256(read(@__FILE__))),sources=Interior.arena_sources(),optical_transfer="Actual 3D phase interface normalized by its maximum; density combines interface squared, film and deposit. Markers are all surviving agents at actual xyz; lineage colors are categorical.")),joinpath(out,"run.json"),"w")
 args=["--palette",palette,"--population","population.json","--population-radius","0.009","--emission","0.12","--absorption","0.45","--scattering","0.1","--density-scale","0.65","--polished-carrier","--carrier-shape","box","--ior","1.38","--orbit","0.7","--elevation","0.35"]
 open(io->JSON3.pretty(io,(id=slug,family="ecology",title=title,description=description,volume="volume",renderer="volume",args=args,size=1600,samples=96)),joinpath(out,"render-job.json"),"w")
 println("READY $slug");flush(stdout)
end
for (i,spec) in enumerate(SPECS)
 isempty(ARGS) || string(i) in split(ARGS[1],',') || continue
 runone(spec;steps=parse(Int,get(ENV,"ECOLOGY_STEPS","1000")));GC.gc()
end
selected=[spec[1] for spec in SPECS if isfile(joinpath(ROOT,"outputs","collection-20261015",spec[1],"run.json"))]
open(io->JSON3.pretty(io,selected),joinpath(ROOT,"outputs","collection-20261015","collection-ecology-selected.json"),"w")
