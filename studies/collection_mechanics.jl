using Interior, JSON3, TOML, SHA, LinearAlgebra
BLAS.set_num_threads(1)
isdefined(Interior,:BondNetwork) || Base.include(Interior,joinpath(@__DIR__,"..","src","mechanics.jl"))
const ROOT=normpath(joinpath(@__DIR__,".."))
const SPECS=[
("fracture-golden-sail","Golden sail",(20,20,3),(1.05,1.25,.13),.33,.18,5.,.20,.30,.14,1.35,"amber-gold"),
("fracture-ice-cleavage","Ice cleavage",(16,16,16),(1.0,1.0,1.0),.34,.015,11.,.14,.23,.08,1.30,"glass-ice"),
("fracture-shear-veil","Shear veil",(20,18,4),(1.1,1.0,.22),.08,.65,2.,.18,.40,.22,1.45,"glass-ice"),
("fracture-copper-lace","Copper lace",(20,16,5),(1.0,1.15,.30),.27,.32,1.,.26,.52,.27,1.65,"amber-gold"),
("fracture-parted-choir","Parted choir",(18,18,8),(1.05,1.0,.52),.43,-.24,7.,.12,.18,.18,1.50,"acoustic-pearl"),
("fracture-held-breath","Held breath",(18,18,10),(1.0,1.0,.65),.19,.27,3.,.32,.43,.24,1.35,"amber-cyan")]
function runone(s;pilot=false)
 slug,title,dims,extent,stretch,shear,damping,failure,seam,jitter,duration,palette=s
 pilot && (dims=(8,8,3);duration=.06)
 out=joinpath(ROOT,"outputs","collection-20261015",pilot ? "pilot-"*slug : slug);mkpath(out)
 w=World(Grid(24,24,24);seed=420+findfirst(==(s),SPECS))
 p=Interior.BondNetwork(w;size=dims,extent=extent,stiffness=8.,mass=.02,damping=damping,failure_strain=failure,weak_seam=seam,stretch_rate=stretch,shear_rate=shear,jitter=jitter)
 initial=Interior.mechanical_energy(p).total;dt=.001;steps=ceil(Int,duration/dt)
 elapsed=@elapsed for k in 1:steps
  step!(w,p,dt)
  k%250==0 && println(slug," step ",k,"/",steps," broken ",count(==(1),p.damage))
 end
 Interior.validate(w,p)
 Interior.export_bond_graph(p,joinpath(out,"channels.obj"))
 record_state(w,p,joinpath(out,"record_state");step=steps)
 open(io->JSON3.pretty(io,(array_order="column-major",position_shape=size(p.position),bonds_shape=size(p.bonds),position=p.position,reference=p.reference,displacement=p.position-p.reference,velocity=p.velocity,node_ids=p.node_ids,bonds=p.bonds,bond_ids=p.bond_ids,surviving_bond_ids=p.bond_ids[p.damage.==0],damage=p.damage,peak_strain=p.peak_strain,threshold=p.threshold,rest_length=p.rest_length)),joinpath(out,"graph.json"),"w")
 # The transparent constant carrier is optical staging; the OBJ is the actual surviving graph.
 z=zeros(Float32,24,24,24);export_volume((density=fill(.004f0,size(z)),activity=z,structure=z,mask=ones(Float32,size(z))),joinpath(out,"volume"))
 recipe=Dict("id"=>slug,"title"=>title,"study"=>"studies/collection_mechanics.jl","nodes_shape"=>collect(dims),"extent"=>collect(extent),"seed"=>w.seed,"stiffness"=>8.,"mass"=>.02,"stretch_rate"=>stretch,"shear_rate"=>shear,"damping"=>damping,"failure_strain"=>failure,"weak_seam"=>seam,"jitter"=>jitter,"dt"=>dt,"steps"=>steps)
 for path in (joinpath(out,"recipe.toml"),joinpath(ROOT,"experiments","collection","mechanics",slug*".toml"))
  pilot && path!=joinpath(out,"recipe.toml") && continue
  mkpath(dirname(path));open(io->TOML.print(io,recipe),path,"w")
 end
 energy=Interior.mechanical_energy(p);residual=energy.total-initial-energy.external_work
 metrics=(nodes=length(p.node_ids),bonds=length(p.bond_ids),broken=count(==(1),p.damage),surviving=count(==(0),p.damage),energy=energy,energy_balance_residual=residual,relative_work_residual=residual/max(abs(energy.external_work),eps()),max_displacement=maximum(norm(p.position[:,i]-p.reference[:,i]) for i in eachindex(p.node_ids)))
 description="Disordered central-force spring network under tensile and shear grip loading; measured tensile strain irreversibly removes bonds. Transparent carrier is optical staging. No healing or continuum calibration."
 open(io->JSON3.pretty(io,(id=slug,family="mechanics",title=title,description=description,time=w.time,params=recipe,metrics=metrics,elapsed_seconds=elapsed,sourcehash=bytes2hex(sha256(read(@__FILE__))),optical_transfer="Exact deformed surviving OBJ segments in gallery filament material; constant near-transparent density carrier, no density-derived fracture geometry.")),joinpath(out,"run.json"),"w")
 args=["--bond-graph","channels.obj","--graph-radius","0.0012","--palette",palette,"--emission","0","--absorption","0.12","--scattering","0.02","--orbit",slug=="fracture-shear-veil" ? "0.8" : "0.55","--elevation","0.36"]
 open(io->JSON3.pretty(io,(id=slug,family="mechanics",title=title,description=description,volume="volume",renderer="volume",args=args,size=1600,samples=96)),joinpath(out,"render-job.json"),"w")
 println("READY ",slug," ",metrics)
end
if "pilot" in ARGS
 runone(SPECS[1];pilot=true)
else
 for (i,s) in enumerate(SPECS)
  !isempty(ARGS) && !(string(i) in split(ARGS[1],',')) && continue
  runone(s);GC.gc()
 end
end
