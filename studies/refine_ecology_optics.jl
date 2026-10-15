# Exact-state optical diagnostics: no interpolation and no new dynamics.
using Interior, Statistics, SHA
using Interior: JSON3
const ROOT=normpath(joinpath(@__DIR__,".."))
const BASE=joinpath(ROOT,"outputs","refinement-20261015","ecology")
const IDS=["ecology-gyroid-gardens","ecology-root-canals","ecology-finite-orchard"]
function components(mask)
 seen=falses(size(mask));sizes=Int[]
 for start in CartesianIndices(mask)
  (!mask[start] || seen[start]) && continue
  queue=[start];seen[start]=true;count=0
  while !isempty(queue)
   cell=pop!(queue);count+=1
   for d=1:3,sgn in (-1,1)
    next=CartesianIndex(ntuple(j->mod1(cell[j]+(j==d ? sgn : 0),size(mask,j)),3))
    mask[next]&&!seen[next]&&(push!(queue,next);seen[next]=true)
   end
  end
  push!(sizes,count)
 end
 (components=length(sizes),occupied_cells=sum(mask),largest_cells=maximum(sizes;init=0),largest_fraction=isempty(sizes) ? 0. : maximum(sizes)/sum(sizes))
end
function runone(slug)
 src=joinpath(ROOT,"outputs","collection-20261015",slug);out=joinpath(BASE,slug);mkpath(out)
 manifest=JSON3.read(read(joinpath(src,"state","state.json"),String),Dict{String,Any});dims=Tuple(Int.(manifest["grid"]["shape"]))
 function field(desc)
  T=desc["dtype"]=="float64" ? Float64 : Float32
  reshape(open(io->read!(io,Vector{T}(undef,prod(desc["shape"]))),joinpath(src,"state",desc["file"])),Tuple(Int.(desc["shape"])))
 end
 dep=field(manifest["fields"]["deposit"]);bs=manifest["processes"][3]["state"];ps=manifest["processes"][4]["state"]
 film=field(bs["film"]);perm=field(bs["permeability"]);signal=field(bs["signal"])
 # Phase state has a wrapped reference; use exact already exported raw channels.
 raw(name)=reshape(open(io->read!(io,Vector{Float32}(undef,prod(dims))),joinpath(src,"volume",name*".f32")),dims)
 phi=raw("phase");interface=raw("interface");bio=film.+dep
 scale=Float32(quantile(vec(bio),.99)); scaled=bio./max(scale,1f-8)
 iscale=Float32(quantile(vec(interface),.99));iface=clamp.(interface./max(iscale,1f-8),0,1)
 threshold=.2scale;comp=components(bio.>threshold)
 met=(source=src,native_shape=dims,simulation_step=manifest["step"],boundary_mass=sum(Float64,film),deposit_mass=sum(Float64,dep),material_p99=scale,threshold_for_components=threshold,periodic_face_connected_components=comp,film_interface_correlation=cor(vec(film),vec(interface)),deposit_interface_correlation=cor(vec(dep),vec(interface)),film_deposit_correlation=cor(vec(film),vec(dep)),biological_material_fraction_above_threshold=sum(bio.>threshold)/length(bio),description="Funded film and controller deposit are actual material fields. Their connected components are geometric observations, not organisms. Phase geometry was imposed at initialization. No spatial upsampling or added simulation.")
 open(io->JSON3.pretty(io,met),joinpath(out,"diagnostics.json"),"w")
 # Native unsmoothed central slices expose the coarse grid rather than disguising it.
 open(joinpath(out,"native-sections.svg"),"w") do io
  print(io,"<svg xmlns='http://www.w3.org/2000/svg' width='1100' height='330' viewBox='0 0 1100 330'><rect width='1100' height='330' fill='#10151a'/>")
  for (q,(label,data)) in enumerate((("Phase (imposed scaffold)",(phi.+1)./2),("Funded film",film./max(quantile(vec(film),.99),1f-8)),("Controller deposit",dep./max(quantile(vec(dep),.99),1f-8)),("Transport permeability",perm)))
   ox=10+(q-1)*275;print(io,"<text x='$ox' y='22' fill='white' font-size='14'>$label — native 32³</text>")
   for j=1:dims[2],i=1:dims[1]
    v=round(Int,255clamp(data[i,j,cld(dims[3],2)],0,1));print(io,"<rect x='$(ox+(i-1)*8)' y='$(38+(j-1)*8)' width='8' height='8' fill='rgb($v,$v,$v)'/>")
   end
  end
  print(io,"</svg>")
 end
 for mode in ("funded-material","scaffold-context")
  dest=joinpath(out,mode);mkpath(dest)
  density=mode=="funded-material" ? .03f0.+1.6f0.*scaled : .02f0.+1.3f0.*scaled.+.08f0.*iface
  # Nonclipped optical scaling preserves relative high-density structures.
  export_volume((density=density,activity=.05f0.*scaled,structure=scaled,boundary_material=scaled,film=film,deposit=dep,phase=phi,interface=interface,permeability=perm,signal=signal,mask=ones(Float32,dims)),joinpath(dest,"volume"))
  desc="Exact1000-step native32³ material observation: $(mode). Opacity and frosted isosurface follow actual funded film + controller deposit divided by its99th percentile; no invented organization, interpolation, or agent markers. Scaffold context additionally shows phase-interface density at low opacity."
  args=["--palette","amber-gold","--emission","0.05","--absorption","1.3","--scattering",".4","--density-scale","1","--inclusion-surface","boundary_material","--inclusion-isovalue",".2","--inclusion-transmission",".25","--surface-roughness",".28","--inclusion-tint","jade","--orbit",".65","--elevation",".4","--frame-content"]
  job=(id=slug*"-"*mode,family="ecology-refinement",title=replace(slug,"ecology-"=>"")*" / "*mode,description=desc,volume="volume",renderer="volume",args=args,size=1600,samples=96)
  open(io->JSON3.pretty(io,job),joinpath(dest,"render-job.json"),"w")
  open(io->JSON3.pretty(io,(source_state=joinpath(src,"state","state.json"),state_sha256=bytes2hex(sha256(read(joinpath(src,"state","state.json")))),study_sha256=bytes2hex(sha256(read(@__FILE__))),native_shape=dims,simulation_steps=1000,optical_only=true,biological_material_p99=scale,mode=mode,description=desc)),joinpath(dest,"provenance.json"),"w")
 end
 println(slug," ",met)
end
for id in IDS;runone(id);end
