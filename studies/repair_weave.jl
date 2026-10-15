# A CPU graph sketch of irreversible fracture and separately funded repair bonds.
# Original bonds are never reset. New bonds have new IDs, a stress-free rest length
# at insertion, and their own later fracture history.
using Interior, LinearAlgebra, SHA, TOML
using Interior: JSON3
function repair_weave(config,out)
 c=TOML.parsefile(config);ispath(out)&&error("Refusing existing run");mkpath(out);cp(config,joinpath(out,"recipe.toml"));cp(@__FILE__,joinpath(out,"study-source.jl"))
 source=read(joinpath(dirname(@__DIR__),"src","mechanics.jl"));write(joinpath(out,"mechanics-source.jl"),source)
 w=World(Grid(1,1,1;extent=(2,2,2));seed=911)
 dims=Tuple(c["dimensions"]);400*13prod(dims)<get(c,"memory_budget_mb",512)*1024^2 || error("Graph memory estimate exceeds study budget");p=BondNetwork(w;size=dims,extent=(1.6,1.35,.075),stiffness=5,mass=get(c,"node_mass",.018),damping=7,failure_strain=.17,weak_seam=.22,jitter=.10,stretch_rate=0,shear_rate=0)
 for q in eachindex(p.node_ids)
  x,y,z=p.position[:,q];p.position[3,q]+=get(c,"fold_amplitude",.15)*sin(get(c,"fold_frequency",2.8)*x+1.5y)+get(c,"ripple_amplitude",.035)*cos(7y)
 end
 p.reference.=p.position;initial_reference=copy(p.reference)
 for e in eachindex(p.bond_ids);i,j=p.bonds[:,e];p.rest_length[e]=norm(p.position[:,j]-p.position[:,i]);end
 original=length(p.bond_ids);birth=zeros(length(p.bond_ids));broken_at=fill(NaN,length(p.bond_ids));parent=zeros(Int,length(p.bond_ids));generation=zeros(Int,length(p.bond_ids));repaired=falses(length(p.bond_ids));stock=fill(Float64(c["repair_stock"]),length(p.node_ids));funded=0.;initialstock=sum(stock);initialenergy=Interior.mechanical_energy(p).total
 dt=c["dt"];frames=Any[];history=Any[];dt<stability_limit(w,p)||error("Initial timestep too large")
 function snapshot(step)
  directory=joinpath(out,"frames",lpad(step,5,'0'));mkpath(directory)
  energy=Interior.mechanical_energy(p);strain=[norm(p.position[:,p.bonds[2,e]]-p.position[:,p.bonds[1,e]])/p.rest_length[e]-1 for e in eachindex(p.bond_ids)]
  data=(time=p.elapsed,step=step,dimensions=dims,array_order="F",positions=permutedims(p.position),velocities=permutedims(p.velocity),reference=permutedims(p.reference),node_ids=p.node_ids,bonds=permutedims(p.bonds),bond_ids=p.bond_ids,rest_length=p.rest_length,stiffness=p.stiffness,damage=p.damage,strain=strain,peak_strain=p.peak_strain,birth_time=birth,parent_bond=parent,generation=generation,repair_stock=stock,original_bonds=original)
  open(io->JSON3.pretty(io,data),joinpath(directory,"geometry.json"),"w")
  residual=energy.total-initialenergy-energy.external_work
  metric=(step=step,time=p.elapsed,nodes=length(p.node_ids),original_bonds=original,new_bonds=length(p.bond_ids)-original,live_bonds=count(==(0),p.damage),broken_bonds=count(==(1),p.damage),broken_repairs=count(e->e>original&&p.damage[e]==1,eachindex(p.bond_ids)),stock=sum(stock),funded_stock=funded,stock_error=sum(stock)+funded-initialstock,energy=energy,energy_residual=residual,finite=all(isfinite,p.position)&&all(isfinite,p.velocity))
  metric.finite&&minimum(stock)>=0&&abs(metric.stock_error)<1e-8 || error("Repair invariant failed")
  push!(frames,(step=step,time=p.elapsed,directory=relpath(directory,out)));push!(history,metric);println(JSON3.write(metric));flush(stdout)
 end
 snapshot(0)
 elapsed=@elapsed for step=1:c["steps"]
  t=p.elapsed+dt;phase=2pi*t/c["cycle_time"]
  extension=.30*(1-cos(phase))/2+.063t/c["cycle_time"];shear=.16sin(phase)
  for q in eachindex(p.grip)
   p.grip[q]==0&&continue
   angle=p.grip[q]*get(c,"grip_twist",0.)*sin(phase)
   y,z=initial_reference[2,q],initial_reference[3,q]
   p.reference[2,q]=cos(angle)*y-sin(angle)*z;p.reference[3,q]=sin(angle)*y+cos(angle)*z
  end
  p.stretch_rate=extension/t;p.shear_rate=shear/t
  advance!(w,p,dt);w.time=p.elapsed
  for e in eachindex(p.bond_ids)
   p.damage[e]==1&&isnan(broken_at[e])&&(broken_at[e]=p.elapsed)
  end
  if step%c["repair_every"]==0
   edges=Tuple{Int,Int}[];parents=Int[];lengths=Float64[];costs=Float64[]
   for e in eachindex(p.bond_ids)
    p.damage[e]==1&&!repaired[e]&&generation[e]<2&&p.elapsed-broken_at[e]>=c["repair_delay"] || continue
    i,j=p.bonds[:,e];l=norm(p.position[:,j]-p.position[:,i]);cost=.18l
    1e-8<l<=1.8p.rest_length[e]&&min(stock[i],stock[j])>=cost/2 || continue
    stock[i]-=cost/2;stock[j]-=cost/2;funded+=cost;repaired[e]=true
    push!(edges,(i,j));push!(parents,e);push!(lengths,l);push!(costs,cost)
   end
   if !isempty(edges)
    firstid=maximum(p.bond_ids)+1;newcount=length(edges)
    p.bonds=hcat(p.bonds,hcat(([i,j] for (i,j) in edges)...));append!(p.bond_ids,firstid:firstid+newcount-1)
    append!(p.rest_length,lengths);append!(p.stiffness,[.7p.stiffness[e] for e in parents]);append!(p.threshold,fill(.22,newcount));append!(p.damage,zeros(newcount));append!(p.peak_strain,zeros(newcount))
    all(lengths .> 0) && all(generation[parents] .< 2) || error("Invalid repair ancestry")
    # Current separation is the new rest length: insertion has zero elastic energy.
    all(q->p.damage[parents[q]]==1 && p.rest_length[end-newcount+q]==norm(p.position[:,edges[q][2]]-p.position[:,edges[q][1]]),1:newcount) || error("Repair must preserve broken parent and begin stress free")
    append!(birth,fill(p.elapsed,newcount));append!(broken_at,fill(NaN,newcount));append!(parent,[p.bond_ids[e] for e in parents]);append!(generation,[generation[e]+1 for e in parents]);append!(repaired,falses(newcount))
   end
  end
  (step%c["capture_every"]==0||step==c["steps"])&&snapshot(step)
 end
 report=(id=c["id"],title=c["title"],family="interfaces",description=c["description"],steps=c["steps"],time=p.elapsed,dt=dt,elapsed_seconds=elapsed,params=c,frames=frames,history=history,source_sha256=Dict("src/mechanics.jl"=>bytes2hex(sha256(source)),"studies/repair_weave.jl"=>bytes2hex(sha256(read(@__FILE__)))),simulation="Existing central-force BondNetwork with cyclic externally driven stretching, shearing and opposing grip rotation, irreversible tensile failure, and study-local delayed repair. Finite node repair stock funds new lower-stiffness bonds at their current separation; old damage and bond IDs are preserved. Repair stock is a normalized construction budget, not a calibrated mass model; node inertia stays fixed.",initialization="Deliberately curved, disordered layered lattice and weakened notch. All later positions and bond failures follow mechanical dynamics. No simulated optical feedback.")
 open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w");println("READY ",out)
end
if abspath(PROGRAM_FILE)==@__FILE__;repair_weave(ARGS[1],ARGS[2]);end
