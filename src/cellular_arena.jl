"""Cellular arena recipes use resource concentrations and physical rate densities.
Only initial/source/sink amounts are multiplied by voxel volume; source budgets
already name total material. Other arena models retain their historical units.
"""
function build_cellular_arena(cfg)
    grid=get(cfg,"grid",Dict()); dims=Tuple(Int.(get(grid,"shape",[48,48,1])))
    extent=Tuple(get(grid,"extent",[24.,24.,1.])); seed=Int(get(cfg,"seed",71))
    w=World(Grid(dims...;extent);seed); fill!(w.fields.nutrient,0)
    phase=nothing; medium=nothing
    pcfg=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in get(cfg,"phase",Dict()))
    enabled=Bool(pop!(pcfg,:enabled,false))
    if enabled
        phase=ImplicitPhaseSeparation(w;merge(Dict{Symbol,Any}(:mineral_pinning=>0,:nutrient_diffusion=>0,:interface_activity=>0,:activity_relaxation=>0),pcfg)...)
        mcfg=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in get(cfg,"medium",Dict()))
        medium=ActiveMedium(w,phase;mcfg...)
    end
    rcfg=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in get(cfg,"resources",Dict()))
    spatial=Bool(pop!(rcfg,:spatial_niches,true))
    # Smooth complementary supply is an imposed boundary condition; a second
    # channel is still available in each region, allowing solitary strategies.
    contrast=Float64(pop!(rcfg,:niche_contrast,.9))
    0<=contrast<=1 || throw(ArgumentError("niche contrast must be in [0,1]"))
    masks=[ones(Float32,dims) for _ in 1:3]
    if spatial
        for I in CartesianIndices(masks[1])
            value=Float32(1+contrast*cos(2pi*(I[1]-.5)/dims[1]))
            masks[1][I]=value; masks[2][I]=2-value
        end
    end
    volume=prod(Float64(w.grid.spacing[d]) for d in 1:(dims[3]==1 ? 2 : 3))
    defaults=Dict{Symbol,Any}(:names=>[:precursor_a,:precursor_b,:waste],:potentials=>[4.,4.,0.],
        :initial=>[.8,.8,0.],:diffusion=>[.06,.06,.015],:source_rates=>[.02,.02,0.],
        :source_budgets=>[1e6,1e6,0.],:sink_rates=>[0.,0.,0.],:sink_budgets=>[0.,0.,0.])
    options=merge(defaults,rcfg)
    for key in (:initial,:source_rates,:sink_rates); options[key]=Float64.(options[key]).*volume; end
    options[:names]=Symbol.(options[:names])
    r=ResourcePool(w;source_masks=masks,options...)
    spatial && foreach(i->r.fields[i].*=masks[i],1:2)
    ccfg=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in get(cfg,"cellular",Dict()))
    n=Int(pop!(ccfg,:population,64)); rng=MersenneTwister(seed+37)
    positions=rand(rng,3,n).*reshape(Float64.(collect(extent)),3,1)
    dims[3]==1 && (positions[3,:].=0)
    p=CellularEcology(w,positions,r;medium,seed,ccfg...)
    cellular_porosity!(w,p)
    processes=phase===nothing ? (r,p) : (r,p,medium,phase)
    Arena(w,processes;dt=get(cfg,"dt",.1),config=cfg)
end

function cellular_arena_snapshot!(fields,agents,metrics,p,w,xs,ys,k,max_agents)
    for (name,value) in pairs(cellular_statistics(p)); metrics[string(name)]=value; end
    n=length(p.ids); radii=cellular_radii(w,p)
    for i in 1:min(n,max_agents-length(agents))
        xyz=collect(view(p.positions,:,i)); actions=collect(view(p.last_actions,:,i))
        enzyme=exp.(2 .* actions[4:7]); enzyme./=sum(enzyme)
        push!(agents,Dict{String,Any}("id"=>p.ids[i],"parent"=>p.parents[i],"lineage"=>p.lineages[i],
            "generation"=>p.generations[i],"cell_position"=>xyz./collect(w.grid.spacing),
            "position"=>xyz.-collect(w.grid.extent)./2,"radius"=>radii[i],"biomass"=>p.biomass[i],
            "energy"=>p.energy[i],"storage"=>sum(view(p.metabolites,:,i)),"metabolites"=>collect(view(p.metabolites,:,i)),
            "volume"=>(p.biomass[i]+sum(view(p.metabolites,:,i)))/p.parameters.density,
            "energy_capacity"=>p.parameters.energy_capacity*p.biomass[i],"storage_capacity"=>p.parameters.store_capacity*p.biomass[i],
            "age"=>p.age[i],"plastic_norm"=>norm(view(p.plastic_marks,:,:,i)),
            "hidden"=>collect(view(p.hidden,:,i)),"adaptation"=>collect(view(p.adaptation,:,i)),"actions"=>actions,
            "genome_signature"=>bytes2hex(SHA.sha256(reinterpret(UInt8,vec(p.weights[:,:,i]))))[1:12],
            "metabolism"=>Dict("allocation"=>enzyme,"adhesion"=>(actions[8]+1)/2,"junction"=>(actions[9]+1)/2,"reserve"=>(actions[11]+1)/2)))
    end
    # Stable-ID edges are live graph state, distinct from proximity clusters.
    living=Set(p.ids)
    metrics["cell_bonds"]=[collect(edge) for edge in p.bonds.endpoints if all(id->id in living,edge)]
    metrics["cellular_model"]="complementary precursors / growing bodies / developmental bonds"
    biomass=zeros(Float32,w.grid.size); s=p.stencil
    for i in 1:n, e in s.offsets[i]:s.offsets[i+1]-1
        biomass[s.sites[e]]+=Float32(p.biomass[i]*s.weights[e]/s.cell_volume)
    end
    fields["biomass"]=Dict("shape"=>[length(xs),length(ys)],"values"=>vec(biomass[xs,ys,k:k]),"units"=>"structural mass / physical volume")
    fields["permeability"]=Dict("shape"=>[length(xs),length(ys)],"values"=>vec(p.resources.permeability[xs,ys,k:k]),"units"=>"relative transport")
    p.mechanical===nothing || (metrics["mechanics"]=p.mechanical)
    nothing
end
