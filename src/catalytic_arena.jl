"""Build cells in a funded habitat whose raw food requires constructed catalysts.
Resource and feedstock initial/source/sink values are physical concentrations or
rate densities. Source budgets remain total amounts, as in cellular recipes.
"""
function build_catalytic_arena(cfg)
    base=deepcopy(cfg);base["model"]="cellular"
    base["resources"]=merge(Dict{String,Any}("initial"=>[.12,.12,0.],"source_rates"=>[0.,0.,0.],
        "source_budgets"=>[0.,0.,0.],"diffusion"=>[.12,.12,.03]),get(base,"resources",Dict()))
    original=build_cellular_arena(base);w=original.world
    c=only(filter(p->p isa CellularEcology,original.processes));r=c.resources
    options=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in get(cfg,"feedstock",Dict()))
    defaults=Dict{Symbol,Any}(:names=>[:bound_a,:bound_b],:potentials=>copy(r.potentials[1:2]),
        :initial=>[1.5,1.5],:diffusion=>[0.,0.],:source_rates=>[.0025,.0025],:source_budgets=>[1e5,1e5],
        :sink_rates=>[0.,0.],:sink_budgets=>[0.,0.])
    options=merge(defaults,options);volume=c.stencil.cell_volume
    for key in (:initial,:source_rates,:sink_rates);options[key]=Float64.(options[key]).*volume;end
    options[:names]=Symbol.(options[:names]);options[:source_masks]=deepcopy(r.source_masks[1:2])
    feedstock=ResourcePool(w;options...)
    for k in 1:2;feedstock.fields[k].*=feedstock.source_masks[k];end
    arguments=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in get(cfg,"catalysis",Dict()))
    habitat=CatalyticHabitat(w,c,feedstock;arguments...)
    processes=(feedstock,original.processes...,habitat)
    if haskey(cfg,"matrix") && get(cfg["matrix"],"enabled",true)
        arguments=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in cfg["matrix"] if k!="enabled")
        matrix=SorptiveMatrix(w,habitat;arguments...)
        processes=(processes...,matrix)
    end
    Arena(w,processes;dt=get(cfg,"dt",.1),config=cfg)
end

function catalytic_arena_snapshot!(fields,agents,metrics,p,w,xs,ys,k)
    metrics["catalytic"]=Dict(string(key)=>value for (key,value) in pairs(catalytic_statistics(p)))
    metrics["feedstock"]=Dict(string(key)=>value for (key,value) in pairs(resource_statistics(p.substrates)))
    for (name,arr,units) in (("enzyme_a",p.enzymes[1],"catalyst mass / physical volume"),
            ("enzyme_b",p.enzymes[2],"catalyst mass / physical volume"),
            ("release_a",p.release[1],"precursor mass / physical volume / time"),
            ("release_b",p.release[2],"precursor mass / physical volume / time"))
        factor=startswith(name,"enzyme") ? 1/p.cells.stencil.cell_volume : 1.
        fields[name]=Dict("shape"=>[length(xs),length(ys)],"values"=>vec(arr[xs,ys,k:k]).*factor,"units"=>units)
    end
    index=Dict(id=>i for (i,id) in enumerate(p.ids))
    for agent in agents
        i=get(index,agent["id"],0);i==0 && continue
        agent["catalysis"]=Dict("investment"=>p.expression[1,i],"allocation"=>p.expression[2:3,i],
            "last_secretion"=>p.last_secretion[:,i],
            "genome_signature"=>bytes2hex(SHA.sha256(reinterpret(UInt8,vec(p.weights[:,:,i]))))[1:12])
    end
    metrics["cellular_model"]="growing bodies / inherited development / constructed catalytic habitat"
    nothing
end
