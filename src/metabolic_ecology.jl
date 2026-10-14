"""A multisubstrate metabolism composed with the existing recurrent ecology.

Traits are inherited genotype, not lifetime neural marks. Allocation columns are
simplexes: specializing in one substrate removes capacity from the others.
Efficiency increases yield but decreases catalytic throughput. Structural matter
returned at death is inert; only the ResourcePool contains chemical fuel.
"""
mutable struct MetabolicEcology
    controller::NeuralEcology
    resources::ResourcePool
    boundary::Any
    workspace::Union{Nothing,NeuralWorkspace}
    allocation::Matrix{Float32}
    assimilation::Vector{Float32}
    efficiency::Vector{Float32}
    cooperation::Vector{Float32}
    byproduct_fraction::Vector{Float32}
    trait_mutation::Float32
    uptake_rate::Float32
    boundary_rate::Float32
    chemical_heat::Float64
    consumed::Vector{Float64}
    contributed_material::Dict{Int,Float64}
    contributed_work::Dict{Int,Float64}
    recent_contributions::Dict{Int,Float64}
end
function MetabolicEcology(positions::AbstractMatrix,resources::ResourcePool;
        boundary=nothing,workspace=nothing,allocation=nothing,assimilation=.4f0,efficiency=.6f0,
        cooperation=.25f0,byproduct_fraction=.2f0,trait_mutation=.04f0,uptake_rate=.5f0,
        boundary_rate=.04f0,kwargs...)
    n=size(positions,2); k=length(resources.names)
    a=allocation===nothing ? fill(1f0/k,k,n) : Float32.(allocation)
    size(a)==(k,n) && all(x->isfinite(x)&&x>=0,a) && all(x->x>0,vec(sum(a;dims=1))) || throw(ArgumentError("invalid enzyme allocations"))
    a./=sum(a;dims=1)
    all(x->isfinite(x)&&0<=x<=1,(assimilation,efficiency,cooperation,byproduct_fraction)) && all(x->isfinite(x)&&x>=0,(trait_mutation,uptake_rate,boundary_rate)) || throw(ArgumentError("invalid metabolic traits/rates"))
    c=NeuralEcology(positions;uptake_rate=0,kwargs...)
    MetabolicEcology(c,resources,boundary,workspace,a,fill(Float32(assimilation),n),fill(Float32(efficiency),n),fill(Float32(cooperation),n),fill(Float32(byproduct_fraction),n),Float32(trait_mutation),Float32(uptake_rate),Float32(boundary_rate),0.,zeros(k),Dict{Int,Float64}(),Dict{Int,Float64}(),Dict{Int,Float64}())
end
function reorder_ecology!(p::MetabolicEcology,order)
    reorder_ecology!(p.controller,order)
    p.allocation=p.allocation[:,order]
    for name in (:assimilation,:efficiency,:cooperation,:byproduct_fraction); setfield!(p,name,getfield(p,name)[order]); end
    p
end
function validate(w::World,p::MetabolicEcology)
    validate(w,p.controller); validate(w,p.resources)
    n=length(p.controller.ids); k=length(p.resources.names)
    size(p.allocation)==(k,n) && all(x->isfinite(x)&&x>=0,p.allocation) && all(x->isapprox(x,1;atol=2e-6),vec(sum(p.allocation;dims=1))) || throw(ArgumentError("invalid metabolic allocation"))
    all(v->length(v)==n&&all(x->isfinite(x)&&0<=x<=1,v),(p.assimilation,p.efficiency,p.cooperation,p.byproduct_fraction)) || throw(ArgumentError("invalid metabolic traits"))
    p.controller.uptake_rate==0 || throw(ArgumentError("metabolic controller must not catabolize inert nutrient"))
    all(x->isfinite(x)&&x>=0,(p.trait_mutation,p.uptake_rate,p.boundary_rate,p.chemical_heat)) || throw(ArgumentError("invalid metabolic rates/ledger"))
    nothing
end
stability_limit(w,p::MetabolicEcology)=min(1f0,stability_limit(w,p.resources))

# Keep 32-sensor recurrent architecture/checkpoint unchanged. The first two local
# channels and gradients become individual species; additional species enter a
# trait-weighted resource channel. Other medium/social sensors retain meaning.
function _metabolic_inputs!(w,p,cells,inputs)
    c=p.controller; h=size(c.hidden,1); dims=w.grid.size; k=length(p.resources.names)
    for i in eachindex(c.ids)
        for channel in 1:2
            concentration(q)=p.resources.fields[channel][q] + (channel==1 ? sum((p.allocation[s,i]*p.resources.fields[s][q] for s in 3:k);init=0f0) : 0f0)
            inputs[h+channel,i]=_neural_squash(concentration(cells[i]))
            for d in 1:3
                lo=CartesianIndex(ntuple(j->j==d ? mod1(cells[i][j]-1,dims[j]) : cells[i][j],3))
                hi=CartesianIndex(ntuple(j->j==d ? mod1(cells[i][j]+1,dims[j]) : cells[i][j],3))
                inputs[h+9+(channel-1)*3+d,i]=dims[d]==1 ? 0f0 : _neural_squash((concentration(hi)-concentration(lo))/2)
            end
        end
    end
    inputs
end

function ecology_step!(w::World,p::MetabolicEcology,dt::Real)
    validate(w,p); isfinite(dt)&&0<=dt<=1 || throw(ArgumentError("ecology dt must lie in [0,1]")); dt==0 && return p
    issorted(p.controller.ids) || reorder_ecology!(p,sortperm(p.controller.ids))
    c=p.controller; n=length(c.ids); h=size(c.hidden,1)
    for i in 1:n,d in 1:3; c.positions[d,i]=w.grid.size[d]==1 ? 0f0 : mod(c.positions[d,i],Float32(w.grid.size[d])); end
    p.workspace===nothing && (p.workspace=NeuralWorkspace(w,c))
    frame=neural_frame!(p.workspace,w,c;evaluate=false)
    cells=frame.cells; inputs=_metabolic_inputs!(w,p,cells,frame.inputs); outputs=frame.outputs
    neural_controller!(outputs,c.weights,c.plastic_marks,inputs)
    requests=zeros(Float32,length(p.resources.names),n)
    for i in 1:n,s in axes(requests,1)
        stock=p.resources.fields[s][cells[i]]
        # High efficiency has an explicit throughput cost; bounded allocation
        # and saturating availability prevent arbitrary uptake capacity.
        requests[s,i]=Float32(dt)*p.uptake_rate*p.allocation[s,i]*(1f0-.7f0*p.efficiency[i])*stock/(.2f0+stock)*((outputs[h+4,i]+1f0)/2)*_neural_medium(c,:permeability,cells[i],1f0)
    end
    amounts=allocate_resources!(p.resources,cells,requests)
    credit=zeros(Float32,n)
    empty!(p.recent_contributions)
    for i in 1:n
        taken=sum(Float64,view(amounts,:,i)); potential=sum(Float64(amounts[s,i])*p.resources.potentials[s] for s in axes(amounts,1))
        exported=0.; retained=0.
        if length(p.resources.names)>=3
            for s in 1:2
                p.resources.potentials[s]>p.resources.potentials[3] || continue
                desired=Float64(amounts[s,i])*(1-p.assimilation[i])*p.byproduct_fraction[i]
                actual=resource_add!(p.resources,3,cells[i],desired)
                exported+=actual; retained+=actual*p.resources.potentials[3]
            end
        end
        # All export is lower-potential matter; remaining oxidation can only
        # deliver the input potential not retained in structural/exported matter.
        useful=max(0.,potential*(1-p.assimilation[i])-retained)*p.efficiency[i]
        c.storage[i]+=Float32(taken*p.assimilation[i]); c.waste+=taken*(1-p.assimilation[i])-exported
        credit[i]=Float32(useful)
        c.energy[i]+=Float32(useful); c.energy_input+=useful; p.chemical_heat+=potential-retained-useful
        p.consumed.+=view(amounts,:,i)
        if p.boundary!==nothing
            gate=(outputs[h+5,i]+1f0)/2
            funded=maintain_boundary!(p.boundary,Tuple(cells[i]),c.storage[i],c.energy[i],Float32(dt)*p.boundary_rate*p.cooperation[i]*gate)
            id=c.ids[i]
            p.contributed_material[id]=get(p.contributed_material,id,0.)+funded.material
            p.contributed_work[id]=get(p.contributed_work,id,0.)+funded.work
            p.recent_contributions[id]=Float64(funded.material)
            c.storage[i]-=funded.material; c.energy[i]-=funded.work; c.heat+=funded.work
        end
    end
    oldids=copy(c.ids)
    oldallocation=p.allocation; oldassimilation=p.assimilation; oldefficiency=p.efficiency
    oldcooperation=p.cooperation; oldbyproduct=p.byproduct_fraction
    ecology_step!(w,c,dt;frame=frame,reward_credit=credit)
    c.ids==oldids && return p
    living=Set(c.ids)
    for history in (p.contributed_material,p.contributed_work,p.recent_contributions)
        for id in collect(keys(history)); id in living || delete!(history,id); end
    end
    oldindices=Dict(id=>i for (i,id) in enumerate(oldids))
    k=length(p.resources.names); p.allocation=zeros(Float32,k,length(c.ids))
    p.assimilation=zeros(Float32,length(c.ids)); p.efficiency=similar(p.assimilation); p.cooperation=similar(p.assimilation); p.byproduct_fraction=similar(p.assimilation)
    for (i,id) in enumerate(c.ids)
        born=!haskey(oldindices,id); parent=born ? c.parents[i] : id; j=oldindices[parent]
        a=view(oldallocation,:,j); assim=oldassimilation[j]; eff=oldefficiency[j]
        coop=oldcooperation[j]; byproduct=oldbyproduct[j]
        if born
            rng=MersenneTwister(string(c.seed,":metabolic:",id,":",parent))
            a=max.(0f0,a.+p.trait_mutation.*randn(rng,Float32,k)); sum(a)==0 && (a.=1f0); a./=sum(a)
            assim=clamp(assim+p.trait_mutation*randn(rng,Float32),0f0,1f0)
            eff=clamp(eff+p.trait_mutation*randn(rng,Float32),0f0,1f0)
            coop=clamp(coop+p.trait_mutation*randn(rng,Float32),0f0,1f0)
            byproduct=clamp(byproduct+p.trait_mutation*randn(rng,Float32),0f0,1f0)
        end
        p.allocation[:,i].=a; p.assimilation[i]=assim; p.efficiency[i]=eff; p.cooperation[i]=coop; p.byproduct_fraction[i]=byproduct
    end
    p
end
advance!(w::World,p::MetabolicEcology,dt)=ecology_step!(w,p,dt)
process_children(p::MetabolicEcology)=(p.controller,)
ecology_mass(w,p::MetabolicEcology)=ecology_mass(w,p.controller)+resource_mass(p.resources)+(p.boundary===nothing ? 0. : boundary_material_mass(p.boundary))
function ecology_statistics(p::MetabolicEcology)
    merge(ecology_statistics(p.controller),(chemical_heat=p.chemical_heat,consumed=copy(p.consumed),mean_efficiency=isempty(p.efficiency) ? 0. : sum(p.efficiency)/length(p.efficiency),mean_cooperation=isempty(p.cooperation) ? 0. : sum(p.cooperation)/length(p.cooperation)))
end
record_process_fields(p::MetabolicEcology)=(controller=record_process_fields(p.controller),allocation=p.allocation,assimilation=p.assimilation,efficiency=p.efficiency,cooperation=p.cooperation,byproduct_fraction=p.byproduct_fraction,chemical_heat=p.chemical_heat,consumed=p.consumed)
