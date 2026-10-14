"""Finite CPU chemical fields, independent of World's structural nutrient/deposit.
Potentials are energy per unit mass; downhill conversion cannot recharge substrate.
"""
mutable struct ResourcePool
    names::Vector{Symbol}
    potentials::Vector{Float32}
    fields::Vector{Array{Float32,3}}
    diffusion::Vector{Float32}
    source_rates::Vector{Float32}
    sink_rates::Vector{Float32}
    source_masks::Vector{Array{Float32,3}}
    sink_masks::Vector{Array{Float32,3}}
    source_remaining::Vector{Float64}
    sink_remaining::Vector{Float64}
    supplied::Vector{Float64}
    removed::Vector{Float64}
    withdrawn::Vector{Float64}
    returned::Vector{Float64}
    released_energy::Float64
    medium::Any
    scratch::Array{Float32,3}
    permeability::Array{Float32,3}
end
function ResourcePool(w::World; names=(:substrate_a,:substrate_b,:byproduct),
        potentials=(4,3,1), initial=(1,1,0), diffusion=(.02,.01,.01),
        source_rates=ntuple(_->0,length(names)), sink_rates=ntuple(_->0,length(names)),
        source_budgets=ntuple(_->0,length(names)), sink_budgets=ntuple(_->0,length(names)),
        source_masks=nothing, sink_masks=nothing, medium=nothing)
    w.backend isa CPU || throw(ArgumentError("ResourcePool currently requires CPU"))
    n=length(names)
    n>=2 && length(unique(names))==n || throw(ArgumentError("need distinct resource names"))
    all(x->length(x)==n,(potentials,initial,diffusion,source_rates,sink_rates,source_budgets,sink_budgets)) || throw(DimensionMismatch("species parameters"))
    all(x->isfinite(x)&&x>=0,Iterators.flatten((potentials,initial,diffusion,source_rates,sink_rates,source_budgets,sink_budgets))) || throw(ArgumentError("resource parameters must be finite nonnegative"))
    masks(x)=x===nothing ? [ones(Float32,w.grid.size) for _ in 1:n] : [Float32.(a) for a in x]
    sm=masks(source_masks); km=masks(sink_masks)
    length(sm)==n && length(km)==n || throw(DimensionMismatch("resource masks"))
    p=ResourcePool(collect(Symbol,names),Float32.(collect(potentials)),
        [fill(Float32(x),w.grid.size) for x in initial],Float32.(collect(diffusion)),
        Float32.(collect(source_rates)),Float32.(collect(sink_rates)),sm,km,Float64.(collect(source_budgets)),Float64.(collect(sink_budgets)),
        zeros(n),zeros(n),zeros(n),zeros(n),0.0,medium,zeros(Float32,w.grid.size),ones(Float32,w.grid.size))
    validate(w,p); p
end
function validate(w::World,p::ResourcePool)
    w.backend isa CPU || throw(ArgumentError("ResourcePool currently requires CPU"))
    all(a->size(a)==w.grid.size,vcat(p.fields,p.source_masks,p.sink_masks)) || throw(DimensionMismatch("resource grid"))
    all(a->all(x->isfinite(x)&&x>=0,a),vcat(p.fields,p.source_masks,p.sink_masks)) || throw(ArgumentError("resources must be finite nonnegative"))
    all(x->isfinite(x)&&x>=0,Iterators.flatten((p.potentials,p.diffusion,p.source_rates,p.sink_rates,p.source_remaining,p.sink_remaining))) || throw(ArgumentError("invalid resource parameters"))
    if p.medium!==nothing
        size(p.medium.permeability)==w.grid.size || throw(DimensionMismatch("resource permeability"))
        all(x->isfinite(x)&&0<=x<=1,p.medium.permeability) || throw(ArgumentError("invalid permeability"))
    end
    nothing
end
resource_mass(p::ResourcePool)=sum(a->sum(Float64,a),p.fields)
resource_potential(p::ResourcePool)=sum(i->Float64(p.potentials[i])*sum(Float64,p.fields[i]),eachindex(p.fields))
resource_index(p::ResourcePool,s::Integer)=Int(s)
resource_index(p::ResourcePool,s::Symbol)=something(findfirst(==(s),p.names),0)>0 ? findfirst(==(s),p.names) : throw(ArgumentError("unknown resource $s"))
function resource_add!(p::ResourcePool,s,cell,amount::Real; external=false)
    isfinite(amount)&&amount>=0 || throw(ArgumentError("addition must be finite nonnegative"))
    i=resource_index(p,s); a=p.fields[i]; before=a[cell]; after=Float32(Float64(before)+amount)
    isfinite(after) || throw(ArgumentError("resource overflow"))
    Float64(after)-Float64(before)>amount && (after=prevfloat(after))
    a[cell]=after; actual=Float64(after)-Float64(before)
    (external ? p.supplied : p.returned)[i]+=actual
    actual
end
"""Allocate [species, consumer] requests proportionally within each shared cell.
Returns actual Float32 field withdrawals as Float64 values, so consumer ledgers
can account for representable field deltas rather than ideal requested amounts.
"""
function allocate_resources!(p::ResourcePool,cells,requests::AbstractMatrix)
    size(requests)==(length(p.fields),length(cells)) || throw(DimensionMismatch("resource requests"))
    all(x->isfinite(x)&&x>=0,requests) || throw(ArgumentError("requests must be finite nonnegative"))
    grants=zeros(Float64,size(requests)); groups=Dict{eltype(cells),Vector{Int}}()
    for (j,c) in enumerate(cells); push!(get!(groups,c,Int[]),j); end
    for (c,consumers) in groups, s in eachindex(p.fields)
        total=sum(j->Float64(requests[s,j]),consumers); total==0 && continue
        a=p.fields[s]; before=a[c]
        after=Float32(max(0.0,Float64(before)-min(Float64(before),total)))
        actual=Float64(before)-Float64(after); a[c]=after; p.withdrawn[s]+=actual
        for j in consumers; grants[s,j]=actual*(Float64(requests[s,j])/total); end
    end
    grants
end
"""Convert equal mass downhill locally; energy is bounded by actual lost potential.
The efficiency fraction is usable work; the remainder is dissipated heat.
"""
function resource_convert!(p::ResourcePool,from,to,cell,amount::Real; efficiency=1)
    f=resource_index(p,from); t=resource_index(p,to)
    f!=t && p.potentials[f]>=p.potentials[t] || throw(ArgumentError("conversion must go downhill between distinct species"))
    isfinite(amount)&&amount>=0 && isfinite(efficiency)&&0<=efficiency<=1 || throw(ArgumentError("invalid conversion"))
    a=p.fields[f]; b=p.fields[t]; olda=a[cell]; oldb=b[cell]
    take=min(Float64(olda),Float64(amount)); newa=Float32(Float64(olda)-take)
    actual=Float64(olda)-Float64(newa); newb=Float32(Float64(oldb)+actual)
    isfinite(newb) || throw(ArgumentError("resource overflow"))
    release=Float64(p.potentials[f])*actual-Float64(p.potentials[t])*(Float64(newb)-Float64(oldb))
    release>=0 || return (mass=0.0,energy=0.0,heat=0.0)
    a[cell]=newa; b[cell]=newb; p.released_energy+=release
    (mass=actual,energy=release*efficiency,heat=release*(1-efficiency))
end
function advance!(w::World,p::ResourcePool,dt::Real)
    check_dt(w,p,dt)
    perm=p.medium===nothing ? p.permeability : p.medium.permeability
    for s in eachindex(p.fields)
        a=p.fields[s]
        # Rates are per cell per time; budgets cap total external mass.
        sm=p.source_masks[s]; km=p.sink_masks[s]
        st=sum(Float64,sm); kt=sum(Float64,km)
        source_scale=st==0 ? 0.0 : min(p.source_remaining[s],Float64(dt)*p.source_rates[s]*st)/st
        sink_scale=kt==0 ? 0.0 : min(p.sink_remaining[s],Float64(dt)*p.sink_rates[s]*kt)/kt
        for c in eachindex(a)
            source=min(p.source_remaining[s],source_scale*sm[c]); sink=min(p.sink_remaining[s],sink_scale*km[c])
            before=a[c]
            # Round source toward the budget, never exceeding its finite amount.
            next=Float32(Float64(before)+source)
            Float64(next)-before>source && (next=prevfloat(next))
            supplied=Float64(next)-before; a[c]=next
            p.supplied[s]+=supplied; p.source_remaining[s]-=supplied
            after=Float32(max(0.0,Float64(next)-min(Float64(next),sink)))
            Float64(next)-after>sink && (after=nextfloat(after))
            removed=Float64(next)-after; a[c]=after
            p.removed[s]+=removed; p.sink_remaining[s]-=removed
        end
        D=p.diffusion[s]; D==0 && continue
        limit=diffusion_limit(w,D); steps=max(1,ceil(Int,Float64(dt)/limit))
        for _ in 1:steps
            medium_diffuse!(p.scratch,a,perm,w.grid,D,Float32(dt/steps))
        end
    end
    w
end
resource_statistics(p::ResourcePool)=(names=copy(p.names),mass=[sum(Float64,a) for a in p.fields],potential=resource_potential(p),supplied=copy(p.supplied),removed=copy(p.removed),withdrawn=copy(p.withdrawn),returned=copy(p.returned),released_energy=p.released_energy,supplied_potential=sum(p.supplied .* p.potentials),removed_potential=sum(p.removed .* p.potentials),withdrawn_potential=sum(p.withdrawn .* p.potentials),returned_potential=sum(p.returned .* p.potentials))
function record_process_fields(p::ResourcePool)
    metadata=(names=Tuple(string.(p.names)),potentials=copy(p.potentials),
        diffusion=copy(p.diffusion),source_rates=copy(p.source_rates),sink_rates=copy(p.sink_rates),
        source_remaining=copy(p.source_remaining),sink_remaining=copy(p.sink_remaining),
        supplied=copy(p.supplied),removed=copy(p.removed),withdrawn=copy(p.withdrawn),
        returned=copy(p.returned),released_energy=p.released_energy,
        statistics=merge(resource_statistics(p),(names=Tuple(string.(p.names)),)))
    keys=Tuple(Symbol(prefix,"_",i) for prefix in ("field","source_mask","sink_mask") for i in eachindex(p.fields))
    values=Tuple(copy(a) for collection in (p.fields,p.source_masks,p.sink_masks) for a in collection)
    merge(metadata,NamedTuple{keys}(values))
end
