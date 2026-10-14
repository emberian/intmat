"""Finite-material ecology with heritable recurrent controllers and a separate energy ledger.
Positions are periodic cell coordinates. Energy is an abstract work budget, not material mass.
"""
mutable struct Ecology
    ids::Vector{Int}
    parents::Vector{Int}
    lineages::Vector{Int}
    generations::Vector{Int}
    positions::Matrix{Float32}
    hidden::Matrix{Float32}
    storage::Vector{Float32}
    energy::Vector{Float32}
    weights::Array{Float32,3}
    age::Vector{Float32}
    seed::Int
    events::Int
    next_id::Int
    capacity::Int
    births::Int
    deaths::Int
    blocked_births::Int
    heat::Float64
    energy_input::Float64
    waste::Float64
    uptake_rate::Float32
    basal_cost::Float32
    work_rate::Float32
    secretion_rate::Float32
    birth_material::Float32
    birth_energy::Float32
    mutation::Float32
    lifespan::Float32
    speed::Float32
    response_time::Float32
end
function Ecology(positions::AbstractMatrix; seed=17, capacity=2000, hidden=4,
                 storage=0.4f0, energy=1f0, uptake_rate=0.5f0, basal_cost=0.025f0,
                 work_rate=0.05f0, secretion_rate=0.025f0, birth_material=0.7f0,
                 birth_energy=0.7f0, mutation=0.08f0, lifespan=80f0, speed=0.3f0, response_time=1f0)
    n=size(positions,2)
    size(positions,1)==3 || throw(DimensionMismatch("positions require three rows"))
    0<=n<=capacity && hidden>=2 || throw(ArgumentError("invalid capacity or hidden size"))
    pars=Float32.((storage,energy,uptake_rate,basal_cost,work_rate,secretion_rate,birth_material,birth_energy,mutation,lifespan,speed,response_time))
    all(x->isfinite(x)&&x>=0,pars) && birth_material>0 && birth_energy>0 && lifespan>0 && response_time>0 && all(isfinite,positions) || throw(ArgumentError("invalid ecology parameters"))
    rng=MersenneTwister(seed)
    # hidden inputs: previous hidden + nutrient/deposit/storage/energy/local density + bias.
    weights=0.4f0.*randn(rng,Float32,hidden+5,hidden+6,n)
    # Splitting is selectable but not artificially rare in founders.
    weights[hidden+5,end,:].=1f0
    Ecology(collect(1:n),zeros(Int,n),collect(1:n),zeros(Int,n),Float32.(positions),zeros(Float32,hidden,n),fill(Float32(storage),n),fill(Float32(energy),n),weights,zeros(Float32,n),Int(seed),0,n+1,Int(capacity),0,0,0,0.0,Float64(n*Float32(energy)),0.0,Float32(uptake_rate),Float32(basal_cost),Float32(work_rate),Float32(secretion_rate),Float32(birth_material),Float32(birth_energy),Float32(mutation),Float32(lifespan),Float32(speed),Float32(response_time))
end
function reorder_ecology!(p::Ecology, order)
    sort(collect(order))==collect(eachindex(p.ids)) || throw(ArgumentError("invalid permutation"))
    for name in (:ids,:parents,:lineages,:generations,:storage,:energy,:age)
        setfield!(p,name,getfield(p,name)[order])
    end
    p.positions=p.positions[:,order]; p.hidden=p.hidden[:,order]; p.weights=p.weights[:,:,order]; p
end
ecology_mass(w,p::Ecology)=Float64(sum(Float64,w.fields.nutrient))+sum(Float64,w.fields.deposit)+sum(Float64,p.storage;init=0.0)+p.waste
ecology_statistics(p::Ecology)=(population=length(p.ids),births=p.births,deaths=p.deaths,blocked_births=p.blocked_births,max_generation=isempty(p.ids) ? 0 : maximum(p.generations),surviving_lineages=length(unique(p.lineages)),heat=p.heat,energy_input=p.energy_input,stored_energy=sum(Float64,p.energy;init=0.0),energy_error=sum(Float64,p.energy;init=0.0)+p.heat-p.energy_input,stored_material=sum(Float64,p.storage;init=0.0),waste=p.waste)
record_process_fields(p::Ecology)=NamedTuple{fieldnames(Ecology)}(Tuple(getfield(p,n) for n in fieldnames(Ecology)))
dimsafe(x,n)=n==1 ? 1 : floor(Int,mod(x,Float32(n)))+1
function validate(w::World,p::Ecology)
    n=length(p.ids); h=size(p.hidden,1)
    rates=(p.uptake_rate,p.basal_cost,p.work_rate,p.secretion_rate,p.birth_material,p.birth_energy,p.mutation,p.lifespan,p.speed,p.response_time)
    all(x->isfinite(x)&&x>=0,rates) && p.birth_material>0 && p.birth_energy>0 && p.lifespan>0 && p.response_time>0 && p.capacity>=0 && h>=2 && all(x->isfinite(x)&&x>=0,(p.heat,p.energy_input,p.waste)) || throw(ArgumentError("invalid ecology parameters or ledgers"))
    size(p.positions)==(3,n) && size(p.hidden)==(h,n) && size(p.weights)==(h+5,h+6,n) && all(v->length(v)==n,(p.parents,p.lineages,p.generations,p.storage,p.energy,p.age)) || throw(DimensionMismatch("ecology state"))
    n<=p.capacity && length(unique(p.ids))==n && all(<(p.next_id),p.ids) || throw(ArgumentError("ecology IDs or capacity"))
    all(a->all(isfinite,a),(p.positions,p.hidden,p.weights)) && all(a->all(x->isfinite(x)&&x>=0,a),(p.storage,p.energy,p.age)) || throw(ArgumentError("invalid ecology state"))
    all(a->a isa Array{Float32,3} && size(a)==w.grid.size,(w.fields.nutrient,w.fields.deposit,w.fields.vx,w.fields.vy,w.fields.vz)) || throw(ArgumentError("ecology requires matching CPU grids"))
    # Inspect only sampled exchange sites; do not rescan every shared grid each step.
    for i in 1:n
        c=CartesianIndex(ntuple(d->dimsafe(p.positions[d,i],w.grid.size[d]),3))
        all(x->isfinite(x)&&x>=0,(w.fields.nutrient[c],w.fields.deposit[c])) && all(isfinite,(w.fields.vx[c],w.fields.vy[c],w.fields.vz[c])) || throw(ArgumentError("invalid sampled ecology grids"))
    end
    nothing
end
stability_limit(w,p::Ecology)=1f0
function ecology_step!(w::World,p::Ecology,dt::Real)
    validate(w,p); isfinite(dt)&&0<=dt<=1 || throw(ArgumentError("ecology dt must lie in [0,1]")); dt==0 && return p
    reorder_ecology!(p,sortperm(p.ids))
    n=length(p.ids); h=size(p.hidden,1); dims=w.grid.size
    cells=Vector{CartesianIndex{3}}(undef,n); counts=zeros(Int,dims)
    for i in 1:n
        for d in 1:3; p.positions[d,i]=dims[d]==1 ? 0f0 : mod(p.positions[d,i],Float32(dims[d])); end
        cells[i]=CartesianIndex(ntuple(d->floor(Int,p.positions[d,i])+1,3)); counts[cells[i]]+=1
    end
    children=Int[]; dead=Int[]
    for i in 1:n
        c=cells[i]
        neighbors=0
        for dy in -1:1, dx in -1:1
            neighbors+=counts[mod1(c[1]+dx,dims[1]),mod1(c[2]+dy,dims[2]),c[3]]
        end
        input=vcat(p.hidden[:,i],Float32[w.fields.nutrient[c],w.fields.deposit[c],p.storage[i],p.energy[i],neighbors/9,1])
        out=tanh.(p.weights[:,:,i]*input)
        alpha=Float32(-expm1(-Float64(dt)/p.response_time))
        p.hidden[:,i].+=(out[1:h].-p.hidden[:,i]).*alpha
        gates=(out[h+1:h+5].+1f0)./2f0
        uptake=min(w.fields.nutrient[c],Float32(dt)*p.uptake_rate*(0.2f0+0.8f0*gates[4]))
        w.fields.nutrient[c]-=uptake; p.storage[i]+=uptake
        # One unit of externally supplied chemical work per unit uptake; matter is retained.
        p.energy[i]+=uptake; p.energy_input+=uptake
        basal=min(p.energy[i],Float32(dt)*p.basal_cost); p.energy[i]-=basal; p.heat+=basal
        requested=Float32(dt)*p.work_rate
        work=min(p.energy[i],requested); p.energy[i]-=work; p.heat+=work
        funded=requested>0 ? work/requested : 0f0
        secretion=min(p.storage[i],Float32(dt)*p.secretion_rate*gates[4]*funded)
        p.storage[i]-=secretion; w.fields.deposit[c]+=secretion
        for d in 1:3
            drift=getproperty(w.fields,(:vx,:vy,:vz)[d])[c]
            p.positions[d,i]=dims[d]==1 ? 0f0 : mod(p.positions[d,i]+Float32(dt)*(p.speed*out[h+d]*funded+drift)/w.grid.spacing[d],Float32(dims[d]))
        end
        p.age[i]+=Float32(dt)
        if p.energy[i]<=0 || p.age[i]>=p.lifespan
            # Death returns all material to the occupied cell and dissipates remaining work.
            w.fields.nutrient[c]+=p.storage[i]; p.storage[i]=0; p.heat+=p.energy[i]; p.energy[i]=0; push!(dead,i); p.deaths+=1
        elseif p.storage[i]>=p.birth_material && p.energy[i]>=2p.birth_energy && gates[5]>0.5
            if n-length(dead)+length(children)>=p.capacity
                p.blocked_births+=1
            else
                p.energy[i]-=p.birth_energy; p.heat+=p.birth_energy
                p.storage[i]/=2; p.energy[i]/=2; push!(children,i)
            end
        end
    end
    for i in children
        p.events+=1
        # Birth event counter, seed, and persistent parent ID completely determine mutation.
        rng=MersenneTwister(string(p.seed,":",p.events,":",p.ids[i]))
        child=p.weights[:,:,i]+p.mutation.*randn(rng,Float32,h+5,h+6)
        p.weights=cat(p.weights,reshape(child,h+5,h+6,1);dims=3)
        p.positions=hcat(p.positions,p.positions[:,i]); p.hidden=hcat(p.hidden,p.hidden[:,i])
        push!(p.ids,p.next_id); p.next_id+=1; push!(p.parents,p.ids[i]); push!(p.lineages,p.lineages[i]); push!(p.generations,p.generations[i]+1)
        push!(p.storage,p.storage[i]); push!(p.energy,p.energy[i]); push!(p.age,0f0); p.births+=1
    end
    if !isempty(dead)
        keep=setdiff(collect(eachindex(p.ids)),dead)
        for name in (:ids,:parents,:lineages,:generations,:storage,:energy,:age); setfield!(p,name,getfield(p,name)[keep]); end
        p.positions=p.positions[:,keep]; p.hidden=p.hidden[:,keep]; p.weights=p.weights[:,:,keep]
    end
    p
end
advance!(w::World,p::Ecology,dt)=ecology_step!(w,p,dt)
