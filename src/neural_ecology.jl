"""Agents with inherited recurrent controllers, lifetime plasticity and spatial feedback.

Coordinates are periodic grid-cell coordinates. `weights` are the genotype;
`plastic_marks` and recurrent `hidden` state are reset at birth. Nine outputs are
movement xyz, uptake, membrane deposition, membrane erosion, signaling, division,
and signed phase remodeling. All active outputs share a paid work budget. A configurable uptake fraction
is catabolized to inert waste and supplies abstract chemical work; the remainder
is retained as structural storage. Recycling spent fuel cannot supply fresh work.
Attach an `ActiveMedium` through `medium` for phase, permeability and signal cues.
"""
mutable struct NeuralEcology
    ids::Vector{Int}
    parents::Vector{Int}
    lineages::Vector{Int}
    generations::Vector{Int}
    positions::Matrix{Float32}
    hidden::Matrix{Float32}
    storage::Vector{Float32}
    energy::Vector{Float32}
    weights::Array{Float32,3}
    plastic_marks::Array{Float32,3}
    last_actions::Matrix{Float32}
    age::Vector{Float32}
    medium::Any
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
    movement_work::Float64
    signal_work::Float64
    remodeling_work::Float64
    deposited::Float64
    eroded::Float64
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
    plasticity_rate::Float32
    plasticity_decay::Float32
    signal_rate::Float32
    remodeling_rate::Float32
    catabolic_fraction::Float32
    energy_yield::Float32
    neighbor_radius::Int
end
const NEURAL_ECOLOGY_SENSORS = 32
function NeuralEcology(positions::AbstractMatrix; medium=nothing, seed=17, capacity=2000,
        hidden=8, storage=.4f0, energy=1f0, uptake_rate=.5f0, basal_cost=.025f0,
        work_rate=.05f0, secretion_rate=.025f0, birth_material=.7f0, birth_energy=.3f0,
        mutation=.08f0, lifespan=80f0, speed=.3f0, response_time=1f0,
        plasticity_rate=.03f0, plasticity_decay=.02f0, signal_rate=.1f0,
        remodeling_rate=.03f0, catabolic_fraction=.5f0, energy_yield=2f0, neighbor_radius=2)
    n=size(positions,2); h=Int(hidden)
    size(positions,1)==3 || throw(DimensionMismatch("positions require three rows"))
    0<=n<=capacity && h>=2 && neighbor_radius>=1 || throw(ArgumentError("invalid capacity, hidden size or neighbor radius"))
    rates=Float32.((uptake_rate,basal_cost,work_rate,secretion_rate,birth_material,birth_energy,mutation,lifespan,speed,response_time,plasticity_rate,plasticity_decay,signal_rate,remodeling_rate,catabolic_fraction,energy_yield))
    all(x->isfinite(x)&&x>=0,rates) && 0<=catabolic_fraction<=1 && birth_material>0 && birth_energy>0 && lifespan>0 && response_time>0 && all(isfinite,positions) && all(x->isfinite(x)&&x>=0,(storage,energy)) || throw(ArgumentError("invalid neural ecology parameters"))
    rng=MersenneTwister(seed)
    weights=.18f0.*randn(rng,Float32,h+9,h+NEURAL_ECOLOGY_SENSORS,n)
    # Evolvable founder priors encourage feeding, turnover and directional food sensing.
    weights[h+4,end,:].=.8f0; weights[h+8,end,:].=.7f0
    for d in 1:3; weights[h+d,h+9+d,:].+=1f0; end
    NeuralEcology(collect(1:n),zeros(Int,n),collect(1:n),zeros(Int,n),Float32.(positions),
        zeros(Float32,h,n),fill(Float32(storage),n),fill(Float32(energy),n),weights,
        zeros(Float32,size(weights)),zeros(Float32,9,n),zeros(Float32,n),medium,
        Int(seed),0,n+1,Int(capacity),0,0,0,0.,Float64(n*Float32(energy)),0.,0.,0.,0.,0.,0.,rates...,Int(neighbor_radius))
end
function reorder_ecology!(p::NeuralEcology,order)
    sort(collect(order))==collect(eachindex(p.ids)) || throw(ArgumentError("invalid permutation"))
    _neural_select!(p,order)
end
function _neural_select!(p,indices)
    for name in (:ids,:parents,:lineages,:generations,:storage,:energy,:age)
        setfield!(p,name,getfield(p,name)[indices])
    end
    p.positions=p.positions[:,indices]; p.hidden=p.hidden[:,indices]
    p.weights=p.weights[:,:,indices]; p.plastic_marks=p.plastic_marks[:,:,indices]
    p.last_actions=p.last_actions[:,indices]; p
end
ecology_mass(w,p::NeuralEcology)=Float64(sum(Float64,w.fields.nutrient))+sum(Float64,w.fields.deposit)+sum(Float64,p.storage;init=0.)+p.waste
function ecology_statistics(p::NeuralEcology)
    (population=length(p.ids),births=p.births,deaths=p.deaths,blocked_births=p.blocked_births,
     max_generation=isempty(p.ids) ? 0 : maximum(p.generations),surviving_lineages=length(unique(p.lineages)),
     heat=p.heat,energy_input=p.energy_input,stored_energy=sum(Float64,p.energy;init=0.),
     energy_error=sum(Float64,p.energy;init=0.)+p.heat-p.energy_input,
     stored_material=sum(Float64,p.storage;init=0.),waste=p.waste,
     movement_work=p.movement_work,signal_work=p.signal_work,remodeling_work=p.remodeling_work,
     deposited=p.deposited,eroded=p.eroded,plastic_norm=norm(p.plastic_marks),
     mean_age=isempty(p.age) ? 0. : sum(Float64,p.age)/length(p.age))
end
function record_process_fields(p::NeuralEcology)
    names=Tuple(n for n in fieldnames(NeuralEcology) if n!=:medium)
    NamedTuple{names}(Tuple(getfield(p,n) for n in names))
end
function validate(w::World,p::NeuralEcology)
    n=length(p.ids); h=size(p.hidden,1)
    size(p.positions)==(3,n) && size(p.hidden)==(h,n) && size(p.weights)==(h+9,h+NEURAL_ECOLOGY_SENSORS,n) && size(p.plastic_marks)==size(p.weights) && size(p.last_actions)==(9,n) || throw(DimensionMismatch("neural ecology state"))
    n<=p.capacity && p.neighbor_radius>=1 && length(unique(p.ids))==n && all(<(p.next_id),p.ids) || throw(ArgumentError("invalid population IDs/capacity"))
    all(v->length(v)==n,(p.parents,p.lineages,p.generations,p.storage,p.energy,p.age)) || throw(DimensionMismatch("neural ecology vectors"))
    all(a->all(isfinite,a),(p.positions,p.hidden,p.weights,p.plastic_marks,p.last_actions)) && all(a->all(x->isfinite(x)&&x>=0,a),(p.storage,p.energy,p.age)) || throw(ArgumentError("invalid neural state"))
    rates=(p.uptake_rate,p.basal_cost,p.work_rate,p.secretion_rate,p.birth_material,p.birth_energy,p.mutation,p.lifespan,p.speed,p.response_time,p.plasticity_rate,p.plasticity_decay,p.signal_rate,p.remodeling_rate,p.catabolic_fraction,p.energy_yield)
    all(x->isfinite(x)&&x>=0,rates) && 0<=p.catabolic_fraction<=1 && p.birth_material>0 && p.birth_energy>0 && p.lifespan>0 && p.response_time>0 || throw(ArgumentError("invalid neural rates"))
    all(a->a isa Array{Float32,3}&&size(a)==w.grid.size,(w.fields.nutrient,w.fields.deposit)) || throw(ArgumentError("NeuralEcology requires matching CPU fields"))
    if p.medium!==nothing
        all(a->size(a)==w.grid.size,(_neural_medium_field(p,:phi),_neural_medium_field(p,:interface),p.medium.permeability,p.medium.signal,p.medium.signal_source,p.medium.remodeling)) || throw(DimensionMismatch("medium grids"))
    end
    all(x->isfinite(x)&&x>=0,(p.heat,p.energy_input,p.waste,p.movement_work,p.signal_work,p.remodeling_work,p.deposited,p.eroded)) || throw(ArgumentError("invalid neural ledgers"))
    for i in 1:n
        c=_neural_cell(view(p.positions,:,i),w.grid.size)
        all(x->isfinite(x)&&x>=0,(w.fields.nutrient[c],w.fields.deposit[c])) || throw(ArgumentError("invalid sampled resources"))
        if p.medium!==nothing
            perm=p.medium.permeability[c]
            isfinite(perm)&&0<=perm<=1 || throw(ArgumentError("permeability must be in [0,1]"))
        end
    end
    nothing
end
stability_limit(w,p::NeuralEcology)=1f0
@inline _neural_cell(pos,dims)=CartesianIndex(ntuple(d->dimsafe(pos[d],dims[d]),3))
@inline _neural_squash(x)=Float32(x/(1+abs(x)))
@inline function _neural_medium_field(p,name)
    m=p.medium
    name in (:phi,:interface) && !hasproperty(m,name) ? getproperty(m.phase,name) : getproperty(m,name)
end
@inline function _neural_medium(p,name,c,default=0f0)
    p.medium===nothing ? default : _neural_medium_field(p,name)[c]
end
# Sensors: local nutrient/membrane/signal/phase/interface/permeability; storage,
# energy/age; xyz central differences of each local channel; density and neighbor
# xyz displacement; bias. Bins bound search to a local cube, including z in 3D.
function _neural_input(w,p,i,c,bins,snapshot)
    dims=w.grid.size; h=size(p.hidden,1)
    localvalues(q)=Float32[w.fields.nutrient[q],w.fields.deposit[q],_neural_medium(p,:signal,q),_neural_medium(p,:phi,q),_neural_medium(p,:interface,q),_neural_medium(p,:permeability,q,1f0)]
    senses=_neural_squash.(localvalues(c))
    append!(senses,Float32[_neural_squash(p.storage[i]),_neural_squash(p.energy[i]),p.age[i]/p.lifespan])
    grads=zeros(Float32,6,3)
    for d in 1:3
        dims[d]==1 && continue
        lo=CartesianIndex(ntuple(k->k==d ? mod1(c[k]-1,dims[k]) : c[k],3))
        hi=CartesianIndex(ntuple(k->k==d ? mod1(c[k]+1,dims[k]) : c[k],3))
        grads[:,d].=_neural_squash.((localvalues(hi).-localvalues(lo))./2)
    end
    # Channel-major ordering: nutrient xyz first, then membrane xyz, etc.
    append!(senses,vec(permutedims(grads)))
    r=p.neighbor_radius; seen=Set{CartesianIndex{3}}(); displacement=zeros(Float32,3); count=0
    for dz in (dims[3]==1 ? (0:0) : (-r:r)), dy in -r:r, dx in -r:r
        q=CartesianIndex(mod1(c[1]+dx,dims[1]),mod1(c[2]+dy,dims[2]),mod1(c[3]+dz,dims[3]))
        q in seen && continue; push!(seen,q)
        for j in get(bins,q,Int[])
            j==i && continue
            delta=Float32[mod(snapshot[d,j]-snapshot[d,i]+dims[d]/2,dims[d])-dims[d]/2 for d in 1:3]
            dims[3]==1 && (delta[3]=0)
            sum(abs2,delta)>r*r && continue
            displacement.+=delta; count+=1
        end
    end
    push!(senses,_neural_squash(count)); append!(senses,displacement./Float32(max(count,1)*r)); push!(senses,1f0)
    vcat(p.hidden[:,i],senses)
end
# Proportional competition for pre-stage stock. Round each scarce share downward
# so no residual correction can favor a particular ID or overdraw the resource.
function _neural_allocate!(field,bins,requests)
    allocated=zeros(Float32,length(requests))
    for (c,occupants) in bins
        total=sum(i->Float64(requests[i]),occupants;init=0.)
        available=Float64(field[c]); total==0 && continue
        spent=0.
        for i in occupants
            exact=total<=available ? Float64(requests[i]) : available*(Float64(requests[i])/total)
            amount=Float32(exact)
            Float64(amount)>exact && (amount=prevfloat(amount))
            allocated[i]=amount; spent+=Float64(amount)
        end
        field[c]=Float32(max(0.,available-spent))
    end
    allocated
end
function _neural_budget(p,actions,delta,energy,active_dims)
    gates=(actions.+1f0)./2
    basal=min(energy,delta*p.basal_cost)
    remaining=energy-basal
    movement=sum(abs,actions[1:active_dims])/active_dims
    signaling=gates[7]; remodeling=abs(actions[9]); constructing=gates[5]+gates[6]
    demand=movement+signaling+remodeling+constructing
    request=delta*p.work_rate*demand; paid=min(remaining,request)
    funded=request>0 ? paid/request : 0f0
    (;gates,basal,movement,signaling,remodeling,constructing,demand,paid,funded)
end
function ecology_step!(w::World,p::NeuralEcology,dt::Real)
    validate(w,p); isfinite(dt)&&0<=dt<=1 || throw(ArgumentError("ecology dt must lie in [0,1]")); dt==0 && return p
    reorder_ecology!(p,sortperm(p.ids)); delta=Float32(dt); dims=w.grid.size; n=length(p.ids); h=size(p.hidden,1)
    for i in 1:n,d in 1:3; p.positions[d,i]=dims[d]==1 ? 0f0 : mod(p.positions[d,i],Float32(dims[d])); end
    snapshot=copy(p.positions); cells=[_neural_cell(view(snapshot,:,i),dims) for i in 1:n]
    bins=Dict{CartesianIndex{3},Vector{Int}}()
    for i in 1:n; push!(get!(bins,cells[i],Int[]),i); end
    # Observe all environmental cues before exchange: no sensor-order artifacts.
    inputs=[_neural_input(w,p,i,cells[i],bins,snapshot) for i in 1:n]
    outputs=[tanh.((p.weights[:,:,i].+p.plastic_marks[:,:,i])*inputs[i]) for i in 1:n]
    uptake_requests=Float32[delta*p.uptake_rate*((outputs[i][h+4]+1f0)/2)*_neural_medium(p,:permeability,cells[i],1f0) for i in 1:n]
    uptakes=_neural_allocate!(w.fields.nutrient,bins,uptake_requests)
    active_dims=dims[3]==1 ? 2 : 3
    budgets=[_neural_budget(p,outputs[i][h+1:h+9],delta,
        p.energy[i]+(uptakes[i]*p.catabolic_fraction)*p.energy_yield,active_dims) for i in 1:n]
    erosion_requests=Float32[delta*p.secretion_rate*budgets[i].gates[6]*budgets[i].funded for i in 1:n]
    erosions=_neural_allocate!(w.fields.deposit,bins,erosion_requests)
    children=Int[]; dead=Int[]
    for i in 1:n
        c=cells[i]; input=inputs[i]; out=outputs[i]; budget=budgets[i]
        alpha=Float32(-expm1(-Float64(delta)/p.response_time)); p.hidden[:,i].+=(out[1:h].-p.hidden[:,i]).*alpha
        actions=out[h+1:h+9]; p.last_actions[:,i].=actions; gates=budget.gates
        uptake=uptakes[i]
        # Catabolism spends finite chemical fuel into inert, accounted waste.
        # Structural recycling can recover matter but cannot recover spent work.
        catabolized=uptake*p.catabolic_fraction
        assimilated=uptake-catabolized
        chemical_work=catabolized*p.energy_yield
        p.storage[i]+=assimilated; p.waste+=catabolized
        p.energy[i]+=chemical_work; p.energy_input+=chemical_work
        basal=budget.basal; paid=budget.paid; funded=budget.funded
        movement=budget.movement; signaling=budget.signaling; remodeling=budget.remodeling
        constructing=budget.constructing; demand=budget.demand
        p.energy[i]-=basal; p.heat+=basal; p.energy[i]-=paid; p.heat+=paid
        if demand>0
            p.movement_work+=paid*movement/demand; p.signal_work+=paid*signaling/demand; p.remodeling_work+=paid*(remodeling+constructing)/demand
        end
        deposit=min(p.storage[i],delta*p.secretion_rate*gates[5]*funded)
        p.storage[i]-=deposit; w.fields.deposit[c]+=deposit; p.deposited+=deposit
        erode=erosions[i]
        p.storage[i]+=erode; p.eroded+=erode
        if p.medium!==nothing
            p.medium.signal_source[c]+=delta*p.signal_rate*signaling*funded
            p.medium.remodeling[c]+=delta*p.remodeling_rate*actions[9]*funded
        end
        proposal=copy(view(p.positions,:,i))
        for d in 1:3
            dims[d]==1 && continue
            proposal[d]=mod(proposal[d]+clamp(delta*p.speed*actions[d]*funded/w.grid.spacing[d],-.45f0,.45f0),Float32(dims[d]))
        end
        target=_neural_cell(proposal,dims); permeability=min(_neural_medium(p,:permeability,c,1f0),_neural_medium(p,:permeability,target,1f0))
        for d in 1:3
            dims[d]==1 && continue
            drift=getproperty(w.fields,(:vx,:vy,:vz)[d])[c]
            active=clamp(delta*p.speed*actions[d]*funded/w.grid.spacing[d],-.45f0,.45f0)
            flow=delta*drift/w.grid.spacing[d]
            p.positions[d,i]=mod(p.positions[d,i]+clamp(active*permeability+flow,-.45f0,.45f0),Float32(dims[d]))
        end
        # Reward-modulated bounded local Hebbian marks; genotype is untouched.
        reward=clamp(chemical_work-basal-paid,-1f0,1f0)
        decay=exp(-delta*p.plasticity_decay)
        p.plastic_marks[:,:,i].=clamp.(decay.*p.plastic_marks[:,:,i].+delta*p.plasticity_rate*reward.*(out*input'),-.5f0,.5f0)
        p.age[i]+=delta
        if p.energy[i]<=0 || p.age[i]>=p.lifespan
            occupied=_neural_cell(view(p.positions,:,i),dims)
            w.fields.nutrient[occupied]+=p.storage[i]; p.storage[i]=0; p.heat+=p.energy[i]; p.energy[i]=0; push!(dead,i); p.deaths+=1
        elseif p.storage[i]>=p.birth_material && p.energy[i]>=2p.birth_energy && gates[8]>.5f0
            if n-length(dead)+length(children)>=p.capacity
                p.blocked_births+=1
            else
                p.energy[i]-=p.birth_energy; p.heat+=p.birth_energy; p.storage[i]/=2; p.energy[i]/=2; push!(children,i)
            end
        end
    end
    for i in children
        p.events+=1
        rng=MersenneTwister(string(p.seed,":",p.events,":",p.ids[i]))
        child=p.weights[:,:,i]+p.mutation.*randn(rng,Float32,h+9,h+NEURAL_ECOLOGY_SENSORS)
        p.weights=cat(p.weights,reshape(child,size(child)...,1);dims=3)
        p.plastic_marks=cat(p.plastic_marks,zeros(Float32,size(child)...,1);dims=3)
        pos=copy(view(p.positions,:,i))
        for d in 1:3; pos[d]=dims[d]==1 ? 0f0 : mod(pos[d]+.35f0*randn(rng,Float32),Float32(dims[d])); end
        p.positions=hcat(p.positions,pos); p.hidden=hcat(p.hidden,zeros(Float32,h)); p.last_actions=hcat(p.last_actions,zeros(Float32,9))
        push!(p.ids,p.next_id); p.next_id+=1; push!(p.parents,p.ids[i]); push!(p.lineages,p.lineages[i]); push!(p.generations,p.generations[i]+1)
        push!(p.storage,p.storage[i]); push!(p.energy,p.energy[i]); push!(p.age,0f0); p.births+=1
    end
    !isempty(dead) && _neural_select!(p,setdiff(collect(eachindex(p.ids)),dead))
    p
end
advance!(w::World,p::NeuralEcology,dt)=ecology_step!(w,p,dt)
