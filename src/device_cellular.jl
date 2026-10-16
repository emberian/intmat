"""Experimental capacity-major cellular execution on a KA backend.

This uses Float32 control with opt-in Float64 battery and physical inventory
storage on CPU/AMD; it is not bit-identical CPU CellularEcology.
`inventory_type=Float64` widens biomass, metabolites, soluble fields and their
transfer buffers; it requires `energy_type=Float64`. Geometry and controllers
remain Float32. This explicit choice belongs to new source-pinned experiments.
`support_execution=:parallel` explicitly parallelizes materialized support
construction with extra candidate scratch and exactly the serial edge/sum order.
The default `:serial` allocates no support scratch; `:parallel` rejects axis support.
Its physiology/controller contracts are shared; junctions use synchronous donor
limiting, bonds use mutual bounded-degree proposals, and mechanics uses checked
fixed substeps with conservatively funded motors. Cell state, genomes, resource
fields, contact bins and lifecycle metadata stay on the device between reads.
"""
mutable struct DeviceCellular{S,E,B,Q,G,O}
    state::S
    environment::E
    backend::B
    parameters::Q
    geometry::G
    attempted_steps::Int
    ownership::O
end
isdefined(@__MODULE__,:cellular_parameters) || include("cellular_initialization.jl")
include("device_cellular_support.jl")

function device_energy_type(backend,energy_type)
    energy_type in (Float32,Float64) || throw(ArgumentError("energy_type must be Float32 or Float64"))
    energy_type===Float64 && occursin("Metal",string(typeof(backend))) && throw(ArgumentError("Float64 energy is unsupported on Metal"))
    energy_type
end
function device_cellular_arrays(backend,C,h,K,S,V,bins;support_strategy=:materialized,axis_capacity=0,energy_type=Float32,inventory_type=Float32)
    device_energy_type(backend,energy_type)
    device_energy_type(backend,inventory_type)
    inventory_type===Float64 && energy_type!==Float64 && throw(ArgumentError("Float64 inventory requires Float64 batteries"))
    z(T,d...)=KernelAbstractions.zeros(backend,T,d...)
    state=(active=z(Int32,C),ids=z(Int32,C),parents=z(Int32,C),lineages=z(Int32,C),generations=z(Int32,C),
        position=z(Float32,C,3),position_next=z(Float32,C,3),mass=z(inventory_type,C),metabolites=z(inventory_type,C,2),
        metabolites_next=z(inventory_type,C,2),energy=z(energy_type,C),age=z(Float32,C),
        hidden=z(Float32,C,h),adaptation=z(Float32,C,h),weights=z(Float32,C,2h+CELLULAR_SENSORS,h+CELLULAR_ACTIONS),
        marks=z(Float32,C,2h+CELLULAR_SENSORS,h+CELLULAR_ACTIONS),actions=z(Float32,C,CELLULAR_ACTIONS),
        inputs=z(Float32,C,2h+CELLULAR_SENSORS),outputs=z(Float32,C,h+CELLULAR_ACTIONS),
        radii=z(Float32,C),enzymes=z(Float32,C,4),requests=z(inventory_type,C,2),grants=z(inventory_type,C,2),
        sites=z(Int32,C,support_strategy===:axis ? 0 : S),support=z(Float32,C,support_strategy===:axis ? 0 : S),support_count=z(Int32,C),
        bonds=z(Int32,C,K),bond_ids=z(Int32,C,K),rest=z(Float32,C,K),candidates=z(Int32,C,K),candidate_rest=z(Float32,C,K),
        force=z(Float32,C,3),rows=z(Float32,C),exchange_scale=z(inventory_type,C,2),waste=z(inventory_type,C),
        ledger=z(Float32,C,9),ledger_low=z(Float32,C,9),cell_bin=z(Int32,C),next=z(Int32,C),new_parent=z(Int32,C))
    if support_strategy===:axis
        names=filter(k->!(k in (:sites,:support,:support_count)),keys(state))
        state=merge(NamedTuple{names}(state),(axis_indices=z(Int32,C,axis_capacity,3),
            axis_distance2=z(Float32,C,axis_capacity,3),support_normalization=z(Float32,C),
            support_fallback=z(Int32,C),support_count=z(Int32,C)))
    elseif support_strategy!==:materialized
        throw(ArgumentError("unknown device support strategy"))
    end
    env=(fields=z(inventory_type,V,3),scratch=z(inventory_type,V,3),demand=z(inventory_type,V,2),scale=z(inventory_type,V,2),
        returns=z(inventory_type,V),permeability=z(Float32,V),occupancy=z(Float32,V),
        source_masks=z(Float32,V,3),sink_masks=z(Float32,V,3),resource_ledger=z(inventory_type,3,6),resource_ledger_low=z(inventory_type,3,6),
        inert=z(Float32,V),inert_low=z(Float32,V),observation=z(Float32,12),observation_min=z(Int32,2),
        quota=z(Float32,3,2),heads=z(Int32,prod(bins)),metadata=z(Int32,10),status=z(Int32,1),
        heap=z(Int32,C),heap_keys=z(UInt32,C),free=z(Int32,C))
    state,env
end
DeviceCellular(w::World,reference::CellularEcology;kwargs...)=device_cellular_import(w,reference;kwargs...)
function device_cellular_import(w::World,reference;capacity=reference.capacity,
        radius_bound=2f0,bin_width=radius_bound/2,mechanical_substeps=16,seed=reference.seed,support_strategy=:materialized,support_execution=:serial,energy_type=Float32,inventory_type=Float32)
    device_energy_type(w.backend,energy_type)
    support_strategy in (:materialized,:axis) || throw(ArgumentError("unknown device support strategy"))
    device_support_execution_validate(support_execution,support_strategy)
    reference.medium===nothing && reference.resources.medium===nothing || throw(ArgumentError("imported device medium coupling is not implemented; construct the explicit device composition"))
    all(x->length(x)==3,(reference.resources.fields,reference.resources.potentials,reference.resources.diffusion,reference.resources.source_rates,reference.resources.sink_rates,reference.resources.source_masks,reference.resources.sink_masks,reference.resources.source_remaining,reference.resources.sink_remaining)) && reference.resources.potentials[3]==0 || throw(ArgumentError("device cells require two precursors and one zero-potential waste species"))
    (!hasproperty(w.fields,:nutrient)||all(iszero,Array(w.fields.nutrient))) || throw(ArgumentError("DeviceCellular owns its resource pool; zero unrelated world nutrient explicitly"))
    C=Int(capacity);h=size(reference.hidden,1);n=length(reference.ids);K=reference.parameters.max_degree
    0<=n<=C && C>0 && C<typemax(Int32) || throw(ArgumentError("invalid fixed device capacity"))
    radius_bound>0 && mechanical_substeps>0 || throw(ArgumentError("positive radius bound/substeps required"))
    reference.stencil.grid.size==w.grid.size && reference.stencil.grid.extent==w.grid.extent || throw(ArgumentError("reference/device grids differ"))
    all(a->size(a)==w.grid.size && all(x->isfinite(x)&&x>=0&&isfinite(Float32(x)),a),vcat(reference.resources.fields,reference.resources.source_masks,reference.resources.sink_masks)) || throw(ArgumentError("invalid reference resource fields"))
    budget_source=inventory_type.(reference.resources.source_remaining);budget_sink=inventory_type.(reference.resources.sink_remaining)
    all(x->isfinite(x)&&x>=0,vcat(budget_source,budget_sink)) || throw(ArgumentError("resource budgets must remain finite in Float32"))
    dims=w.grid.size;spacing=w.grid.spacing;V=prod(dims);nd=dims[3]==1 ? 2 : 3
    exchange_radius=Float32(reference.parameters.exchange_radius)
    isfinite(exchange_radius)&&exchange_radius>0 || throw(ArgumentError("exchange radius must remain positive finite in Float32"))
    candidates=ntuple(d->d>nd ? 1 : min(dims[d],2ceil(Int,exchange_radius/spacing[d])+2),3)
    S=prod(candidates); binwidth=Float32(bin_width)
    isfinite(binwidth)&&binwidth>0 || throw(ArgumentError("positive finite device bin width required"))
    bins=ntuple(d->d>nd ? 1 : max(1,floor(Int,w.grid.extent[d]/binwidth)),3)
    prod(BigInt.(bins))<typemax(Int32)÷2 || throw(ArgumentError("device bin indexing capacity exceeded"))
    volume=prod(Float32(spacing[d]) for d in 1:nd)
    isfinite(volume)&&volume>0 || throw(ArgumentError("device voxel volume must remain finite and positive"))
    q=merge(map(x->x isa Integer ? Int32(x) : Float32(x),reference.parameters),
        (potentials=Tuple(Float32.(reference.resources.potentials)),diffusion=Tuple(Float32.(reference.resources.diffusion)),
        source_rates=Tuple(Float32.(reference.resources.source_rates)),sink_rates=Tuple(Float32.(reference.resources.sink_rates)),
        source_totals=hasproperty(reference.resources,:source_totals) ? Tuple(Float32.(reference.resources.source_totals)) : Tuple(Float32(sum(Float64,m)) for m in reference.resources.source_masks),
        sink_totals=hasproperty(reference.resources,:sink_totals) ? Tuple(Float32.(reference.resources.sink_totals)) : Tuple(Float32(sum(Float64,m)) for m in reference.resources.sink_masks),
        radius_bound=Float32(radius_bound),bin_width=binwidth,mechanical_substeps=Int32(mechanical_substeps),seed=UInt32(seed)))
    all(x->isfinite(x)&&x>=0,Iterators.flatten(v isa Tuple ? v : (v,) for v in values(q))) || throw(ArgumentError("parameters must remain finite in Float32"))
    all(>(0),(q.exchange_radius,q.density,q.store_capacity,q.energy_capacity,q.half_saturation,q.response_time,q.adaptation_time,q.motor_cost,q.minimum_mass,q.division_mass)) || throw(ArgumentError("required rates/scales must remain positive in Float32"))
    geom=(dims=Int32.(dims),spacing=spacing,extent=w.grid.extent,nd=Int32(nd),volume=volume,
        candidates=Int32.(candidates),bins=Int32.(bins),capacity=Int32(C),hidden=Int32(h),degree=Int32(K),
        support_strategy=support_strategy===:axis ? Int32(1) : Int32(0),support_radius2=exchange_radius^2,energy_bits=Int32(8sizeof(energy_type)),inventory_bits=Int32(8sizeof(inventory_type)))
    state,env=device_cellular_arrays(w.backend,C,h,K,S,V,bins;support_strategy,axis_capacity=maximum(candidates),energy_type,inventory_type)
    # Device buffers are already zeroed; stage founders only, never a full-capacity host replica.
    host=map(a->zeros(eltype(a),n,Base.tail(size(a))...),state)
    host.active[1:n].=1;host.ids[1:n].=Int32.(reference.ids);host.parents[1:n].=Int32.(reference.parents)
    host.lineages[1:n].=Int32.(reference.lineages);host.generations[1:n].=Int32.(reference.generations)
    for (dst,src) in ((:position,:positions),(:metabolites,:metabolites),(:hidden,:hidden),(:adaptation,:adaptation),(:actions,:last_actions))
        getproperty(host,dst)[1:n,:].=eltype(getproperty(host,dst)).(permutedims(getproperty(reference,src)))
    end
    for (dst,src) in ((:mass,:biomass),(:energy,:energy),(:age,:age))
        getproperty(host,dst)[1:n].=eltype(getproperty(host,dst)).(getproperty(reference,src))
    end
    host.weights[1:n,:,:].=permutedims(reference.weights,(3,2,1))
    host.marks[1:n,:,:].=permutedims(reference.plastic_marks,(3,2,1))
    index=Dict(id=>i for (i,id) in enumerate(reference.ids));degree=zeros(Int,n)
    for (edge,rest) in zip(reference.bonds.endpoints,reference.bonds.rest_lengths)
        i=index[edge[1]];j=index[edge[2]]
        degree[i]+=1;degree[j]+=1
        max(degree[i],degree[j])<=K || throw(ArgumentError("initial bonds exceed declared maximum degree"))
        for (a,b) in ((i,j),(j,i))
            host.bonds[a,degree[a]]=b;host.bond_ids[a,degree[a]]=reference.ids[b];host.rest[a,degree[a]]=rest
        end
    end
    for value in values(host)
        eltype(value)<:AbstractFloat && !all(isfinite,value) && throw(ArgumentError("reference state must remain finite in Float32"))
    end
    all(x->isfinite(Float32(x)),host.energy) || throw(ArgumentError("battery must remain representable in Float32 controller observations"))
    all(>(0),host.mass[1:n]) && all(>=(0),host.metabolites) && all(>=(0),host.energy) && all(>=(0),host.age) || throw(ArgumentError("reference physiological state must remain valid after Float32 conversion"))
    0<=q.efficiency<=1 && 0<=q.mutation_probability<=1 || throw(ArgumentError("invalid device probability/efficiency"))
    for name in keys(state)
        dst=getproperty(state,name);src=getproperty(host,name)
        copyto!(view(dst,1:n,ntuple(_->Colon(),ndims(dst)-1)...),src)
    end
    for s in 1:3
        copyto!(view(env.fields,:,s),vec(reference.resources.fields[s]))
        copyto!(view(env.source_masks,:,s),vec(reference.resources.source_masks[s]))
        copyto!(view(env.sink_masks,:,s),vec(reference.resources.sink_masks[s]))
    end
    ledger=zeros(inventory_type,3,6);ledger[:,1].=budget_source;ledger[:,2].=budget_sink
    ledgerlow=zeros(inventory_type,3,6)
    ledgerlow[:,1].=inventory_type.(reference.resources.source_remaining.-Float64.(budget_source))
    ledgerlow[:,2].=inventory_type.(reference.resources.sink_remaining.-Float64.(budget_sink))
    copyto!(env.resource_ledger_low,ledgerlow)
    copyto!(env.resource_ledger,ledger);fill!(env.permeability,1f0)
    metadata=zeros(Int32,10);metadata[1]=reference.next_id;metadata[2]=n;copyto!(env.metadata,metadata)
    p=DeviceCellular(state,env,w.backend,q,geom,0,device_cellular_contract(state,env,q,geom;lifecycle_highwater=n,support_execution))
    device_cellular_launch!(p,dc_geometry!,p.state,p.environment,q,geom;ndrange=C)
    device_cellular_launch!(p,dc_support!,p.state,p.environment,q,geom;ndrange=C)
    device_cellular_launch!(p,dc_porosity_scatter!,p.state,p.environment,geom;ndrange=C)
    device_cellular_launch!(p,dc_porosity_commit!,p.environment,q;ndrange=V)
    device_cellular_status!(p)
    p
end

function device_cellular_contract(state,env,q,g;resource_workspace=nothing,resource_transient_aliases=(),
        lifecycle_workspace=nothing,lifecycle_highwater=nothing,support_execution=:serial,support_workspace=nothing)
    mechanics=(state=NamedTuple{(:active,:ids,:position,:position_next,:mass,:metabolites,:radii,:next,:cell_bin,:bonds,:bond_ids,:rest,:actions,:force,:rows,:energy,:ledger,:ledger_low)}(state),
        environment=NamedTuple{(:status,:heads,:metadata)}(env),
        parameters=NamedTuple{(:density,:radius_bound,:repulsion,:stiffness,:speed,:motor_cost)}(q),
        geometry=NamedTuple{(:nd,:extent,:bins,:degree)}(g))
    controller=(state=NamedTuple{(:active,:weights,:marks,:inputs,:outputs)}(state),
        environment=(status=env.status,),geometry=(hidden=g.hidden,))
    resource=resource_workspace===nothing ? DeviceResourceWorkspace(KernelAbstractions.get_backend(env.fields),size(env.fields,1);inventory_type=eltype(env.fields)) : resource_workspace
    device_resource_workspace_validate(resource)
    size(resource.grants)==(size(env.fields,1),6) || throw(ArgumentError("resource workspace grid differs"))
    device_resource_workspace_nonalias(resource,state,env;allowed_transients=resource_transient_aliases)
    if lifecycle_workspace===nothing
        if lifecycle_highwater===nothing
            host=(ids=Array(state.ids),active=Array(state.active),new_parent=Array(state.new_parent))
            lifecycle_highwater=device_lifecycle_host_highwater(host)
        end
        lifecycle_workspace=DeviceLifecycleWorkspace(KernelAbstractions.get_backend(env.fields),Int(g.capacity);initial=lifecycle_highwater)
    end
    device_lifecycle_workspace_validate(lifecycle_workspace,KernelAbstractions.get_backend(env.fields),Int(g.capacity))
    device_lifecycle_workspace_nonalias(lifecycle_workspace,state,env,resource)
    strategy=g.support_strategy==1 ? :axis : :materialized
    device_support_execution_validate(support_execution,strategy)
    support=support_execution===:parallel ? (support_workspace===nothing ? DeviceSupportWorkspace(KernelAbstractions.get_backend(env.fields),state,g) : support_workspace) : nothing
    support_execution===:serial && support_workspace!==nothing && throw(ArgumentError("serial support cannot own parallel scratch"))
    device_support_workspace_validate(support,KernelAbstractions.get_backend(env.fields),g)
    device_support_workspace_nonalias(support,state,env,resource,lifecycle_workspace)
    (state=state,environment=env,parameters=q,geometry=g,
        state_shapes=map(size,state),environment_shapes=map(size,env),resource_transient_aliases=resource_transient_aliases,
        stages=(mechanics=mechanics,controller=controller,resource=resource,lifecycle=lifecycle_workspace,support=support))
end
"""Cheap identity/shape preflight. Direct array-content edits remain caller responsibility."""
function device_cellular_validate(w::World,p::DeviceCellular)
    o=p.ownership
    p.state===o.state && p.environment===o.environment && p.parameters===o.parameters && p.geometry===o.geometry ||
        throw(ArgumentError("device cellular array roles/parameters were replaced; construct a new process"))
    map(size,p.state)==o.state_shapes && map(size,p.environment)==o.environment_shapes ||
        throw(ArgumentError("device cellular owned array shapes changed"))
    device_resource_workspace_validate(o.stages.resource)
    device_lifecycle_workspace_validate(o.stages.lifecycle,p.backend,Int(p.geometry.capacity))
    device_support_workspace_validate(o.stages.support,p.backend,p.geometry)
    w.backend==p.backend && w.grid.size==Tuple(Int.(p.geometry.dims)) && w.grid.extent==p.geometry.extent ||
        throw(ArgumentError("device cellular world mismatch"))
    true
end
function device_cellular_launch!(p,kernel,args...;ndrange)
    if kernel===dc_support!
        ndrange==Int(p.geometry.capacity) || throw(ArgumentError("support launch capacity differs"))
        device_support_launch!(p,device_support_workspace(p),args...)
    else
        kernel(p.backend,128)(args...;ndrange)
    end
end
function device_cellular_status!(p::DeviceCellular)
    KernelAbstractions.synchronize(p.backend)
    code=only(Array(p.environment.status))
    code==0 || throw(ArgumentError("DeviceCellular halted with status bits $code (1 radius/support, 2 mechanics resolution, 4 identity overflow, 8 nonfinite, 16 accounting precision). State is not a complete step; observe diagnostic arrays; this partial state cannot be resumed."))
    p
end
"""Shared bound applied to the actual rounded Float32 submitted timestep."""
function device_cellular_timestep_limit(p::DeviceCellular)
    q=p.parameters;g=p.geometry
    D=maximum(q.diffusion);S=sum(inv(Float64(g.spacing[d])^2) for d in 1:Int(g.nd))
    min(.25,Float64(q.response_time/2),D>0 ? .45/(D*S) : Inf)
end
# Typed, optional anabolic/habitat extension. Concrete device descriptors and
# host launch callbacks are supplied by the composition owner; Nothing preserves
# the base physiology. Every custom launched kernel must honor the shared latch.
device_anabolic_argument(::Nothing)=nothing
device_anabolic_begin!(::Nothing,w,p,dt)=nothing
device_anabolic_sense!(::Nothing,w,p,dt)=nothing
device_controller_override!(anabolism,w,p,dt)=nothing
device_anabolic_prepare!(::Nothing,w,p,dt)=nothing
device_anabolic_inherit!(::Nothing,w,p,step)=nothing
device_anabolic_finish!(::Nothing,w,p,dt)=nothing
function device_cellular_step!(w::World,p::DeviceCellular,dt::Real;check_status=true,execution=:projected,anabolism=nothing,lifecycle=:prefix)
    device_cellular_validate(w,p)
    execution in (:staged,:projected) || throw(ArgumentError("unknown device execution recipe"))
    lifecycle in (:prefix,:dense) || throw(ArgumentError("unknown device lifecycle recipe"))
    delta=Float32(dt)
    isfinite(delta) && 0<delta<=device_cellular_timestep_limit(p) || throw(ArgumentError("invalid device cellular timestep"))
    p.attempted_steps+=1
    device_cellular_enqueue!(w,p,delta,execution,anabolism,UInt32(p.attempted_steps),lifecycle)
    check_status && device_cellular_status!(p)
    p
end
# Internal fixed-capacity launch schedule. The checked ordinary entry point owns
# host attempt bookkeeping. A captured queue may instead supply a device clock;
# neither capture nor this function advances host counters or World time.
function device_cellular_enqueue!(w,p,delta::Float32,execution::Symbol,anabolism,step_index,lifecycle::Symbol=:prefix)
    g=p.geometry;q=p.parameters;s=p.state;e=p.environment
    C=Int(g.capacity);V=prod(Int.(g.dims))
    # All arrays are fixed capacity. Latched failure causes subsequent kernels to no-op.
    device_anabolic_begin!(anabolism,w,p,delta)
    device_resource_phase!(p,e,q,g,delta)
    device_cellular_launch!(p,dc_diffuse!,e,q,g,delta;ndrange=(V,3))
    device_cellular_launch!(p,dc_commit_resources!,e;ndrange=(V,3))
    device_cellular_launch!(p,dc_clear_grid!,e;ndrange=V)
    device_cellular_launch!(p,dc_clear_heads!,e;ndrange=length(e.heads))
    device_cellular_launch!(p,dc_geometry!,s,e,q,g;ndrange=C)
    device_cellular_launch!(p,dc_support!,s,e,q,g;ndrange=C)
    device_cellular_launch!(p,dc_sense!,s,e,q,g;ndrange=C)
    device_anabolic_sense!(anabolism,w,p,delta)
    if execution===:staged
        device_cellular_launch!(p,dc_controller!,s,e,g;ndrange=(C,Int(g.hidden)+CELLULAR_ACTIONS))
    else
        controller=p.ownership.stages.controller
        device_cellular_launch!(p,dc_controller!,controller.state,controller.environment,controller.geometry;ndrange=(C,Int(g.hidden)+CELLULAR_ACTIONS))
    end
    device_controller_override!(anabolism,w,p,delta)
    device_cellular_launch!(p,dc_regulate!,s,e,q,g,delta;ndrange=C)
    device_cellular_launch!(p,dc_allocate_sites!,e;ndrange=(V,2))
    device_cellular_launch!(p,dc_gather!,s,e,g;ndrange=C)
    # Keep the junction stage ABI small: whole-state tuples trigger a pathological
    # LLVM EarlyCSE alias walk in the pinned AMD compiler for some 3D specializations.
    junction=(active=s.active,ids=s.ids,bonds=s.bonds,bond_ids=s.bond_ids,mass=s.mass,
        metabolites=s.metabolites,metabolites_next=s.metabolites_next,actions=s.actions,
        exchange_scale=s.exchange_scale,ledger=s.ledger,ledger_low=s.ledger_low)
    junction_status=(status=e.status,);junction_q=(junction_rate=q.junction_rate,);junction_g=(degree=g.degree,)
    device_cellular_launch!(p,dc_exchange_scale!,junction,junction_status,junction_q,junction_g,delta;ndrange=C)
    device_cellular_launch!(p,dc_exchange_commit!,junction,junction_status,junction_q,junction_g,delta;ndrange=C)
    device_cellular_launch!(p,dc_commit_metabolites!,s,e;ndrange=C)
    device_anabolic_prepare!(anabolism,w,p,delta)
    device_cellular_launch!(p,dc_metabolism!,s,e,q,g,device_anabolic_argument(anabolism),delta;ndrange=C)
    device_cellular_launch!(p,dc_return_waste!,s,e,g;ndrange=C)
    device_cellular_launch!(p,dc_commit_waste!,e;ndrange=V)
    if execution===:staged
        device_cellular_mechanics!(p,delta,Val(false))
    else
        device_cellular_mechanics!(p,delta,Val(true))
    end
    device_cellular_launch!(p,dc_propose_bonds!,s,e,q,g;ndrange=C)
    device_cellular_launch!(p,dc_commit_bonds!,s,e,q,g;ndrange=C)
    highwater=device_lifecycle_workspace(p).highwater
    if lifecycle===:dense
        device_cellular_launch!(p,dc_lifecycle_dense!,s,e,q,g,step_index,highwater;ndrange=1)
    else
        device_cellular_launch!(p,dc_lifecycle_prefix!,s,e,q,g,step_index,highwater;ndrange=1)
    end
    device_cellular_launch!(p,dc_divide!,s,e,q,g,step_index;ndrange=C)
    device_anabolic_inherit!(anabolism,w,p,step_index)
    device_cellular_launch!(p,dc_clear_occupancy!,e;ndrange=V)
    device_cellular_launch!(p,dc_support!,s,e,q,g;ndrange=C)
    device_cellular_launch!(p,dc_porosity_scatter!,s,e,g;ndrange=C)
    device_cellular_launch!(p,dc_porosity_commit!,e,q;ndrange=V)
    device_anabolic_finish!(anabolism,w,p,delta)
    device_cellular_launch!(p,dc_validate_complete!,s,e,g;ndrange=max(C,V))
    device_cellular_launch!(p,dc_complete_step!,e;ndrange=1)
    p
end
function device_cellular_mechanics!(p,delta,::Val{projected}) where {projected}
    stage=projected ? p.ownership.stages.mechanics : (state=p.state,environment=p.environment,parameters=p.parameters,geometry=p.geometry)
    s=stage.state;e=stage.environment;q=stage.parameters;g=stage.geometry
    C=Int(p.geometry.capacity);substeps=Int(p.parameters.mechanical_substeps);dt=delta/Float32(substeps)
    for _ in 1:substeps
        device_cellular_launch!(p,dc_clear_heads!,e;ndrange=length(e.heads))
        device_cellular_launch!(p,dc_geometry!,s,e,q,g;ndrange=C)
        device_cellular_launch!(p,dc_forces!,s,e,q,g,dt;ndrange=C)
        device_cellular_launch!(p,dc_move!,s,e,q,g,dt;ndrange=C)
        device_cellular_launch!(p,dc_commit_positions!,s,e;ndrange=C)
    end
    device_cellular_launch!(p,dc_clear_heads!,e;ndrange=length(e.heads))
    device_cellular_launch!(p,dc_geometry!,s,e,q,g;ndrange=C)
    nothing
end
advance!(w::World,p::DeviceCellular,dt)=device_cellular_step!(w,p,dt)
function device_cellular_snapshot(p::DeviceCellular;check_status=true)
    check_status && device_cellular_status!(p)
    (state=map(Array,p.state),environment=map(Array,p.environment),parameters=p.parameters,
        geometry=p.geometry,attempted_steps=p.attempted_steps)
end
"""Compact GPU-reduced monitoring, not a full genome/scratch download.
Float32 reduction roundoff applies; use snapshot + Float64 sums for strict audits.
"""
function device_cellular_statistics(p::DeviceCellular)
    device_cellular_status!(p);e=p.environment;s=p.state;g=p.geometry;q=p.parameters
    fill!(e.observation,0f0);fill!(e.observation_min,reinterpret(Int32,Inf32))
    device_cellular_launch!(p,dc_observe_cells!,s,e,q,g;ndrange=Int(g.capacity))
    device_cellular_launch!(p,dc_observe_grid!,e,q;ndrange=prod(Int.(g.dims)))
    KernelAbstractions.synchronize(p.backend)
    o=Float64.(Array(e.observation));m=Array(e.metadata)
    ledger=Float64.(Array(e.resource_ledger)).+Float64.(Array(e.resource_ledger_low))
    mins=reinterpret(Float32,Array(e.observation_min));population=Int(m[2])
    (population=population,births=Int(m[3]),deaths=Int(m[4]),blocked_births=Int(m[5]),completed_steps=Int(m[6]),
        bonds=round(Int,o[7]),material=o[1]+o[3],potential=o[2]+o[4],
        supplied=sum(ledger[:,3]),removed=sum(ledger[:,4]),
        supplied_potential=sum(ledger[k,3]*q.potentials[k] for k in 1:3),
        removed_potential=sum(ledger[k,4]*q.potentials[k] for k in 1:3),work=o[5],heat=o[6],
        synthesized=o[8],oxidized=o[9],finite=all(isfinite,o)&&o[10]==0,minimum_resource=mins[1],
        minimum_energy=population==0 ? 0f0 : mins[2],observation_bytes=4*(12+2+10+36+1))
end

include("device_cellular_kernels.jl")
include("device_resource_exchange.jl")
include("device_cellular_lifecycle.jl")

"""Float64 host audit of an explicit full snapshot; this intentionally transfers all state."""
function device_cellular_audit(snapshot)
    s=snapshot.state;e=snapshot.environment;q=snapshot.parameters
    active=findall(==(1),s.active)
    ledger=Float64.(e.resource_ledger).+Float64.(e.resource_ledger_low)
    cellledger=Float64.(s.ledger).+Float64.(s.ledger_low)
    material=sum(Float64,e.fields)+sum(Float64,e.inert)+sum(Float64,e.inert_low)+
        sum(Float64,s.mass[active];init=0.0)+sum(Float64,s.metabolites[active,:];init=0.0)
    potential=sum(sum(Float64,e.fields[:,k])*q.potentials[k] for k in 1:3)+
        sum(sum(Float64,s.metabolites[active,k];init=0.0)*q.potentials[k] for k in 1:2)+sum(Float64,s.energy[active];init=0.0)
    (material=material,potential=potential,work=sum(cellledger[:,1]),heat=sum(cellledger[:,2]),
        supplied=sum(ledger[:,3]),removed=sum(ledger[:,4]),
        supplied_potential=sum(ledger[k,3]*q.potentials[k] for k in 1:3),
        removed_potential=sum(ledger[k,4]*q.potentials[k] for k in 1:3),
        population=length(active),completed_steps=Int(e.metadata[6]),
        finite=all(a->!(eltype(a)<:AbstractFloat)||all(isfinite,a),values(s))&&
            all(a->!(eltype(a)<:AbstractFloat)||all(isfinite,a),values(e)))
end

include("device_cellular_checkpoint.jl")

include("device_cellular_initialization.jl")
