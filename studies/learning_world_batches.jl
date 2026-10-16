# Independent policies share the existing constructive batch engine and execute
# their original host hooks through a narrow execution-only launch seam.
using Interior, SHA, TOML
using Interior: KernelAbstractions
isdefined(@__MODULE__,:ConstructiveBatchCells)||include("constructive_world_batches.jl")
isdefined(@__MODULE__,:NativeLifetimeCredit)||include("native_lifetime_credit.jl")

struct WorldBatchScalar{A}
    values::A
end
WorldBatchAdapt.adapt_structure(to,x::WorldBatchScalar)=WorldBatchScalar(WorldBatchAdapt.adapt(to,x.values))
@inline batch_world(x::WorldBatchScalar,w)=batch_world(x.values,w)[1]

struct SidecarBatchPolicy{T,L,K,O,P,C,N,R,F,S}
    name::Symbol
    sidecar_type::T
    launch_module::L
    kernels::K
    member_options::O
    member_properties::P
    cohort_properties::C
    cache_names::N
    repack::R
    restore::F
    source::S
end
const WORLD_SIDECAR_POLICIES=Dict{Symbol,Any}()
const WORLD_SIDECAR_KERNELS=IdDict()
const LEARNING_MEMBER_OPTIONS=(:rate,:sigma,:kappa,:seed,:learning,:critic,:critic_learning)
function sidecar_module_ast(ex,name)
    ex isa Expr||return nothing
    ex.head==:module&&ex.args[2]==name&&return ex
    for child in ex.args
        found=sidecar_module_ast(child,name);found===nothing||return found
    end
    nothing
end
function project_sidecar_kernel!(kernel,scope,path)
    haskey(WORLD_SIDECAR_KERNELS,kernel)&&return WORLD_SIDECAR_KERNELS[kernel]
    moduleexpr=sidecar_module_ast(Meta.parseall(read(path,String)),nameof(scope))
    moduleexpr===nothing&&error("sidecar module changed")
    candidates=filter(moduleexpr.args[end].args) do ex
        ex isa Expr&&ex.head==:macrocall&&ex.args[1]==Symbol("@kernel")&&ex.args[end].args[1].args[1]==nameof(kernel)
    end
    ex=only(candidates);original_args=copy(ex.args[end].args[1].args[2:end])
    old,new,generated=batch_kernel_ast(deepcopy(ex))
    fn=generated.args[end];sig=fn.args[1];body=fn.args[2]
    sig.args[1]=Symbol("policybatch_",new,"_",string(objectid(@__MODULE__);base=16))
    insertions=Expr[]
    for arg in original_args
        j=findfirst(==(arg),sig.args)
        if j===nothing
            # Cell kernels treat q/cq as common physical parameters. Policy
            # chemistry can contain private scalars (e.g. a critic gain).
            if arg in (:q,:cq)
                packed=Symbol("packed_",arg)
                for statement in body.args
                    if statement isa Expr&&statement.head==:(=)&&statement.args[1]===arg
                        statement.args[2]=:($(GlobalRef(@__MODULE__,:batch_world))($packed,batch_world_id))
                    end
                end
            end
            continue
        end
        packed=Symbol("policy_packed_",arg);sig.args[j]=packed
        push!(insertions,:($arg=$(GlobalRef(@__MODULE__,:batch_world))($packed,batch_world_id)))
    end
    indexpos=findfirst(x->x isa Expr&&x.head==:(=)&&x.args[2] isa Expr&&x.args[2].head==:macrocall&&x.args[2].args[1]==Symbol("@index"),body.args)
    splice!(body.args,indexpos+1:indexpos,insertions)
    Core.eval(scope,generated)
    result=Base.invokelatest(getproperty,scope,sig.args[1]);WORLD_SIDECAR_KERNELS[kernel]=result
    result
end
struct BatchedPolicyLaunch{B,R,K}
    backend::B
    worlds::Int
    seeds::R
    kernels::K
end
function batched_policy_launch!(launch,h,kernel,args...;ndrange)
    dims=ndrange isa Tuple ? ndrange : (ndrange,)
    launch.kernels[kernel](launch.backend,64)(args...,launch.seeds;ndrange=(dims...,launch.worlds))
end
function register_sidecar_batch_policy!(policy)
    map=IdDict()
    for (kernel,scope,path) in policy.kernels;map[kernel]=project_sidecar_kernel!(kernel,scope,path);end
    scope=policy.launch_module
    Core.eval(scope,:(policy_launch!(x::$(BatchedPolicyLaunch),h,kernel,args...;ndrange)=$(GlobalRef(@__MODULE__,:batched_policy_launch!))(x,h,kernel,args...;ndrange)))
    WORLD_SIDECAR_POLICIES[policy.name]=(policy=policy,kernels=map)
    policy
end
function sidecar_batch_policy(h)
    matches=[entry for entry in values(WORLD_SIDECAR_POLICIES) if h isa entry.policy.sidecar_type]
    length(matches)==1||throw(ArgumentError("register exactly one batch policy for this sidecar type"))
    only(matches)
end
function restore_lifetime_sidecar(p,x)
    q=x.options
    NativeLifetimeCredit.sidecar(p;dt=Float64(q.dt),decision_time=Float64(q.dt)*q.decision_stride,
        reward_time=Float64(q.dt)*q.reward_stride,rate=q.rate,sigma=q.sigma,kappa=q.kappa,seed=q.seed,
        learning=q.learning,critic=q.critic,critic_learning=q.critic_learning)
end
const LIFETIME_BATCH_POLICY=register_sidecar_batch_policy!(SidecarBatchPolicy(:lifetime,NativeLifetimeCredit.CreditHabitat,
    NativeLifetimeCredit,((NativeLifetimeCredit.choose!,NativeLifetimeCredit,joinpath(@__DIR__,"native_lifetime_credit.jl")),
    (NativeLifetimeCredit.reward!,NativeLifetimeCredit,joinpath(@__DIR__,"native_lifetime_credit.jl"))),
    LEARNING_MEMBER_OPTIONS,(),(),(),
    (h,physical,state,backend,options,cache,properties,launcher)->NativeLifetimeCredit.CreditHabitat(physical,state,backend,options,h.steps,launcher),
    restore_lifetime_sidecar,(lifetime=bytes2hex(sha256(read(joinpath(@__DIR__,"native_lifetime_credit.jl")))),)))
const LEARNING_BATCH_SOURCE=(adapter=bytes2hex(sha256(read(@__FILE__))),lifetime=LIFETIME_BATCH_POLICY.source,constructive=CONSTRUCTIVE_BATCH_SOURCE)
struct LearningWorldBatch{P,H,O}
    physical::P
    learner::H
    contract::O
end
function learning_world_batch(ps,learners;backend=nothing,identities=["world-$i" for i in eachindex(ps)])
    length(ps)==length(learners)>0||throw(DimensionMismatch("one paid learner per physical world required"))
    backend===nothing&&(backend=first(ps).cells.backend)
    length(identities)==length(ps)&&length(unique(identities))==length(ps)||throw(ArgumentError("unique experiment identities required"))
    reference=first(learners);entry=sidecar_batch_policy(reference);policy=entry.policy;common=reference.options
    common_properties=NamedTuple{policy.cohort_properties}(Tuple(getproperty(reference,k) for k in policy.cohort_properties))
    invariant(q)=NamedTuple{Tuple(k for k in keys(q) if !(k in policy.member_options))}(Tuple(q[k] for k in keys(q) if !(k in policy.member_options)))
    for (p,h) in zip(ps,learners)
        all(k->getproperty(h,k)==common_properties[k],policy.cohort_properties)||throw(ArgumentError("policy cohort clocks differ"))
        h.habitat===p.habitat||throw(ArgumentError("learner does not own this habitat"))
        sidecar_batch_policy(h).policy===policy&&h.steps==reference.steps&&invariant(h.options)==invariant(common)||throw(ArgumentError("batch requires common policy, physical/decision/reward/forcing clocks and action rows"))
    end
    physical=constructive_world_batch(ps;backend)
    upload(group)=backend==reference.backend ? group : map(a->batch_upload(backend,a),group)
    state=upload(pack_world_arrays(getproperty.(learners,:state)))
    cache=upload(pack_world_arrays([NamedTuple{policy.cache_names}(Tuple(getproperty(h,k) for k in policy.cache_names)) for h in learners]))
    function scalar(k,property=false)
        values=[property ? getproperty(h,k) : h.options[k] for h in learners]
        storage=KernelAbstractions.zeros(backend,eltype(values),1,length(values));copyto!(storage,reshape(values,1,:))
        WorldBatchScalar(WorldBatchArray(storage,(1,)))
    end
    opts=NamedTuple{policy.member_options}(Tuple(scalar(k) for k in policy.member_options))
    properties=NamedTuple{policy.member_properties}(Tuple(scalar(k,true) for k in policy.member_properties))
    launcher=BatchedPolicyLaunch(backend,length(ps),(cells=physical.cells.seeds,habitat=physical.cells.habitat_seeds),entry.kernels)
    mixed=merge(common,opts)
    h=policy.repack(reference,physical.habitat,state,backend,mixed,cache,merge(common_properties,properties),launcher)
    own=(physical=physical,habitat=h.habitat,state=h.state,member_options=opts,properties=properties,common_properties=common_properties,cache=cache,
        shapes=map(size,state),cache_shapes=map(size,cache),options=mixed,common_options=common,policy=policy,launcher=launcher,identities=Tuple(String.(identities)))
    LearningWorldBatch(physical,h,own)
end
function learning_batch_validate(worlds,b)
    p=b.physical;h=b.learner;o=b.contract
    p===o.physical&&h.habitat===p.habitat===o.habitat&&h.state===o.state&&h.options===o.options&&h.launcher===o.launcher||throw(ArgumentError("learning batch ownership changed"))
    map(size,h.state)==o.shapes&&map(size,o.cache)==o.cache_shapes||throw(ArgumentError("learning state/cache shapes changed"))
    all(k->getproperty(h,k)===o.cache[k],o.policy.cache_names)||throw(ArgumentError("policy cache ownership changed"))
    all(k->getproperty(h,k)==o.common_properties[k],o.policy.cohort_properties)||throw(ArgumentError("policy cohort clocks changed"))
    constructive_batch_validate_worlds(worlds,p);constructive_batch_validate_ownership(p)
    nothing
end
"""Advance a bounded block using the production queue, then synchronize once.

Structural ownership is checked at the block boundary. Every numerical kernel
still honors its own world's failure latch, and clocks count completed steps.
"""
function learning_world_batch_advance!(worlds,b,steps::Integer,dt;check_status=true,execution=:projected,lifecycle=:prefix)
    learning_batch_validate(worlds,b)
    p=b.physical;c=p.cells;h=b.learner
    delta=Float32(dt)
    steps>=0&&delta==h.options.dt&&0<delta<=Interior.stability_limit(first(worlds),p)||throw(ArgumentError("invalid paid learning batch block/timestep"))
    execution in (:staged,:projected)&&lifecycle in (:prefix,:dense)||throw(ArgumentError("unknown execution recipe"))
    before=copy(Array(c.environment.metadata.storage)[6,:])
    for _ in 1:steps
        c.attempted_steps+=1
        Interior.device_cellular_enqueue!(first(worlds),c,delta,execution,h,UInt32(c.attempted_steps),lifecycle)
    end
    KernelAbstractions.synchronize(c.backend)
    after=Array(c.environment.metadata.storage)[6,:]
    for i in eachindex(worlds),_ in 1:Int(after[i]-before[i]);worlds[i].time+=Float64(delta);end
    check_status&&Interior.device_cellular_status!(c)
    b
end
learning_world_batch_step!(worlds,b,dt;kwargs...)=learning_world_batch_advance!(worlds,b,1,dt;kwargs...)
# Download each packed parent once, then produce ordinary host snapshots. This
# avoids W full-parent downloads for W observations without changing any state.
function learning_world_batch_snapshot(b;check_status=true)
    p=b.physical;c=p.cells;h=p.habitat
    check_status&&Interior.device_cellular_status!(c)
    host(group)=map(a->Array(a.storage),group)
    state=host(c.state);environment=host(c.environment);arrays=host(h.arrays);feed=host(h.feedstock.arrays);learning=host(b.learner.state)
    opts=map(a->Array(a.values.storage),b.contract.member_options)
    properties=map(a->Array(a.values.storage),b.contract.properties);cache=host(b.contract.cache)
    seeds=Array(c.seeds);habitat_seeds=Array(c.habitat_seeds)
    slice(group,i)=map(a->copy(selectdim(a,ndims(a),i)),group)
    map(1:c.worlds) do i
        cells=(state=slice(state,i),environment=slice(environment,i),parameters=merge(c.parameters,(seed=seeds[i],)),geometry=c.geometry,attempted_steps=c.attempted_steps)
        physical=(;cells,arrays=slice(arrays,i),feedstock=slice(feed,i),parameters=merge(h.parameters,(seed=habitat_seeds[i],)),feed_parameters=h.feedstock.parameters)
        options=merge(b.contract.common_options,map(a->a[1,i],opts))
        (;physical,state=slice(learning,i),options,steps=b.learner.steps,cache=slice(cache,i),properties=merge(b.contract.common_properties,map(a->a[1,i],properties)))
    end
end
function learning_world_batch_checkpoint(worlds,b;members=collect(eachindex(worlds)))
    learning_batch_validate(worlds,b)
    indices=Int.(collect(members));length(unique(indices))==length(indices)&&!isempty(indices)&&all(i->1<=i<=length(worlds),indices)||throw(ArgumentError("invalid checkpoint members"))
    snapshots=learning_world_batch_snapshot(b;check_status=false)
    physical=map(indices) do i
        host=constructive_batch_host_carrier(snapshots[i].physical)
        w=World(worlds[i].grid,CPU(),NamedTuple(),NamedTuple(),worlds[i].time,worlds[i].seed)
        Interior.device_habitat_checkpoint(w,host)
    end
    physical_payload=(format="Interior.ExperimentalConstructiveWorldBatch",version=1,julia_version=VERSION,
        adapter_source=CONSTRUCTIVE_BATCH_SOURCE,source_member_count=length(worlds),member_indices=indices,members=physical)
    journal=map(indices) do i
        x=snapshots[i];metadata=x.physical.cells.environment.metadata
        (;identity=b.contract.identities[i],completed_steps=Int(metadata[6]),attempted_steps=x.physical.cells.attempted_steps,
            learning_steps=x.steps,time=worlds[i].time,seed=worlds[i].seed,status=only(x.physical.cells.environment.status))
    end
    (;format="Interior.ExperimentalLearningWorldBatch",version=2,julia_version=VERSION,source=LEARNING_BATCH_SOURCE,policy=b.contract.policy.name,policy_source=b.contract.policy.source,
        physical=physical_payload,learners=[(;state=snapshots[i].state,options=snapshots[i].options,steps=snapshots[i].steps,cache=snapshots[i].cache,properties=snapshots[i].properties) for i in indices],journal)
end

function learning_sidecar_restore(policy,physical,x)
    h=policy.restore(physical,x);q=x.options
    invariant(options)=NamedTuple{Tuple(k for k in keys(options) if !(k in policy.member_options))}(Tuple(options[k] for k in keys(options) if !(k in policy.member_options)))
    invariant(h.options)==invariant(q)||throw(ArgumentError("learner options differ from current clock/schema contract"))
    h.options=q
    keys(x.cache)==policy.cache_names&&keys(x.properties)==(policy.cohort_properties...,policy.member_properties...)||throw(ArgumentError("policy cache/property roles differ"))
    for key in policy.cache_names
        a=getproperty(h,key);v=x.cache[key]
        size(a)==size(v)&&eltype(a)===eltype(v)&&all(isfinite,v)||throw(ArgumentError("invalid policy cache $key"))
        copyto!(a,v)
    end
    keys(h.state)==keys(x.state)||throw(ArgumentError("learner state roles differ"))
    for k in keys(h.state)
        a=h.state[k];v=x.state[k]
        size(a)==size(v)&&eltype(a)===eltype(v)&&all(isfinite,v)||throw(ArgumentError("invalid learner state $k"))
        copyto!(a,v)
    end
    all(x->x in (0,1),h.state.ready)&&all(>=(0),h.state.paid)&&maximum(abs,h.state.learned)<=.3f0&&maximum(abs,h.state.critic)<=.3f0||throw(ArgumentError("invalid paid learner state"))
    x.steps>=0||throw(ArgumentError("invalid learner clock"));h.steps=x.steps
    h
end
function learning_world_batch_restore(worlds,x)
    x.format=="Interior.ExperimentalLearningWorldBatch"&&x.version==2&&x.julia_version==VERSION&&x.source==LEARNING_BATCH_SOURCE||throw(ArgumentError("learning batch source/schema differs"))
    length(worlds)==length(x.learners)==length(x.journal)>0||throw(DimensionMismatch("learning checkpoint member count differs"))
    haskey(WORLD_SIDECAR_POLICIES,x.policy)||throw(ArgumentError("register the checkpoint policy before restoring"))
    policy=WORLD_SIDECAR_POLICIES[x.policy].policy
    policy.source==x.policy_source||throw(ArgumentError("policy source differs"))
    xp=x.physical
    xp.format=="Interior.ExperimentalConstructiveWorldBatch"&&xp.version==1&&xp.julia_version==VERSION&&xp.adapter_source==CONSTRUCTIVE_BATCH_SOURCE||throw(ArgumentError("physical batch source/schema differs"))
    length(xp.members)==length(worlds)==length(xp.member_indices)&&length(unique(xp.member_indices))==length(worlds)&&all(i->1<=i<=xp.source_member_count,xp.member_indices)||throw(ArgumentError("invalid physical member journal"))
    identities=getproperty.(x.journal,:identity)
    length(unique(identities))==length(identities)||throw(ArgumentError("duplicate experiment identity"))
    # Validate sidecars and per-member journals before mutating caller clocks.
    hosts=map(eachindex(worlds)) do i
        cp=x.physical.members[i]
        w=World(worlds[i].grid,CPU(),NamedTuple(),NamedTuple(),0.,worlds[i].seed)
        Interior.device_habitat_restore(w,cp)
    end
    learners=map(eachindex(hosts)) do i
        h=learning_sidecar_restore(policy,hosts[i],x.learners[i]);j=x.journal[i];cp=x.physical.members[i]
        metadata=hosts[i].cells.environment.metadata
        j.status==0&&j.completed_steps==metadata[6]&&j.attempted_steps==hosts[i].cells.attempted_steps&&
            j.learning_steps==h.steps&&j.time==cp.cells.world_time&&j.seed==cp.cells.world_seed||throw(ArgumentError("member journal differs from complete physical/learning state"))
        h
    end
    batch=learning_world_batch(hosts,learners;backend=first(worlds).backend,identities)
    constructive_batch_validate_worlds(worlds,batch.physical)
    for i in eachindex(worlds);worlds[i].time=x.journal[i].time;worlds[i].seed=x.journal[i].seed;end
    batch
end
learning_world_batch_observation(b)=map(learning_snapshot_observation,learning_world_batch_snapshot(b))

function learning_benchmark_members(profile,count;backend=CPU(),capacity=64,population=4)
    template=TOML.parsefile(profile)
    template["grid"]["shape"]=[8,8,8];template["grid"]["extent"]=[1.,1.,1.]
    template["cellular"]["founders"]=population;template["cellular"]["capacity"]=capacity
    template["cellular"]["growth_rate"]=0.;template["cellular"]["mutation"]=0.
    template["cellular"]["adaptation_strength"]=.3
    for section in ("resources","feedstock"),key in ("source_rates","source_budgets","sink_rates","sink_budgets")
        template[section][key].=0.
    end
    template["native"]["dt"]=.02;template["native"]["energy_type"]="Float64";template["native"]["inventory_type"]="Float64"
    # Configuration only: the shared production initializer owns all arrays,
    # inherited priors, physical budgets and native model construction.
    members=map(1:count) do i
        config=deepcopy(template);config["seed"]=211+i
        config["resources"]["initial"]=isodd(i) ? [.1,.5,0.] : [.5,.1,0.]
        w,p,_=NativeLifetimeCredit.NativeActionCredit.native_constructed_garden(backend,config)
        h=NativeLifetimeCredit.sidecar(p;dt=.02,decision_time=.04,reward_time=.08,seed=6200+i,
            rate=isodd(i) ? .003 : -.003,critic=iseven(i))
        (w,p,h)
    end
    getindex.(members,1),getindex.(members,2),getindex.(members,3)
end
function learning_serial_advance!(worlds,physical,learners,steps,dt)
    # Enqueue every independent owner's complete trajectory before synchronizing.
    # This baseline does not pay a device barrier on each step or each world.
    for (w,p,h) in zip(worlds,physical,learners),_ in 1:steps
        Interior.device_cellular_step!(w,p.cells,dt;anabolism=h,check_status=false)
    end
    KernelAbstractions.synchronize(first(physical).cells.backend)
    for (i,p) in enumerate(physical)
        code=only(Array(p.cells.environment.status));code==0||error("serial member $i failed with status $code")
    end
    for w in worlds,_ in 1:steps;w.time+=Float64(Float32(dt));end
    nothing
end
function learning_serial_snapshot(physical,learners)
    [(physical=Interior.device_habitat_snapshot(p),state=map(Array,h.state),options=h.options,steps=h.steps) for (p,h) in zip(physical,learners)]
end
function learning_snapshot_observation(x)
    a=Interior.device_habitat_audit(x.physical);s=x.physical.cells.state;active=findall(==(1),s.active)
    (;ids=s.ids[active],positions=s.position[active,:],mass=s.mass[active],energy=s.energy[active],
        learning_work=x.state.paid[active],mark_max=maximum(abs,x.state.learned),critic_max=maximum(abs,x.state.critic),
        last_return=x.state.last_reward[active],baseline=x.state.baseline[active],audit=a,steps=x.steps)
end
function learning_array_error(x,y)
    size(x)==size(y)||return Inf
    maximum((isequal(a,b) ? 0. : abs(Float64(a)-Float64(b)) for (a,b) in zip(x,y));init=0.)
end
function learning_batch_parity(serial,batched;tolerance=0.)
    length(serial)==length(batched)||error("member count differs")
    maxerror=0.
    for (a,b) in zip(serial,batched)
        a.options==b.options&&a.steps==b.steps||error("member learning clocks/options differ")
        for key in (:active,:ids,:parents,:lineages,:generations,:bonds,:bond_ids)
            a.physical.cells.state[key]==b.physical.cells.state[key]||error("member identity/topology differs: $key")
        end
        for group in (:state,:environment)
            left=a.physical.cells[group];right=b.physical.cells[group]
            for key in keys(left)
                eltype(left[key])<:AbstractFloat||continue
                err=learning_array_error(left[key],right[key])
                maxerror=max(maxerror,err)
            end
        end
        for key in keys(a.state)
            x=a.state[key];y=b.state[key]
            if eltype(x)<:AbstractFloat
                maxerror=max(maxerror,learning_array_error(x,y))
            else
                x==y||error("member private learner identity differs: $key")
            end
        end
        for group in (:arrays,:feedstock)
            for key in keys(a.physical[group])
                x=a.physical[group][key];y=b.physical[group][key]
                key===:enzyme_bits&&(x=reinterpret(Float32,x);y=reinterpret(Float32,y))
                maxerror=max(maxerror,learning_array_error(x,y))
            end
        end
    end
    maxerror<=tolerance||error("complete member parity failed: $maxerror")
    maxerror
end

"""Bounded end-to-end ensemble benchmark; reports aggregate and per-world costs.

The assay uses 8³ closed, nondividing production worlds with material receptors,
paid actors/critics, and explicit faster learning clocks (.04/.08 physical time).
No claim about the utility of learned behavior follows from this throughput gate.
"""
function benchmark_learning_world_batches(profile,output;backend=CPU(),counts=(8,32),steps=16,repeats=3,tolerance=0.,total_capacity=4096,total_population=32)
    warm_started=time_ns()
    # Exercise compilation, checkpoint, restore and observation outside samples.
    ww,pp,hh=learning_benchmark_members(profile,2)
    warm=learning_world_batch(pp,hh;backend)
    warmworlds=[World(w.grid,backend,NamedTuple(),NamedTuple(),0.,w.seed) for w in ww]
    learning_world_batch_advance!(warmworlds,warm,8,.02)
    cp=learning_world_batch_checkpoint(warmworlds,warm)
    rw=[World(w.grid,backend,NamedTuple(),NamedTuple(),0.,w.seed) for w in ww]
    learning_world_batch_restore(rw,cp);learning_world_batch_observation(warm)
    # Warm standalone device kernels too, using actual paid updates.
    sw,sp,sh=learning_benchmark_members(profile,1;backend)
    learning_serial_advance!(sw,sp,sh,8,.02);learning_serial_snapshot(sp,sh)
    global_warmup_seconds=(time_ns()-warm_started)/1e9
    results=Dict{String,Any}()
    for count in counts
        memory=nothing
        total_capacity%count==0&&total_population%count==0||error("total capacity/population must divide member counts")
        capacity=total_capacity÷count;population=total_population÷count
        # KA launch geometry can specialize by capacity/world partition; warm
        # each exact case rather than assuming one small cohort covers it.
        case_warm_started=time_ns()
        cw,cp,ch=learning_benchmark_members(profile,count;capacity,population)
        casebatch=learning_world_batch(cp,ch;backend)
        caseworlds=[World(w.grid,backend,NamedTuple(),NamedTuple(),0.,w.seed) for w in cw]
        csw,csp,csh=learning_benchmark_members(profile,count;backend,capacity,population)
        learning_serial_advance!(csw,csp,csh,8,.02)
        learning_world_batch_advance!(caseworlds,casebatch,8,.02)
        learning_batch_parity(learning_serial_snapshot(csp,csh),learning_world_batch_snapshot(casebatch);tolerance)
        case_warmup_seconds=(time_ns()-case_warm_started)/1e9
        memory=learning_world_batch_memory(casebatch)
        println("warmed worlds=",count," capacity=",capacity," population=",population," owned_bytes=",sum(values(memory)));flush(stdout)
        samples=NamedTuple[]
        for repeat in 1:repeats
            arms=Dict{Symbol,Any}()
            for arm in (isodd(repeat) ? (:serial,:batch) : (:batch,:serial))
                GC.gc();start=time_ns()
                if arm===:serial
                    worlds,physical,learners=learning_benchmark_members(profile,count;backend,capacity,population)
                    setup=(time_ns()-start)/1e9;GC.gc();start=time_ns()
                    learning_serial_advance!(worlds,physical,learners,steps,.02)
                    stepping=(time_ns()-start)/1e9;start=time_ns()
                    snapshot=learning_serial_snapshot(physical,learners)
                    observations=map(learning_snapshot_observation,snapshot)
                    observing=(time_ns()-start)/1e9
                else
                    cpuworlds,physical,learners=learning_benchmark_members(profile,count;capacity,population)
                    batch=learning_world_batch(physical,learners;backend)
                    memory=learning_world_batch_memory(batch)
                    worlds=[World(w.grid,backend,NamedTuple(),NamedTuple(),0.,w.seed) for w in cpuworlds]
                    setup=(time_ns()-start)/1e9;GC.gc();start=time_ns()
                    learning_world_batch_advance!(worlds,batch,steps,.02)
                    stepping=(time_ns()-start)/1e9;start=time_ns()
                    snapshot=learning_world_batch_snapshot(batch)
                    observations=map(learning_snapshot_observation,snapshot)
                    observing=(time_ns()-start)/1e9
                end
                arms[arm]=(;setup,stepping,observing,total=setup+stepping+observing,snapshot,observations)
            end
            err=learning_batch_parity(arms[:serial].snapshot,arms[:batch].snapshot;tolerance)
            compact(x)=(;setup=x.setup,stepping=x.stepping,observing=x.observing,total=x.total,
                paid=sum(sum(o.learning_work) for o in x.observations),finite=all(o->o.audit.finite,x.observations))
            push!(samples,(serial=compact(arms[:serial]),batch=compact(arms[:batch]),maxerror=err))
        end
        med(xs)=sort(xs)[cld(length(xs),2)]
        serialstep=med([x.serial.stepping for x in samples]);batchstep=med([x.batch.stepping for x in samples])
        serialtotal=med([x.serial.total for x in samples]);batchtotal=med([x.batch.total for x in samples])
        setup_difference=med([x.batch.setup for x in samples])-med([x.serial.setup for x in samples])
        saved_per_shared_step=(serialstep-batchstep)/steps
        setup_break_even_steps=saved_per_shared_step>0 ? max(0,ceil(Int,setup_difference/saved_per_shared_step)) : nothing
        results[string(count)]=(;worlds=count,case_warmup_seconds,setup_difference,setup_break_even_steps,capacity_per_world=capacity,population_per_world=population,total_capacity,total_population,total_voxels=512count,owned_bytes=memory,total_owned_bytes=sum(values(memory)),steps,repeats,samples,aggregate_step_speedup=serialstep/batchstep,
            end_to_end_speedup=serialtotal/batchtotal,batch_world_steps_per_second=count*steps/batchstep,
            batch_wall_ms_per_shared_step=1000batchstep/steps,batch_amortized_ms_per_world_step=1000batchstep/(steps*count))
    end
    report=(;schema="Interior.LearningWorldBatchBenchmark.v1",julia=string(VERSION),threads=Threads.nthreads(),backend=string(typeof(backend)),
        source=LEARNING_BATCH_SOURCE,global_warmup_seconds,warmup_scope="initialization and compile/gates of both serial and batch modes; not an isolated incremental batch compile penalty",serial_sync_policy="all independent owners enqueue before one backend barrier; one status read per standalone owner",profile_sha256=bytes2hex(sha256(read(profile))),counts,steps,repeats,tolerance,
        scope="independent coupled 8^3 worlds, production physics + paid lifetime actor/critic, .02 physics/.04 decision/.08 reward clocks; aggregate throughput, not faster single-world dynamics",results)
    open(io->Interior.JSON3.pretty(io,report),output,"w")
    report
end
function learning_world_batch_memory(b)
    c=b.physical.cells;h=b.physical.habitat;work=Interior.device_resource_workspace(c)
    groups=(cells=values(c.state),environment=values(c.environment),habitat=values(h.arrays),feedstock=values(h.feedstock.arrays),
        learning=values(b.learner.state),options=values(b.contract.member_options),properties=values(b.contract.properties),cache=values(b.contract.cache),
        workspace=(work.grants,(a for level in work.levels for a in (level.high,level.low))...,Interior.device_lifecycle_workspace(c).highwater,c.seeds,c.habitat_seeds))
    seen=Set{UInt}()
    map(groups) do arrays
        bytes=0
        for a in arrays
            parent=a isa WorldBatchScalar ? a.values.storage : a isa WorldBatchArray ? a.storage : a
            key=objectid(parent);key in seen&&continue;push!(seen,key)
            bytes+=length(parent)*sizeof(eltype(parent))
        end
        bytes
    end
end

construction_sidecar_keywords(q)=(dt=parse(Float64,string(q.dt)),rate=q.rate,sigma=q.sigma,seed=q.seed,
    learning=q.learning,critic=q.critic,enabled=q.enabled,phase_offset=q.phase_offset)
function register_construction_world_batches!(C;name=:construction)
    path=joinpath(@__DIR__,"native_construction_credit.jl");lifetime=joinpath(@__DIR__,"native_lifetime_credit.jl")
    kernels=((C.choose!,C,path),(C.actuate!,C,path),(C.feed_phase!,C,path),(C.L.reward!,C.L,lifetime))
    policy=SidecarBatchPolicy(name,C.ConstructionHabitat,C.L,kernels,(LEARNING_MEMBER_OPTIONS...,:enabled),(),(),(:original_masks,),
        (h,physical,state,backend,options,cache,properties,launcher)->C.ConstructionHabitat(physical,state,backend,options,cache.original_masks,h.steps,launcher),
        (p,x)->C.sidecar(p;construction_sidecar_keywords(x.options)...),
        (construction=bytes2hex(sha256(read(path))),lifetime=bytes2hex(sha256(read(lifetime)))))
    register_sidecar_batch_policy!(policy)
end
function register_joint_construction_world_batches!(J;name=:joint_construction)
    C=J.C;path=joinpath(@__DIR__,"native_construction_credit.jl");joint=joinpath(@__DIR__,"native_joint_construction_credit.jl");lifetime=joinpath(@__DIR__,"native_lifetime_credit.jl")
    kernels=((J.choose_joint!,J,joint),(J.reserve!,J,joint),(C.actuate!,C,path),(C.feed_phase!,C,path),(J.L.reward!,J.L,lifetime))
    policy=SidecarBatchPolicy(name,J.JointHabitat,J.L,kernels,(LEARNING_MEMBER_OPTIONS...,:enabled),(:memory,),(),(:original_masks,),
        (h,physical,state,backend,options,cache,properties,launcher)->J.JointHabitat(C.ConstructionHabitat(physical,state,backend,options,cache.original_masks,h.steps,launcher),properties.memory),
        (p,x)->J.sidecar(p;memory=x.properties.memory,construction_sidecar_keywords(x.options)...),
        (joint=bytes2hex(sha256(read(joint))),construction=bytes2hex(sha256(read(path))),lifetime=bytes2hex(sha256(read(lifetime)))))
    register_sidecar_batch_policy!(policy)
end
