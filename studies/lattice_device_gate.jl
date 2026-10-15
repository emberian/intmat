# Baseline finite-genotype lattice law only. GPU execution requires an allocated window.
using Interior,SHA,Serialization
using Interior:JSON3,KernelAbstractions
include("lattice_ecology.jl")
using .LatticeEcologyStudy
const LD=LatticeEcologyStudy
const EXPECTED_LATTICE_SHA="e6d6cab2b0624965a5e7ec66b7f9f9ce44a01f9ded5653db660de61368cde9b5"
@assert LD.MODEL_SHA==EXPECTED_LATTICE_SHA "gate is pinned to the reviewed baseline law"
const LATTICE_DEVICE_PROTOCOL=(float_atol=1e-4,float_rtol=2e-4,ledger_relative=2e-5,
    capacity_atol=2e-6,hidden_bound=1.00001,dt=.1,parity_steps=200,continuation_steps=17,
    shapes=((8,8,1),(4,4,4)),benchmark_shapes=((32,32,32),(64,64,64)),
    warmup_steps=128,benchmark_batch_steps=128,benchmark_min_seconds=3.,benchmark_max_batches=128,
    benchmark_repeats=3,minimum_benchmark_occupancy=.5,
    benchmark_initial_occupancy=.65,benchmark_initial_mass=.7,benchmark_precursor=2.,benchmark_reserve=1024.,
    benchmark_construction=.006,benchmark_matrix_decay=.02)
const PERSISTENT=(:environment,:cells,:hidden,:actions,:heat,:work,:genotype,:lineage,:generation,:events,:ids,:parents,:weights)
function gate_backend(name)
    name=="cpu"&&return CPU()
    name=="amd"||error("backend must be cpu or amd")
    Core.eval(@__MODULE__,:(using AMDGPU))
    Core.eval(@__MODULE__,:(AMDGPU.functional()))||error("native AMD runtime unavailable")
    Core.eval(@__MODULE__,:(AMDGPU.allowscalar(false)))
    Core.eval(@__MODULE__,:(AMDGPU.ROCBackend()))
end
function valid_lattice(w,p,initial)
    s=lattice_statistics(w,p);q=LATTICE_DEVICE_PROTOCOL
    s.finite&&s.minimum_physical>=0&&s.capacity_excess<=q.capacity_atol&&s.hidden_maxabs<=q.hidden_bound&&
        abs(s.material-initial.material)<=q.ledger_relative*max(1,initial.material)&&
        abs(s.energy_accounted-initial.energy_accounted)<=q.ledger_relative*max(1,initial.energy_accounted)
end
function compare_lattice(w,p,v,q;exact=false,scratch=true)
    errors=Dict{String,Any}();passed=w.time==v.time&&p.tick==q.tick
    for name in (scratch ? (PERSISTENT...,:scratch) : PERSISTENT)
        a=Array(getfield(p,name));b=Array(getfield(q,name))
        good=eltype(a)<:Integer||exact ? a==b : all(abs.(a.-b).<=LATTICE_DEVICE_PROTOCOL.float_atol.+LATTICE_DEVICE_PROTOCOL.float_rtol.*abs.(a))
        passed&=good
        errors[string(name)]=(passed=good,maxabs=maximum(abs.(Float64.(a).-Float64.(b))),different=count(a.!=b))
    end
    (passed=passed,fields=errors)
end
function expected_rejection(f)
    try;f();false;catch e;e isa ArgumentError||e isa DimensionMismatch;end
end
function forced_contact(backend,kind)
    w=lattice_world(Grid(2,2,2;extent=(2,2,2));seed=11,backend)
    p=LatticeEcology(w;loci=4,occupancy=0,founder_genotypes=1,
        parameters=lattice_parameters(move_rate=1,division_rate=1,mutation_probability=1))
    # Stage on host; scalar device indexing is forbidden.
    g=zeros(Int32,8);g[1]=1;ids=zeros(UInt64,8);ids[1]=1;c=zeros(Float32,8,3);c[1,:].=(1,1,7)
    a=zeros(Float32,8,3);a[1,kind==:birth ? 3 : 2]=1
    h=zeros(Float32,8,4);h[1,:].=.4
    for (dest,src) in ((p.genotype,g),(p.lineage,g),(p.ids,ids),(p.cells,c),(p.actions,a),(p.hidden,h));copyto!(dest,src);end
    initial=lattice_statistics(w,p)
    LD.matching_kernel!(backend,128)(p.cells,p.hidden,p.actions,p.environment,p.genotype,p.lineage,p.generation,p.ids,p.parents,p.events,p.heat,p.work,p.dims,p.parameters,1,1,p.seed,1,p.loci,1f0;ndrange=8)
    KernelAbstractions.synchronize(backend)
    final=lattice_statistics(w,p)
    event=kind==:birth ? final.births==1&&final.mutations==1&&final.work[5]>0 : final.hops==1&&final.births==0&&final.work[4]>0
    @assert event&&valid_lattice(w,p,initial)
    w,p
end
function gate_fixture(backend,dims,output)
    mkpath(output);w=lattice_world(Grid(dims...;extent=Float32.(dims));seed=29)
    p=LatticeEcology(w;loci=6,founder_genotypes=4,occupancy=.25)
    initial=lattice_statistics(w,p);initialpath=joinpath(output,"initial.jls");lattice_checkpoint(w,p,initialpath)
    v,q=restore_lattice(initialpath;backend);@assert compare_lattice(w,p,v,q;exact=true).passed
    comparisons=Any[]
    for step in 1:LATTICE_DEVICE_PROTOCOL.parity_steps
        step!(w,p,.1);step!(v,q,.1)
        if step%10==0
            KernelAbstractions.synchronize(backend);comparison=compare_lattice(w,p,v,q)
            @assert comparison.passed&&valid_lattice(w,p,initial)&&valid_lattice(v,q,initial)
            push!(comparisons,(step=step,comparison=comparison))
        end
    end
    stats=lattice_statistics(v,q)
    @assert stats.births>0&&stats.mutations>0&&stats.hops>0&&stats.matrix>0&&stats.work[3]>0
    record_state(w,p,joinpath(output,"cpu-native");step=200);record_state(v,q,joinpath(output,"target-native");step=200)
    checkpoint=joinpath(output,"target.jls");lattice_checkpoint(v,q,checkpoint)
    x,r=restore_lattice(checkpoint;backend);restored=compare_lattice(v,q,x,r;exact=true,scratch=false);@assert restored.passed
    for _ in 1:LATTICE_DEVICE_PROTOCOL.continuation_steps;step!(v,q,.1);step!(x,r,.1);end
    KernelAbstractions.synchronize(backend);continued=compare_lattice(v,q,x,r;exact=true);@assert continued.passed
    saved=q.scratch;q.scratch=q.environment
    alias=expected_rejection(()->Interior.validate(v,q));q.scratch=saved;@assert alias
    before=(time=v.time,tick=q.tick,cells=Array(q.cells))
    bad_dt=expected_rejection(()->step!(v,q,-.1));@assert bad_dt&&v.time==before.time&&q.tick==before.tick&&Array(q.cells)==before.cells
    clobber=expected_rejection(()->lattice_checkpoint(v,q,checkpoint));@assert clobber
    c=deserialize(checkpoint)
    sourcebad=joinpath(output,"bad-source.jls");serialize(sourcebad,merge(c,(source_sha="invalid",)))
    source_rejected=expected_rejection(()->restore_lattice(sourcebad;backend));@assert source_rejected
    sharedbad=joinpath(output,"bad-shared.jls");serialize(sharedbad,merge(c,(shared_sources=Dict("invalid"=>"invalid"),)))
    shared_rejected=expected_rejection(()->restore_lattice(sharedbad;backend));@assert shared_rejected
    (dimensions=dims,passed=true,initial=initial,final=stats,comparisons=comparisons,
        exact_restore=restored,exact_continuation=continued,alias_rejected=alias,bad_dt_rejected=bad_dt,
        checkpoint_clobber_rejected=clobber,source_rejected=source_rejected,shared_source_rejected=shared_rejected)
end
function sustained_benchmark(backend,dims)
    w=lattice_world(Grid(dims...;extent=Float32.(dims));seed=29,backend)
    p=LatticeEcology(w;loci=12,founder_genotypes=8,occupancy=.65,initial_mass=.7,precursor=2,reserve=1024,
        parameters=lattice_parameters(construction=.006,matrix_decay=.02))
    for _ in 1:LATTICE_DEVICE_PROTOCOL.warmup_steps;step!(w,p,.1);end
    KernelAbstractions.synchronize(backend);initial=lattice_statistics(w,p)
    @assert initial.population>=LATTICE_DEVICE_PROTOCOL.minimum_benchmark_occupancy*prod(dims)
    warmed=map(name->Array(getfield(p,name)),PERSISTENT);warmed_time=w.time;warmed_tick=p.tick
    samples=Any[]
    for repeat in 1:LATTICE_DEVICE_PROTOCOL.benchmark_repeats
        # Each sustained repeat begins with exactly the same active warmed population.
        for (name,data) in zip(PERSISTENT,warmed);copyto!(getfield(p,name),data);end
        w.time=warmed_time;p.tick=warmed_tick;KernelAbstractions.synchronize(backend)
        seconds=0.;steps=0;batches=0
        while seconds<LATTICE_DEVICE_PROTOCOL.benchmark_min_seconds&&batches<LATTICE_DEVICE_PROTOCOL.benchmark_max_batches
            KernelAbstractions.synchronize(backend)
            seconds+=@elapsed begin
                for _ in 1:LATTICE_DEVICE_PROTOCOL.benchmark_batch_steps;step!(w,p,.1);end
                KernelAbstractions.synchronize(backend)
            end
            steps+=LATTICE_DEVICE_PROTOCOL.benchmark_batch_steps;batches+=1
        end
        stats=lattice_statistics(w,p)
        active=stats.population>=LATTICE_DEVICE_PROTOCOL.minimum_benchmark_occupancy*prod(dims)
        @assert valid_lattice(w,p,initial)&&active
        push!(samples,(repeat=repeat,steps=steps,seconds=seconds,batches=batches,
            reached_minimum_duration=seconds>=LATTICE_DEVICE_PROTOCOL.benchmark_min_seconds,
            voxel_steps_per_second=prod(dims)*steps/seconds,active_dense=active,statistics=stats))
    end
    (dimensions=dims,initial=initial,samples=samples,resident_bytes=initial.resident_bytes,
        sustained_duration_passed=all(s->s.reached_minimum_duration,samples))
end
function lattice_device_main(kind,output;benchmark=false)
    benchmark&&kind!="amd"&&error("volume benchmarks require an allocated AMD window")
    Threads.nthreads()<=2||error("at most two host threads");output=abspath(output)
    ispath(output)&&error("fresh output directory required");mkpath(output)
    cp(@__FILE__,joinpath(output,"gate-source.jl"));write(joinpath(output,"model-source.jl"),LD.MODEL_SOURCE)
    # Write protocol before loading/launching an accelerator or observing outcomes.
    provenance=(model_sha256=LD.MODEL_SHA,gate_sha256=bytes2hex(sha256(read(@__FILE__))),julia=string(VERSION),
        shared_sources=Interior.arena_sources(),protocol=LATTICE_DEVICE_PROTOCOL,backend_requested=kind)
    open(io->JSON3.pretty(io,provenance),joinpath(output,"protocol.json"),"w")
    backend=gate_backend(kind)
    Base.invokelatest(lattice_device_runs,backend,kind,output,benchmark,provenance)
end
function lattice_device_runs(backend,kind,output,benchmark,provenance)
    contacts=Any[]
    for event in (:birth,:hop)
        w,p=forced_contact(CPU(),event);v,q=forced_contact(backend,event)
        comparison=compare_lattice(w,p,v,q);@assert comparison.passed
        push!(contacts,(kind=string(event),comparison=comparison,statistics=lattice_statistics(v,q)))
    end
    fixtures=Any[]
    for dims in LATTICE_DEVICE_PROTOCOL.shapes
        push!(fixtures,gate_fixture(backend,dims,joinpath(output,join(dims,"x"))))
    end
    result=(passed=true,backend=kind,provenance=provenance,contacts=contacts,fixtures=fixtures)
    open(io->JSON3.pretty(io,result),joinpath(output,"gate.json"),"w")
    if benchmark
        kind=="amd"||error("sustained volume benchmark requires allocated AMD window")
        rows=Any[]
        for dims in LATTICE_DEVICE_PROTOCOL.benchmark_shapes
            push!(rows,sustained_benchmark(backend,dims))
            open(io->JSON3.pretty(io,(passed=all(r->r.sustained_duration_passed,rows),provenance=provenance,rows=rows)),joinpath(output,"benchmark.json"),"w")
        end
    end
    println("READY ",output," backend=",kind," benchmark=",benchmark)
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    length(ARGS) in (2,3)||error("cpu|amd fresh-output [--benchmark]")
    length(ARGS)==3&&ARGS[3]!="--benchmark"&&error("unknown option")
    lattice_device_main(ARGS[1],ARGS[2];benchmark=length(ARGS)==3&&ARGS[3]=="--benchmark")
end
