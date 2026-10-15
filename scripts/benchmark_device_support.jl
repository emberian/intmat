#!/usr/bin/env julia
# One sequential GPU process: direct initialization, complete laws, alternating pairs.
include("device_cellular_setup.jl")
include("device_habitat_gate.jl")
using Statistics
using Interior.KernelAbstractions: @kernel, @index
steps=parse(Int,get(ENV,"DEVICE_SUPPORT_STEPS","30"))
pairs=parse(Int,get(ENV,"DEVICE_SUPPORT_PAIRS","3"))
sizes=parse.(Int,split(get(ENV,"DEVICE_SUPPORT_SIZES","64,128,256"),','))
founders=parse(Int,get(ENV,"DEVICE_SUPPORT_FOUNDERS","512"))
capacity=parse(Int,get(ENV,"DEVICE_CAPACITY",string(4founders)))
all(>(0),(steps,pairs,founders)) && all(>(0),sizes) && capacity>=2founders || error("invalid work/capacity")
out=isempty(ARGS) ? "outputs/device-support-pairs" : ARGS[1]
ispath(out) && error("fresh output required")
mkpath(dirname(abspath(out)));mkdir(out)
@kernel function support_benchmark_prior!(s,h)
    i=@index(Global,Linear)
    if s.active[i]==1
        s.weights[i,end,h+8]=2f0;s.weights[i,end,h+9]=2f0;s.weights[i,end,h+10]=2f0
        i==1 && (s.mass[i]=.1f0)
    end
end
function support_population(nvox,which_backend,strategy,n,C)
    grid=Grid(nvox,nvox,nvox;extent=(4,4,4))
    w=World(grid,which_backend,NamedTuple(),NamedTuple(),0.,17)
    side=ceil(Int,cbrt(n));x=zeros(3,n)
    for (i,xyz) in enumerate(Iterators.product(0:side-1,0:side-1,0:side-1))
        i>n && break
        x[:,i].=4 .* (collect(xyz).+.5)./side
    end
    p=Interior.device_cellular_initialize(w,x;capacity=C,seed=17,support_strategy=strategy,
        exchange_radius=.18,density=128.,biomass=1.85,energy=1.,mutation=.01,
        initial=(.3,.3,0.),diffusion=(.001,.001,.001),mechanical_substeps=16)
    support_benchmark_prior!(which_backend,128)(p.state,p.geometry.hidden;ndrange=C)
    # Initial biological edit is part of this declared fixture, not an integrator stage.
    Interior.device_cellular_launch!(p,Interior.dc_clear_heads!,p.environment;ndrange=length(p.environment.heads))
    Interior.device_cellular_launch!(p,Interior.dc_geometry!,p.state,p.environment,p.parameters,p.geometry;ndrange=C)
    Interior.device_cellular_launch!(p,Interior.dc_clear_occupancy!,p.environment;ndrange=length(p.environment.occupancy))
    Interior.device_cellular_launch!(p,Interior.dc_porosity_scatter!,p.state,p.environment,p.geometry;ndrange=C)
    Interior.device_cellular_launch!(p,Interior.dc_porosity_commit!,p.environment,p.parameters;ndrange=length(p.environment.occupancy))
    Interior.device_cellular_status!(p)
    w,p
end
function support_memory(p)
    bytes(a)=sizeof(eltype(a))*length(a)
    support_keys=filter(k->k in (:sites,:support,:support_count,:axis_indices,:axis_distance2,:support_normalization,:support_fallback),keys(p.state))
    (support_bytes=sum(k->bytes(getproperty(p.state,k)),support_keys),
        state_bytes=sum(bytes,values(p.state)),environment_bytes=sum(bytes,values(p.environment)))
end
function support_errors(a,b)
    (identities_equal=a.state.active==b.state.active && a.state.ids==b.state.ids,
        position=maximum(abs.(a.state.position.-b.state.position)),
        hidden=maximum(abs.(a.state.hidden.-b.state.hidden)),
        metabolites=maximum(abs.(a.state.metabolites.-b.state.metabolites)),
        fields=maximum(abs.(a.environment.fields.-b.environment.fields)))
end
function support_gate()
    states=Dict{Symbol,Any}();errors=Any[]
    for strategy in (:materialized,:axis)
        wc,pc=support_population(16,CPU(),strategy,8,32)
        wg,pg=support_population(16,backend,strategy,8,32)
        for _ in 1:20
            Interior.device_cellular_step!(wc,pc,.001)
            Interior.device_cellular_step!(wg,pg,.001)
        end
        a=Interior.device_cellular_snapshot(pc);b=Interior.device_cellular_snapshot(pg)
        err=support_errors(a,b)
        err.identities_equal && all(<=(1e-4),(err.position,err.hidden,err.metabolites,err.fields)) || error("CPU/AMD parity failed")
        b.environment.metadata[3]>0 && b.environment.metadata[4]>0 || error("fixture missed births/deaths")
        any(>(0),b.state.bonds) || error("fixture missed contacts")
        push!(errors,(strategy=strategy,errors=err,stats=Interior.device_cellular_statistics(pg)))
        states[strategy]=b
        if strategy===:axis
            cp=Interior.device_cellular_checkpoint(wg,pg)
            restored=Interior.device_cellular_restore(wg,cp)
            Interior.device_cellular_snapshot(restored)==cp.snapshot || error("checkpoint bits differ")
            for _ in 1:20
                Interior.device_cellular_step!(wg,pg,.001)
                Interior.device_cellular_step!(wg,restored,.001)
            end
            continuation=support_errors(Interior.device_cellular_snapshot(pg),Interior.device_cellular_snapshot(restored))
            continuation.identities_equal && all(<=(1e-4),(continuation.position,continuation.hidden,continuation.metabolites,continuation.fields)) || error("continuation failed")
            push!(errors,(strategy=:axis_restart,errors=continuation))
        end
        wc=nothing;pc=nothing;wg=nothing;pg=nothing;GC.gc()
    end
    cross=support_errors(states[:materialized],states[:axis])
    cross.identities_equal && all(<=(1e-4),(cross.position,cross.hidden,cross.metabolites,cross.fields)) || error("dense/axis parity failed")
    (checks=errors,dense_axis=cross)
end
function support_sample(nvox,strategy)
    construction=@elapsed w,p=support_population(nvox,backend,strategy,founders,capacity)
    memory=support_memory(p);before=Interior.device_cellular_statistics(p);after=nothing
    allocated=@allocated elapsed=@elapsed for k in 1:steps
        Interior.device_cellular_step!(w,p,.001;check_status=false)
        if k%8==0 || k==steps;after=Interior.device_cellular_statistics(p);end
    end
    # Biological reads are capacity-sized, not full-volume checkpoint copies.
    biology=(active=Array(p.state.active),ids=Array(p.state.ids),position=Array(p.state.position),
        hidden=Array(p.state.hidden),bond_ids=Array(p.state.bond_ids))
    result=(nvox=nvox,strategy=strategy,steps=steps,dt=.001,seconds=elapsed,
        full_steps_per_second=steps/elapsed,construction_seconds=construction,
        host_allocation_bytes=allocated,memory=memory,initial=before,final=after)
    w=nothing;p=nothing;GC.gc()
    result,biology
end
gate=support_gate();GC.gc()
habitat_gates=Any[]
for strategy in (:materialized,:axis)
    push!(habitat_gates,device_habitat_3d_gate(backend;support_strategy=strategy));GC.gc()
end
samples=Any[];comparisons=Any[]
for nvox in sizes,pair in 1:pairs
    results=Dict{Symbol,Any}();biology=Dict{Symbol,Any}()
    for strategy in (isodd(pair) ? (:materialized,:axis) : (:axis,:materialized))
        sample,state=support_sample(nvox,strategy)
        push!(samples,merge((pair=pair,),sample));results[strategy]=sample;biology[strategy]=state
    end
    a=biology[:materialized];b=biology[:axis]
    same=a.active==b.active && a.ids==b.ids
    push!(comparisons,(nvox=nvox,pair=pair,paired_speedup=results[:materialized].seconds/results[:axis].seconds,
        identities_equal=same,position_error=maximum(abs.(a.position.-b.position)),
        hidden_error=maximum(abs.(a.hidden.-b.hidden)),bonds_equal=a.bond_ids==b.bond_ids))
    empty!(biology);GC.gc()
    open(joinpath(out,"support-pairs.json"),"w") do io
        JSON3.pretty(io,(source=Interior.DEVICE_CELLULAR_SOURCE_PIN,benchmark_sha256=bytes2hex(Interior.SHA.sha256(read(@__FILE__))),source_commit=get(ENV,"DEVICE_SOURCE_COMMIT","unspecified"),runtime=device_cellular_runtime(),
            law=(extent=4,exchange_radius=.18,density=128.,founders=founders,capacity=capacity,mechanical_substeps=16),
            timing_boundary="complete projected steps including birth/death/contact stages and batch8 compact failure/observation downloads; excludes construction, compile gate, final biological reads",gate=gate,habitat_gates=habitat_gates,samples=samples,comparisons=comparisons))
    end
end
println(comparisons)
