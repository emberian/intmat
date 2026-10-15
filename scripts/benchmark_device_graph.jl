#!/usr/bin/env julia
# Exclusive AMD window only: alternating complete native/graph trajectories.
include("device_cellular_setup.jl")
using Statistics
using Interior.KernelAbstractions: @kernel, @index
backend isa AMDGPU.ROCBackend || error("native AMD backend required")
isdefined(Interior,:DeviceCellularGraph) || Base.include(Interior,joinpath(@__DIR__,"..","studies","amd_cellular_graph.jl"))
length(ARGS)==1 || error("fresh result JSON path required")
ispath(ARGS[1]) && error("fresh result path required")
N=parse(Int,get(ENV,"DEVICE_GRAPH_FOUNDERS","512"))
C=parse(Int,get(ENV,"DEVICE_CAPACITY","8192"))
G=parse(Int,get(ENV,"DEVICE_GRAPH_GRID","32"))
steps=parse(Int,get(ENV,"DEVICE_GRAPH_STEPS","30"))
pairs=parse(Int,get(ENV,"DEVICE_GRAPH_PAIRS","3"))
N>=2 && C>=2N && G>=8 && 1<=steps<=100 && 1<=pairs<=5 || error("invalid bounded graph benchmark")
# Closed concentration .3/species in physical volume64 plus founder biomass/stores.
# This is a conservative lifetime population bound, not merely C=2N headroom.
material_forecast=64*.6+N*(1.85+.24)
C*.25>=material_forecast || error("capacity cannot hold conservative closed material population bound")
@kernel function graph_population_prior!(s,h,n)
    i=@index(Global,Linear)
    if i<=n
        s.weights[i,end,h+8]=2f0;s.weights[i,end,h+9]=2f0;s.weights[i,end,h+10]=2f0
        i==1 && (s.mass[i]=.1f0)
    end
end
function graph_population(n)
    grid=Grid(G,G,G;extent=(4,4,4));w=World(grid,backend,NamedTuple(),NamedTuple(),0.,911)
    side=ceil(Int,cbrt(n));positions=zeros(3,n)
    for (i,xyz) in enumerate(Iterators.product(0:side-1,0:side-1,0:side-1))
        i>n && break
        positions[:,i].=4 .* (collect(xyz).+.5)./side
    end
    p=Interior.device_cellular_initialize(w,positions;capacity=C,seed=911,support_strategy=:axis,
        exchange_radius=.18,density=128.,biomass=1.85,energy=1.,mutation=.01,
        initial=(.3,.3,0.),diffusion=(.001,.001,.001),mechanical_substeps=16)
    graph_population_prior!(backend,128)(p.state,p.geometry.hidden,n;ndrange=C)
    AMDGPU.synchronize()
    w,p
end
# Warm both scalar and dynamic-clock signatures without touching timed states.
ww,wp=graph_population(min(N,8))
Interior.device_cellular_step!(ww,wp,.001)
wp.attempted_steps+=1
Interior.device_cellular_enqueue!(ww,wp,.001f0,:projected,nothing,Interior.DeviceGraphStepIndex(wp.environment.metadata))
Interior.device_cellular_status!(wp)
ww=nothing;wp=nothing;GC.gc();AMDGPU.synchronize()
function graph_nodes(owner)
    count=Ref{Csize_t}(0)
    AMDGPU.HIP.hipGraphGetNodes(owner.graph,C_NULL,count)
    Int(count[])
end
function run_one(mode)
    w,p=graph_population(N)
    initial=Interior.device_cellular_audit(Interior.device_cellular_snapshot(p))
    capture_seconds=0.;owner=nothing;nodes=0
    if mode===:graph
        started=time_ns();owner=Interior.device_cellular_capture_graph(w,p,.001,AMDGPU)
        capture_seconds=(time_ns()-started)/1e9;nodes=graph_nodes(owner)
    end
    AMDGPU.synchronize()
    begin_event=AMDGPU.HIP.HIPEvent(AMDGPU.stream();do_record=false,timing=true)
    end_event=AMDGPU.HIP.HIPEvent(AMDGPU.stream();do_record=false,timing=true)
    wall_start=time_ns();AMDGPU.HIP.record(begin_event);enqueue_start=time_ns()
    if mode===:graph
        Interior.device_cellular_graph_replay!(owner,steps;check_status=false)
    else
        Interior.device_graph_preflight(w,p,nothing)
        for _ in 1:steps
            Interior.device_cellular_step!(w,p,.001;check_status=false)
        end
    end
    enqueue_seconds=(time_ns()-enqueue_start)/1e9
    AMDGPU.HIP.record(end_event);AMDGPU.HIP.synchronize(end_event)
    Interior.device_cellular_status!(p)
    wall_seconds=(time_ns()-wall_start)/1e9
    stream_seconds=Float64(AMDGPU.HIP.elapsed(begin_event,end_event))
    snapshot=Interior.device_cellular_snapshot(p);final=Interior.device_cellular_audit(snapshot)
    material_error=final.material-initial.material-final.supplied+initial.supplied+final.removed-initial.removed
    potential_error=final.potential-initial.potential+final.work-initial.work+final.heat-initial.heat-
        final.supplied_potential+initial.supplied_potential+final.removed_potential-initial.removed_potential
    final.population>0 && snapshot.environment.metadata[5]==0 || error("extinction or capacity block invalidates declared case")
    abs(material_error)<.002 && abs(potential_error)<.02 || error("graph benchmark accounting failed")
    if mode===:graph
        checkpoint=Interior.device_cellular_checkpoint(w,p)
        rw=World(w.grid,backend,NamedTuple(),NamedTuple(),0.,911)
        restored=Interior.device_cellular_restore(rw,checkpoint)
        Interior.device_cellular_snapshot(restored)==snapshot || error("graph strict restore failed")
    end
    record=(;mode,steps,wall_seconds,enqueue_seconds,stream_seconds,capture_seconds,nodes,
        full_steps_per_second=steps/wall_seconds,initial,final,material_error,potential_error)
    record,snapshot
end
function compare_trajectories(a,b)
    for key in (:active,:ids,:parents,:lineages,:generations,:bonds,:bond_ids)
        a.state[key]==b.state[key] || error("graph trajectory identity/contact mismatch: $key")
    end
    a.environment.metadata==b.environment.metadata || error("graph lifecycle counters differ")
    errors=Dict(string(key)=>maximum(abs,Float64.(a.state[key]).-Float64.(b.state[key]);init=0.)
        for key in (:position,:mass,:metabolites,:energy,:age,:hidden,:adaptation,:weights,:marks,:actions))
    errors["fields"]=maximum(abs,Float64.(a.environment.fields).-Float64.(b.environment.fields);init=0.)
    maximum(values(errors))<1e-4 || error("graph/native primary trajectory tolerance exceeded")
    errors
end
records=Any[];comparisons=Any[]
for pair in 1:pairs
    order=isodd(pair) ? (:native,:graph) : (:graph,:native)
    states=Dict{Symbol,Any}()
    for mode in order
        record,state=run_one(mode);push!(records,merge(record,(;pair)));states[mode]=state
        GC.gc();AMDGPU.synchronize()
    end
    push!(comparisons,(;pair,errors=compare_trajectories(states[:native],states[:graph])))
    empty!(states);GC.gc();AMDGPU.synchronize()
end
native=median([r.wall_seconds for r in records if r.mode===:native])
graph=median([r.wall_seconds for r in records if r.mode===:graph])
source=Dict(name=>bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"..",name))))
    for name in ("studies/device_graph_clock.jl","studies/amd_cellular_graph.jl", "scripts/benchmark_device_graph.jl"))
result=(runtime=device_cellular_runtime(),source,base_source=Interior.DEVICE_CELLULAR_SOURCE_PIN,
    founders=N,capacity=C,grid=G,steps,pairs,dt=.001,mechanical_substeps=16,material_forecast,
    median_native_seconds=native,median_graph_seconds=graph,whole_step_speedup=native/graph,
    records,comparisons,
    timing_boundary="complete native/graph batches including one equal preflight and final status transfer; construction, capture, compilation and full audits excluded",
    stream_timing_note="HIP event interval includes device work, preflight transfer and inter-kernel gaps; it is not an isolated kernel-time or launch-gap attribution")
mkpath(dirname(abspath(ARGS[1])));open(io->JSON3.pretty(io,result),ARGS[1],"w");println(JSON3.write(result))
