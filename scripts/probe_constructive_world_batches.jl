# Bounded native composition gate. Launch AMD only in a root-granted GPU window.
include("device_cellular_setup.jl")
include(joinpath(@__DIR__,"..","studies","constructive_world_batches.jl"))
include(joinpath(@__DIR__,"..","studies","constructive_batch_fixture.jl"))
kind=="amd"&&include(joinpath(@__DIR__,"..","studies","cellular_world_batches_amdgpu.jl"))
using Statistics, Serialization
function constructive_probe_compare(a,b,errors=Dict{String,Float64}(),path="";atol=2e-5,rtol=2e-5)
    if a isa NamedTuple
        keys(a)==keys(b)||error("role mismatch $path")
        for key in keys(a);constructive_probe_compare(getproperty(a,key),getproperty(b,key),errors,"$path.$key";atol,rtol);end
    elseif a isa AbstractArray
        size(a)==size(b)&&eltype(a)==eltype(b)||error("array schema mismatch $path")
        path in (".cells.state.next",".cells.environment.heads")&&return errors
        if path==".arrays.enzyme_bits"
            return constructive_probe_compare(reinterpret(Float32,a),reinterpret(Float32,b),errors,path*"[Float32]";atol,rtol)
        elseif path==".cells.environment.metadata"
            a[[1:6;8:length(a)]]==b[[1:6;8:length(b)]]||error("identity/counter mismatch $path")
            return constructive_probe_compare([reinterpret(Float32,a[7])],[reinterpret(Float32,b[7])],errors,path*".radius";atol,rtol)
        elseif eltype(a)<:AbstractFloat
            all(isfinite,a)&&all(isfinite,b)||error("nonfinite $path")
            delta=maximum(abs.(Float64.(a).-Float64.(b));init=0.)
            scale=max(maximum(abs,a;init=0f0),maximum(abs,b;init=0f0))
            errors[path]=max(get(errors,path,0.),delta)
            delta<=atol+rtol*scale||error("Float32 mismatch $path delta=$delta")
        else
            a==b||error("integer state mismatch $path")
        end
    else
        isequal(a,b)||error("parameter/state mismatch $path")
    end
    errors
end
function constructive_probe_worlds(batch)
    g=batch.cells.geometry
    [World(Grid(Tuple(Int.(g.dims))...;extent=g.extent),batch.cells.backend,NamedTuple(),NamedTuple(),0.,i) for i in 1:batch.cells.worlds]
end
function constructive_probe_all(serial,batch)
    errors=Dict{String,Float64}()
    for (i,(_,p)) in enumerate(serial)
        a=Interior.device_habitat_snapshot(p;check_status=false)
        b=constructive_world_batch_snapshot(batch,i;check_status=false)
        constructive_probe_compare(a,b,errors)
        s=b.cells.state;e=b.cells.environment;seen=zeros(Int,length(s.active))
        for head in e.heads
            cell=Int(head);visits=0
            while cell!=0
                1<=cell<=length(seen)||error("out-of-world constructive bin link")
                seen[cell]+=1;visits+=1;visits<=length(seen)||error("cyclic constructive bin link")
                cell=Int(s.next[cell])
            end
        end
        seen==Int.((s.active.==1).&(s.new_parent.==0))||error("constructive bin ownership/coverage mismatch")
    end
    errors
end
function constructive_probe_case(n,side;worlds=4,steps=4,samples=3,hidden=8)
    serial=[constructive_batch_fixture(100+i;flux=true,backend,shape=(side,side,side),extent=(8.,8.,8.),
        population=n,capacity=max(16,4n),hidden,mechanics=true) for i in 1:worlds]
    batch=constructive_world_batch(last.(serial));clocks=constructive_probe_worlds(batch)
    initial=[Interior.device_habitat_audit(constructive_world_batch_snapshot(batch,i)) for i in 1:worlds]
    serial_step()=foreach(wp->Interior.device_cellular_step!(wp[1],wp[2].cells,.01;check_status=false,anabolism=wp[2].habitat),serial)
    batch_step()=Interior.device_cellular_step!(clocks[1],batch.cells,.01;check_status=false,anabolism=batch.habitat)
    serial_status()=foreach(wp->Interior.device_cellular_status!(wp[2].cells),serial)
    serial_step();serial_status();batch_step();Interior.device_cellular_status!(batch.cells)
    constructive_probe_all(serial,batch)
    st=Float64[];bt=Float64[]
    for sample in 1:samples
        push!(st,@elapsed(begin for _ in 1:steps;serial_step();end;serial_status() end))
        push!(bt,@elapsed(begin for _ in 1:steps;batch_step();end;Interior.device_cellular_status!(batch.cells) end))
    end
    errors=constructive_probe_all(serial,batch)
    final=[Interior.device_habitat_audit(constructive_world_batch_snapshot(batch,i)) for i in 1:worlds]
    material=[b.material-a.material-b.supplied+a.supplied+b.removed-a.removed for (a,b) in zip(initial,final)]
    potential=[b.potential-a.potential+b.work-a.work+b.heat-a.heat-b.supplied_potential+a.supplied_potential+b.removed_potential-a.removed_potential for (a,b) in zip(initial,final)]
    maximum(abs,material)<max(2e-5,2e-6*maximum(a.material for a in initial))||error("constructive material gate")
    maximum(abs,potential)<max(2e-5,2e-6*maximum(a.potential for a in initial))||error("constructive potential/work/heat gate")
    # Checkpoint/restore and one post-restart matched step are outside timing.
    for i in 1:worlds;clocks[i].time=Float64(.01f0)*final[i].completed_steps;end
    checkpoint=constructive_world_batch_checkpoint(clocks,batch)
    io=IOBuffer();serialize(io,checkpoint);seekstart(io)
    restored=constructive_world_batch_restore(constructive_probe_worlds(batch),deserialize(io))
    for i in 1:worlds
        constructive_probe_compare(constructive_world_batch_snapshot(batch,i),constructive_world_batch_snapshot(restored,i),errors)
    end
    serial_step();serial_status();batch_step();Interior.device_cellular_status!(batch.cells)
    Interior.device_cellular_step!(clocks[1],restored.cells,.01;anabolism=restored.habitat)
    constructive_probe_all(serial,batch);constructive_probe_all(serial,restored)
    (worlds=worlds,cells_per_world=n,capacity_per_world=Int(batch.cells.geometry.capacity),grid=(side,side,side),extent=(8,8,8),
        exchange_radius=batch.cells.parameters.exchange_radius,hidden=hidden,mechanical_substeps=16,
        cell_seeds=Int.(Array(batch.cells.seeds)),habitat_seeds=Int.(Array(batch.cells.habitat_seeds)),
        steps_per_sample=steps,serial_seconds=st,batch_seconds=bt,aggregate_ratio=median(st)/median(bt),
        field_max_errors=errors,material_errors=material,potential_work_heat_errors=potential,
        audits=final,serialized_checkpoint_bytes=position(io),restart_parity=true,
        timing="complete coupled enqueue + endpoint status; excludes construction, JIT, snapshots, checkpoint and restart")
end
# First isolate failure states in the coupled raw/free/bound/matrix composition.
failed_serial=[constructive_batch_fixture(200+i;backend,flux=true,shape=(6,6,6),extent=(6.,6.,6.),mechanics=true) for i in 1:2]
fill!(failed_serial[1][2].cells.environment.status,Int32(8))
failed_batch=constructive_world_batch(last.(failed_serial));failed_worlds=constructive_probe_worlds(failed_batch)
for (w,p) in failed_serial;Interior.device_cellular_step!(w,p.cells,.01;check_status=false,anabolism=p.habitat);end
constructive_world_batch_step!(failed_worlds,failed_batch,.01;check_status=false)
constructive_probe_all(failed_serial,failed_batch)
healthy_checkpoint=constructive_world_batch_checkpoint(failed_worlds,failed_batch;members=[2])
length(healthy_checkpoint.members)==1&&healthy_checkpoint.member_indices==[2]||error("healthy member recovery failed")
small=constructive_probe_case(2,6;worlds=2,steps=1,samples=1,hidden=2)
moderate=constructive_probe_case(parse(Int,get(ENV,"CONSTRUCTIVE_PROBE_CELLS","64")),16;
    worlds=parse(Int,get(ENV,"CONSTRUCTIVE_PROBE_WORLDS","4")),steps=parse(Int,get(ENV,"CONSTRUCTIVE_PROBE_STEPS","4")))
result=(runtime=device_cellular_runtime(),source=CONSTRUCTIVE_BATCH_SOURCE,gpu=kind!="cpu",failure_isolation=true,
    healthy_member_checkpoint=true,small=small,moderate=moderate)
out=get(ENV,"CONSTRUCTIVE_PROBE_RESULT","outputs/constructive-world-batches-native.json")
ispath(out)&&error("fresh result path required");mkpath(dirname(out))
open(out,"w") do io;JSON3.pretty(io,result);println(io);end
println((gpu=result.gpu,small_ratio=small.aggregate_ratio,moderate_ratio=moderate.aggregate_ratio,
    maximum_float_error=maximum(values(moderate.field_max_errors)),restart=true))
