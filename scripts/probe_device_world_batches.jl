# Run only in a coordinated GPU window. Julia/KA physical kernels are shared.
# DEVICE_BACKEND=amd DEVICE_AMD_OPTLEVEL=1 julia --project=PINNED_ENV this_file
include("device_cellular_setup.jl")
include(joinpath(@__DIR__,"..","studies","cellular_world_batches.jl"))
kind=="amd" && include(joinpath(@__DIR__,"..","studies","cellular_world_batches_amdgpu.jl"))
using Statistics
function batch_probe_fixture(seed,n,side;extent=8.)
    grid=Grid(side,side,side;extent=(extent,extent,extent))
    w=World(grid,backend,NamedTuple(),NamedTuple(),0.,Int(seed))
    k=ceil(Int,cbrt(n));x=zeros(3,n)
    for i in 1:n
        j=i-1;x[:,i].=extent.*((mod(j,k)+.43,mod(div(j,k),k)+.43,div(j,k*k)+.43)./k)
    end
    n>=2 && (x[:,2].=x[:,1].+[1.,0.,0.]) # an actual contact, not empty launches
    p=Interior.device_cellular_initialize(w,x;capacity=max(16,2n),seed=seed,
        support_strategy=:axis,exchange_radius=.8,initial=(.4,.4,0.),diffusion=(.02,.02,.01),
        source_rates=(.01,.01,0.),source_remaining=(1.,1.,0.),mechanical_substeps=16)
    mass=Array(p.state.mass);mass[1]=2f0;n>=3&&(mass[3]=.1f0);copyto!(p.state.mass,mass)
    w,p
end
function batch_probe_parity(serial,batch;atol=2e-5,rtol=2e-5)
    errors=Dict{String,Float64}()
    for (world,(_,p)) in enumerate(serial)
        a=Interior.device_cellular_snapshot(p;check_status=false)
        b=cellular_world_batch_snapshot(batch,world)
        for role in (:state,:environment),key in keys(getproperty(a,role))
            x=getproperty(getproperty(a,role),key);y=getproperty(getproperty(b,role),key)
            # Atomic linked-list insertion order is scratch, not identity.
            (role===:state&&key===:next || role===:environment&&key===:heads)&&continue
            if eltype(x)<:Integer
                if role===:environment&&key===:metadata
                    x=x[1:6];y=y[1:6] # slot7 encodes a Float32 maximum radius
                end
                x==y||error("integer/identity mismatch: world=$world $role.$key")
            else
                all(isfinite,x)&&all(isfinite,y)||error("nonfinite $role.$key")
                err=maximum(abs.(Float64.(x).-Float64.(y));init=0.)
                scale=max(maximum(abs,x;init=0f0),maximum(abs,y;init=0f0))
                errors["$role.$key"]=max(get(errors,"$role.$key",0.),err)
                err<=atol+rtol*scale||error("Float32 parity mismatch: world=$world $role.$key error=$err")
            end
        end
        # Rebuilt bins cover pre-lifecycle survivors; newborns are inserted
        # by the next step geometry stage, as in the ordinary schedule.
        s=b.state;e=b.environment;seen=zeros(Int,length(s.active))
        for head in e.heads
            cell=Int(head);visits=0
            while cell!=0
                1<=cell<=length(seen)||error("cross-world/out-of-range bin link")
                seen[cell]+=1;visits+=1;visits<=length(seen)||error("cyclic bin chain")
                cell=Int(s.next[cell])
            end
        end
        seen==Int.((s.active.==1).&(s.new_parent.==0))||error("bin list does not cover pre-lifecycle world cells")
    end
    errors
end
function batch_probe_case(n,side;worlds=4,steps=8,samples=3)
    serial=[batch_probe_fixture(UInt32(100+i),n,side) for i in 1:worlds]
    batch=CellularWorldBatch(last.(serial));w=first(serial)[1]
    initial=[Interior.device_cellular_audit(cellular_world_batch_snapshot(batch,i)) for i in 1:worlds]
    serial_step()=foreach(wp->Interior.device_cellular_step!(wp[1],wp[2],.01;check_status=false),serial)
    batch_step()=cellular_world_batch_step!(w,batch,.01;check_status=false)
    serial_status()=foreach(wp->Interior.device_cellular_status!(wp[2]),serial)
    serial_step();serial_status();batch_step();Interior.device_cellular_status!(batch)
    batch_probe_parity(serial,batch) # full native-backend parity before timing
    st=Float64[];bt=Float64[]
    for sample in 1:samples
        push!(st,@elapsed(begin for _ in 1:steps;serial_step();end;serial_status() end))
        push!(bt,@elapsed(begin for _ in 1:steps;batch_step();end;Interior.device_cellular_status!(batch) end))
    end
    errors=batch_probe_parity(serial,batch)
    final=[Interior.device_cellular_audit(cellular_world_batch_snapshot(batch,i)) for i in 1:worlds]
    material=[b.material-a.material-b.supplied+a.supplied+b.removed-a.removed for (a,b) in zip(initial,final)]
    potential=[b.potential-a.potential+b.work-a.work+b.heat-a.heat-b.supplied_potential+a.supplied_potential+b.removed_potential-a.removed_potential for (a,b) in zip(initial,final)]
    maxmat=maximum(abs,material);maxenergy=maximum(abs,potential)
    maxmat<max(2e-5,2e-6*maximum(a.material for a in initial))||error("material gate")
    maxenergy<max(2e-5,2e-6*maximum(a.potential for a in initial))||error("potential/work/heat gate")
    (cells_per_world=n,capacity_per_world=Int(batch.geometry.capacity),worlds=worlds,
        grid=Tuple(Int.(batch.geometry.dims)),extent=batch.geometry.extent,exchange_radius=batch.parameters.exchange_radius,
        packed_resident_bytes=sum(sizeof(a.storage) for a in values(batch.state))+sum(sizeof(a.storage) for a in values(batch.environment))+sizeof(batch.seeds),
        seeds=Int.(Array(batch.seeds)),steps_per_sample=steps,samples=samples,
        serial_seconds=st,batch_seconds=bt,aggregate_ratio=median(st)/median(bt),
        material_errors=material,potential_work_heat_errors=potential,field_max_errors=errors,
        completed_steps=[a.completed_steps for a in final],
        timing="all full steps plus endpoint status synchronization; excludes compilation, packing and full audit copies")
end
# Isolate a latched failed world from a healthy peer, including its random seed.
failure_serial=[batch_probe_fixture(UInt32(201+i),4,8) for i in 1:2]
fill!(failure_serial[1][2].environment.status,Int32(8))
failure_batch=CellularWorldBatch(last.(failure_serial))
for (w,p) in failure_serial;Interior.device_cellular_step!(w,p,.01;check_status=false);end
cellular_world_batch_step!(failure_serial[1][1],failure_batch,.01;check_status=false)
Interior.KernelAbstractions.synchronize(backend)
batch_probe_parity(failure_serial,failure_batch)
failed=cellular_world_batch_snapshot(failure_batch,1);healthy=cellular_world_batch_snapshot(failure_batch,2)
failed.environment.status[1]==8&&failed.environment.metadata[6]==0||error("failed-world latch changed")
healthy.environment.status[1]==0&&healthy.environment.metadata[6]==1||error("healthy world froze")
small=batch_probe_case(4,8;worlds=2,steps=2,samples=1)
moderate=batch_probe_case(parse(Int,get(ENV,"BATCH_PROBE_CELLS","64")),16;
    worlds=parse(Int,get(ENV,"BATCH_PROBE_WORLDS","4")),steps=parse(Int,get(ENV,"BATCH_PROBE_STEPS","8")))
result=(runtime=device_cellular_runtime(),source=CELLULAR_WORLD_BATCH_SOURCE,
    failure_isolation=true,small_native_fixture=small,moderate=moderate,gpu=(kind!="cpu"))
path=get(ENV,"BATCH_PROBE_RESULT","outputs/device-world-batches-probe.json")
ispath(path)&&error("fresh result path required")
mkpath(dirname(path));open(path,"w") do io;JSON3.pretty(io,result);println(io);end
println(result)
