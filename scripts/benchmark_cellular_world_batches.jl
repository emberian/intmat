# JULIA_NUM_THREADS=1 julia --project scripts/benchmark_cellular_world_batches.jl
using Interior, Statistics, TOML
isdefined(Main,:CellularWorldBatch) || include(joinpath(@__DIR__,"..","studies","cellular_world_batches.jl"))
function fixture(i,n,side,extent)
    w=World(Grid(side,side;extent=(extent,extent,1)));fill!(w.fields.nutrient,0)
    r=ResourcePool(w;potentials=(4,4,0),initial=(2.,2.,0.),diffusion=(0,0,0))
    positions=zeros(3,n)
    for j in 1:n
        positions[1,j]=mod(.31j+.1i,extent);positions[2,j]=mod(.71j+.07i,extent)
    end
    c=CellularEcology(w,positions,r;capacity=max(16,2n),seed=UInt32(11+i))
    w,Interior.DeviceCellular(w,c;mechanical_substeps=16)
end
function run_benchmark(;worlds=4,n=8,side=8,steps=10,samples=3,extent=4.)
    fixtures=[fixture(i,n,side,extent) for i in 1:worlds];batch=CellularWorldBatch(last.(fixtures));w=first(fixtures)[1]
    initial=map(p->Interior.device_cellular_audit(Interior.device_cellular_snapshot(p)),last.(fixtures))
    serial_step()=foreach(wp->Interior.device_cellular_step!(wp[1],wp[2],.01),fixtures)
    batch_step()=cellular_world_batch_step!(w,batch,.01)
    serial_step();batch_step() # full schedule compilation, excluded
    serial=Float64[];batched=Float64[]
    for sample in 1:samples
        push!(serial,@elapsed(begin for _ in 1:steps;serial_step();end end))
        push!(batched,@elapsed(begin for _ in 1:steps;batch_step();end end))
    end
    parity=all(isequal(getproperty(getproperty(p,role),k),Array(batch_world(getproperty(getproperty(batch,role),k),i))) for (i,(_,p)) in enumerate(fixtures),role in (:state,:environment) for k in keys(getproperty(p,role)))
    serial_allocated=@allocated serial_step();batch_allocated=@allocated batch_step()
    final=[Interior.device_cellular_audit(cellular_world_batch_snapshot(batch,i)) for i in 1:worlds]
    material_errors=[b.material-a.material-(b.supplied-a.supplied)+(b.removed-a.removed) for (a,b) in zip(initial,final)]
    energy_errors=[b.potential-a.potential+(b.work-a.work)+(b.heat-a.heat)-(b.supplied_potential-a.supplied_potential)+(b.removed_potential-a.removed_potential) for (a,b) in zip(initial,final)]
    result=Dict("seed_vector"=>Int.(Array(batch.seeds)),"physical_extent"=>[extent,extent,1.],"completed_steps_per_world"=>batch.attempted_steps,"serial_step_allocated"=>serial_allocated,"batch_step_allocated"=>batch_allocated,"material_balance_errors"=>material_errors,"potential_work_heat_errors"=>energy_errors,"worlds"=>worlds,"initial_cells_per_world"=>n,"grid"=>[side,side,1],"steps_per_sample"=>steps,"samples"=>samples,"threads"=>Threads.nthreads(),"serial_seconds"=>serial,"batch_seconds"=>batched,"aggregate_speedup"=>median(serial)/median(batched),"all_arrays_exact"=>parity,"source"=>Dict(string(k)=>v for (k,v) in pairs(CELLULAR_WORLD_BATCH_SOURCE)))
    TOML.print(stdout,result);println();parity||error("complete schedule parity failed")
    result
end
run_benchmark(worlds=parse(Int,get(ENV,"BATCH_WORLDS","4")),n=parse(Int,get(ENV,"BATCH_CELLS","8")),side=parse(Int,get(ENV,"BATCH_SIDE","8")),steps=parse(Int,get(ENV,"BATCH_STEPS","10")))
