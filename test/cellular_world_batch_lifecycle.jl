using Interior,Test,KernelAbstractions,Serialization
isdefined(Main,:CellularWorldBatch)||include(joinpath(@__DIR__,"..","studies","cellular_world_batches.jl"))
isdefined(Main,:constructive_world_batch)||include(joinpath(@__DIR__,"..","studies","constructive_world_batches.jl"))
isdefined(Main,:constructive_batch_fixture)||include(joinpath(@__DIR__,"..","studies","constructive_batch_fixture.jl"))
@testset "Packed lifecycle markers survive birth, death and failure" begin
    fixtures=[constructive_batch_fixture(seed) for seed in (301,302,303)]
    for (w,p) in fixtures;fill!(p.cells.state.weights,0);p.cells.state.weights[1,end,Int(p.cells.geometry.hidden)+10]=10;end
    fixtures[2][2].cells.state.mass[1]=.01f0
    fill!(fixtures[3][2].cells.environment.status,Int32(8))
    batch=constructive_world_batch(last.(fixtures));work=Interior.device_lifecycle_workspace(batch.cells)
    @test batch_lifecycle_validate(batch.cells)
    @test size(work.highwater.storage)==(1,3)
    @test Array(work.highwater.storage)==fill(Int32(2),1,3)
    for step in 1:2
        for (w,p) in fixtures;Interior.device_cellular_step!(w,p.cells,.01;check_status=false,anabolism=p.habitat);end
        constructive_world_batch_step!(first.(fixtures),batch,.01;check_status=false)
        for i in 1:3
            @test only(Array(batch_world(work.highwater,i)))==only(Array(Interior.device_lifecycle_workspace(fixtures[i][2].cells).highwater))
            @test constructive_snapshot_equal(Interior.device_habitat_snapshot(fixtures[i][2];check_status=false),constructive_world_batch_snapshot(batch,i;check_status=false))
        end
    end
    @test work.highwater.storage[1,1]>2
    @test work.highwater.storage[1,3]==2
    checkpoint=constructive_world_batch_checkpoint(first.(fixtures),batch;members=[1,2])
    targets=[World(first(fixtures)[1].grid;seed=i) for i in 1:2]
    for w in targets;fill!(w.fields.nutrient,0);end
    restored=constructive_world_batch_restore(targets,checkpoint)
    @test Array(Interior.device_lifecycle_workspace(restored.cells).highwater.storage)==work.highwater.storage[:,1:2]
    @test_throws ArgumentError batch_lifecycle_workspace(CPU(),8,[-1,1])
    @test_throws ArgumentError batch_lifecycle_workspace(CPU(),8,[9,1])
    bad=deepcopy(batch.cells.base);old=Interior.device_lifecycle_workspace(bad)
    h=bad.environment.status
    alias=Interior.DeviceLifecycleWorkspace(h,bad.backend,old.capacity,(highwater=h,shape=size(h)))
    bad.ownership=merge(bad.ownership,(stages=merge(bad.ownership.stages,(lifecycle=alias,)),))
    @test_throws ArgumentError batch_lifecycle_validate(bad)
end

@testset "Batch carrier preserves Float64 battery representation" begin
    fixtures=map(1:2) do i
        w=World(Grid(3,3;extent=(3,3,1)));fill!(w.fields.nutrient,0)
        r=ResourcePool(w;potentials=(4,4,0),initial=(.2,.3,.4),diffusion=(0,0,0))
        c=CellularEcology(w,reshape([.6+.1i,.7,0.],3,1),r;capacity=4,seed=UInt32(400+i))
        (w,Interior.DeviceCellular(w,c;mechanical_substeps=2,energy_type=Float64))
    end
    batch=CellularWorldBatch(last.(fixtures))
    @test eltype(batch.state.energy)===Float64
    @test batch.geometry.energy_bits==64
    for step in 1:3
        for (w,p) in fixtures;Interior.device_cellular_step!(w,p,.01;lifecycle=:dense);end
        cellular_world_batch_step!(first(fixtures)[1],batch,.01)
        for (i,(_,p)) in enumerate(fixtures)
            @test isequal(Interior.device_cellular_snapshot(p).state,cellular_world_batch_snapshot(batch,i).state)
        end
    end
end

@testset "Lifecycle habitat backing aliases fail before any sibling advances" begin
    fixtures=[constructive_batch_fixture(seed) for seed in (501,502)]
    reference=constructive_world_batch(last.(fixtures));old=reference.cells.base
    for role in (:enzyme_bits,:bound,:raw)
        physical=role===:raw ? reference.habitat.feedstock.arrays.fields : getproperty(reference.habitat.arrays,role)
        parent=eltype(physical.storage)===Int32 ? physical.storage : reinterpret(Int32,physical.storage)
        highwater=WorldBatchArray(view(parent,1:1,1,:),(1,))
        alias=Interior.DeviceLifecycleWorkspace(highwater,old.backend,Int(old.geometry.capacity),(highwater=highwater,shape=size(highwater)))
        ownership=merge(old.ownership,(stages=merge(old.ownership.stages,(lifecycle=alias,)),))
        base=CellularWorldBatch(old.state,old.environment,old.backend,old.parameters,old.geometry,ownership,old.attempted_steps,old.worlds,old.seeds)
        h=reference.habitat
        bad=constructive_batch_carrier(base,h.arrays,h.feedstock.arrays,h.parameters,h.feedstock.parameters,reference.cells.habitat_seeds)
        worlds=first.(fixtures);before=[constructive_world_batch_snapshot(bad,i;check_status=false) for i in 1:2]
        times=getproperty.(worlds,:time);attempts=bad.cells.attempted_steps
        @test batch_lifecycle_validate(base) # This alias is outside base cellular roles.
        @test_throws ArgumentError constructive_world_batch_step!(worlds,bad,.01;check_status=false)
        @test getproperty.(worlds,:time)==times
        @test bad.cells.attempted_steps==attempts
        for i in 1:2
            @test constructive_snapshot_equal(before[i],constructive_world_batch_snapshot(bad,i;check_status=false))
        end
        @test_throws ArgumentError constructive_world_batch_checkpoint(worlds,bad;members=[2])
    end
end
