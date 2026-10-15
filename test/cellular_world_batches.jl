using Interior, Test
include(joinpath(@__DIR__,"..","studies","cellular_world_batches.jl"))
function batch_fixture(offset)
    w=World(Grid(3,3;extent=(3,3,1)));fill!(w.fields.nutrient,0)
    r=ResourcePool(w;potentials=(4,4,0),initial=(.2,.3,.4),diffusion=(0,0,0))
    c=CellularEcology(w,reshape([.6+offset,.7,0.],3,1),r;capacity=4,seed=UInt32(11+round(Int,100offset)))
    w,Interior.DeviceCellular(w,c;mechanical_substeps=2)
end
@testset "Independent worlds: complete schedule identity" begin
    fixtures=[batch_fixture(.1i) for i in 1:3];ps=last.(fixtures)
    batch=CellularWorldBatch(ps)
    @test length(unique(Array(batch.seeds)))==3
    for step in 1:5
        for (w,p) in fixtures;Interior.device_cellular_step!(w,p,.01);end
        cellular_world_batch_step!(first(fixtures)[1],batch,.01)
        for (i,p) in enumerate(ps),role in (:state,:environment)
            serial=getproperty(p,role);packed=getproperty(batch,role)
            @test all(isequal(getproperty(serial,k),Array(batch_world(getproperty(packed,k),i))) for k in keys(serial))
        end
    end
end
@testset "Failure latch cannot leak between worlds" begin
    fixtures=[batch_fixture(.1i) for i in 1:2];ps=last.(fixtures)
    ps[1].environment.status[1]=8
    batch=CellularWorldBatch(ps)
    for (w,p) in fixtures;Interior.device_cellular_step!(w,p,.01;check_status=false);end
    cellular_world_batch_step!(first(fixtures)[1],batch,.01;check_status=false)
    for (i,p) in enumerate(ps),role in (:state,:environment)
        serial=getproperty(p,role);packed=getproperty(batch,role)
        @test all(isequal(getproperty(serial,k),Array(batch_world(getproperty(packed,k),i))) for k in keys(serial))
    end
    @test batch.environment.status.storage[1,2]==0
    @test_throws ArgumentError cellular_world_batch_step!(first(fixtures)[1],batch,.01;anabolism=:unsupported)
end
@testset "Birth, death and compensated ledgers stay world-local" begin
    fixtures=[batch_fixture(.1i) for i in 1:3];ps=last.(fixtures)
    for p in ps
        fill!(p.state.weights,0)
        p.state.weights[1,end,Int(p.geometry.hidden)+10]=10
    end
    ps[1].state.mass[1]=2f0;ps[2].state.mass[1]=.01f0
    batch=CellularWorldBatch(ps)
    initial=[Interior.device_cellular_audit(cellular_world_batch_snapshot(batch,i)) for i in 1:3]
    for (w,p) in fixtures;Interior.device_cellular_step!(w,p,.01);end
    cellular_world_batch_step!(first(fixtures)[1],batch,.01)
    @test ps[1].environment.metadata[3]>=1
    @test ps[2].environment.metadata[4]>=1
    for (i,p) in enumerate(ps)
        snapshot=cellular_world_batch_snapshot(batch,i)
        @test isequal(snapshot.state,Interior.device_cellular_snapshot(p).state)
        @test isequal(snapshot.environment,Interior.device_cellular_snapshot(p).environment)
        a=initial[i]
        if any(==(1),snapshot.state.active)
            b=Interior.device_cellular_audit(snapshot)
        else
            e=snapshot.environment;q=snapshot.parameters
            # Existing host audit rejects empty Float32 active-cell reductions.
            b=(material=sum(Float64,e.fields)+sum(Float64,e.inert)+sum(Float64,e.inert_low),
               potential=sum(sum(Float64,e.fields[:,k])*q.potentials[k] for k in 1:3),
               work=sum(Float64.(snapshot.state.ledger[:,1]).+Float64.(snapshot.state.ledger_low[:,1])),
               heat=sum(Float64.(snapshot.state.ledger[:,2]).+Float64.(snapshot.state.ledger_low[:,2])),
               supplied=0.,removed=0.,supplied_potential=0.,removed_potential=0.)
        end
        @test abs(b.material-a.material-(b.supplied-a.supplied)+(b.removed-a.removed))<1e-5
        @test abs(b.potential-a.potential+(b.work-a.work)+(b.heat-a.heat)-(b.supplied_potential-a.supplied_potential)+(b.removed_potential-a.removed_potential))<1e-5
    end
end
@testset "Contiguous world slices respect bounds and ownership" begin
    a=WorldBatchArray(reshape(collect(1:12),2,3,2),(2,3));x=batch_world(a,2)
    @test Array(x)==reshape(collect(7:12),2,3)
    @test_throws BoundsError x[7]
    @test_throws BoundsError setindex!(x,0,7)
    x[2,3]=99
    @test a.storage[2,3,2]==99
    @test a.storage[2,3,1]==6
end
@testset "Bulk audit snapshots own their world slice" begin
    x=WorldBatchArray(reshape(collect(1:12),2,3,2),(2,3))
    e=WorldBatchArray(zeros(Int32,1,2),(1,))
    p=CellularWorldBatch((x=x,),(status=e,),CPU(),(seed=UInt32(1),),
        (dims=(2,3,1),),(unused=0,),0,2,UInt32[5,9])
    snapshot=cellular_world_batch_snapshot(p,2)
    @test snapshot.state.x==reshape(collect(7:12),2,3)
    @test snapshot.parameters.seed==9
    snapshot.state.x[1]=999
    @test x.storage[1,1,2]==7
    @test snapshot.environment.status==[0]
    @test_throws BoundsError cellular_world_batch_snapshot(p,3)
end

@testset "Resource reduction workspace has disjoint world slices" begin
    fixtures=[batch_fixture(.1i) for i in 1:3]
    batch=CellularWorldBatch(last.(fixtures));work=Interior.device_resource_workspace(batch)
    @test batch_resource_validate(batch)
    arrays=(work.grants,(a for level in work.levels for a in (level.high,level.low))...)
    for a in arrays
        @test size(a.storage)==(size(a)...,3)
        for world in 1:3;fill!(batch_world(a,world),Float32(world));end
        @test all(all(==(Float32(world)),Array(batch_world(a,world))) for world in 1:3)
        @test all(batch_world(a,i).offset+length(a)<=batch_world(a,i+1).offset for i in 1:2)
        fill!(a.storage,0)
    end
    @test all(!Base.mightalias(a.storage,b.storage) for (i,a) in enumerate(arrays) for b in arrays[1:i-1])
end
