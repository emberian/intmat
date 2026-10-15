using Interior,Test,KernelAbstractions
isdefined(@__MODULE__,:constructive_device_fixture) || include(joinpath(@__DIR__,"..","scripts","device_habitat_fixture.jl"))
function lifecycle_workspace_canonical(x)
    names=filter(!=(:free),keys(x.cells.environment))
    merge(x,(cells=merge(x.cells,(environment=NamedTuple{names}(x.cells.environment),)),))
end
function lifecycle_force_division!(p;death_slot=0)
    mass=Array(p.cells.state.mass);mass[1:2].=3f0
    death_slot>0 && (mass[death_slot]=.025f0)
    copyto!(p.cells.state.mass,mass)
    weights=Array(p.cells.state.weights);h=Int(p.cells.geometry.hidden)
    for i in 1:2;weights[i,:,h+10].=0;weights[i,end,h+10]=5f0;end
    copyto!(p.cells.state.weights,weights)
end
@testset "Owned lifecycle watermark preflight and aliases" begin
    w,p=constructive_device_fixture();c=p.cells;work=Interior.device_lifecycle_workspace(c)
    @test work.highwater==Int32[2] && work.capacity==16
    @test Interior.device_lifecycle_workspace_validate(work,c.backend,16)
    @test_throws ArgumentError Interior.DeviceLifecycleWorkspace(CPU(),16;highwater=zeros(Int32,2))
    @test_throws ArgumentError Interior.DeviceLifecycleWorkspace(CPU(),16;initial=17)
    before=Interior.device_habitat_snapshot(p)
    resource=Interior.device_resource_workspace(c)
    for array in (view(c.state.ids,1:1),view(c.environment.status,1:1),view(c.environment.metadata,1:1),
            view(reinterpret(Int32,vec(resource.grants)),1:1))
        alias=Interior.DeviceLifecycleWorkspace(CPU(),16;highwater=array)
        @test_throws ArgumentError Interior.device_cellular_contract(c.state,c.environment,c.parameters,c.geometry;
            resource_workspace=resource,lifecycle_workspace=alias)
    end
    @test before==Interior.device_habitat_snapshot(p)
    forged=Interior.DeviceLifecycleWorkspace(copy(work.highwater),work.backend,work.capacity,work.contract)
    ownership=c.ownership
    c.ownership=merge(ownership,(stages=merge(ownership.stages,(lifecycle=forged,)),))
    @test_throws ArgumentError Interior.device_cellular_step!(w,c,.01)
    c.ownership=ownership
    resize!(work.highwater,2)
    fields=deepcopy(w.fields)
    @test_throws ArgumentError step!(w,(Deposition(threshold=0),c),.01)
    @test fields==w.fields && before==Interior.device_habitat_snapshot(p;check_status=false)
    resize!(work.highwater,1)
    @test_throws ArgumentError Interior.device_cellular_step!(w,c,.01;lifecycle=:invalid)
    @test c.attempted_steps==0
end
@testset "Production prefix and dense oracle preserve both battery precisions" begin
    for precision in (Float32,Float64),strategy in (:materialized,:axis)
        w,p=constructive_device_fixture(;energy_type=precision,support_strategy=strategy)
        lifecycle_force_division!(p;death_slot=2)
        cp=Interior.device_habitat_checkpoint(w,p)
        w2=World(w.grid;seed=w.seed);fill!(w2.fields.nutrient,0)
        p2=Interior.device_habitat_restore(w2,cp)
        @test eltype(p.cells.state.energy)===precision && cp.cells.source.lifecycle==Interior.DEVICE_CELLULAR_SOURCE_PIN.lifecycle
        @test !(:highwater in keys(cp.cells.snapshot.state)) && !(:highwater in keys(cp.cells.snapshot.environment))
        for step in 1:6
            Interior.device_cellular_step!(w,p.cells,.01;anabolism=p.habitat,lifecycle=:dense)
            Interior.device_cellular_step!(w2,p2.cells,.01;anabolism=p2.habitat)
            @test lifecycle_workspace_canonical(Interior.device_habitat_snapshot(p))==lifecycle_workspace_canonical(Interior.device_habitat_snapshot(p2))
            @test Interior.device_lifecycle_workspace(p.cells).highwater==Interior.device_lifecycle_workspace(p2.cells).highwater
        end
        @test p2.cells.environment.metadata[3]>=1 && p2.cells.environment.metadata[4]>=1
        @test p2.cells.state.ids[2]!=2 && p2.cells.state.parents[2]==1
        # Switching oracle recipes must keep the owned watermark current.
        Interior.device_cellular_step!(w,p.cells,.01;anabolism=p.habitat)
        Interior.device_cellular_step!(w2,p2.cells,.01;anabolism=p2.habitat,lifecycle=:dense)
        @test lifecycle_workspace_canonical(Interior.device_habitat_snapshot(p))==lifecycle_workspace_canonical(Interior.device_habitat_snapshot(p2))
    end
end
@testset "Checkpoint reconstructs dead-tail watermark and owns fresh storage" begin
    w,p=constructive_device_fixture(;energy_type=Float64)
    lifecycle_force_division!(p)
    Interior.device_cellular_step!(w,p.cells,.01;anabolism=p.habitat)
    @test p.cells.environment.metadata[3]==2 && only(Interior.device_lifecycle_workspace(p.cells).highwater)==4
    p.cells.state.mass[4]=.025f0
    Interior.device_cellular_step!(w,p.cells,.01;anabolism=p.habitat)
    @test p.cells.state.active[4]==0 && p.cells.state.ids[4]>0
    cp=Interior.device_habitat_checkpoint(w,p)
    w2=World(w.grid;seed=w.seed);fill!(w2.fields.nutrient,0)
    p2=Interior.device_habitat_restore(w2,cp)
    a=Interior.device_lifecycle_workspace(p.cells);b=Interior.device_lifecycle_workspace(p2.cells)
    @test a.highwater==b.highwater==Int32[4] && !Base.mightalias(a.highwater,b.highwater)
    @test Interior.device_habitat_snapshot(p)==Interior.device_habitat_snapshot(p2)
    for _ in 1:3
        Interior.device_cellular_step!(w,p.cells,.01;anabolism=p.habitat)
        Interior.device_cellular_step!(w2,p2.cells,.01;anabolism=p2.habitat)
    end
    @test Interior.device_habitat_snapshot(p)==Interior.device_habitat_snapshot(p2)
    p2.cells.environment.status[1]=8
    before=Interior.device_habitat_snapshot(p2;check_status=false);watermark=copy(b.highwater)
    Interior.device_cellular_step!(w2,p2.cells,.01;anabolism=p2.habitat,check_status=false)
    after=Interior.device_habitat_snapshot(p2;check_status=false)
    @test before.cells.state==after.cells.state && before.cells.environment==after.cells.environment && b.highwater==watermark
    @test_throws ArgumentError Interior.device_habitat_checkpoint(w2,p2)
end
