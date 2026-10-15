using Interior, Test, KernelAbstractions
isdefined(Interior,:DeviceGraphStepIndex) || Base.include(Interior,joinpath(@__DIR__,"..","studies","device_graph_clock.jl"))
@kernel function graph_clock_rng_probe!(out,metadata,clock,seed)
    i=@index(Global,Linear)
    out[i,1]=Interior.dc_random(seed,i,clock,11)
    out[i,2]=Interior.dc_normal(seed,i,clock,19)
end
@testset "Device graph RNG clock reads the current complete-step count" begin
    metadata=zeros(Int32,7);clock=Interior.DeviceGraphStepIndex(metadata)
    output=zeros(Float32,17,2)
    for complete in (0,1,19,10000,typemax(Int32)-1)
        metadata[6]=complete
        graph_clock_rng_probe!(CPU(),64)(output,metadata,clock,UInt32(911);ndrange=17)
        KernelAbstractions.synchronize(CPU())
        for i in 1:17
            @test output[i,1]===Interior.dc_random(UInt32(911),i,UInt32(complete+1),11)
            @test output[i,2]===Interior.dc_normal(UInt32(911),i,UInt32(complete+1),19)
        end
    end
end

isdefined(Interior,:device_graph_enqueue_prototype!) || Base.include(Interior,joinpath(@__DIR__,"..","studies","device_graph_enqueue.jl"))
function graph_cpu_fixture()
    w=World(Grid(8,8,8;extent=(4,4,4));seed=911);fill!(w.fields.nutrient,0)
    positions=Float64[1 3;2 2;1 3]
    p=Interior.device_cellular_initialize(w,positions;capacity=16,seed=911,biomass=1.85,
        exchange_radius=.18,density=128.,mutation=.03,mutation_probability=1.,
        initial=(.3,.3,0.),diffusion=(.001,.001,.001),mechanical_substeps=16,
        source_remaining=(0,0,0),sink_remaining=(0,0,0),source_rates=(0,0,0),sink_rates=(0,0,0))
    p.state.weights[:,end,Int(p.geometry.hidden)+10].=2f0
    w,p
end
@testset "Private captured scheduling preserves full native trajectory" begin
    wr,reference=graph_cpu_fixture();wg,queued=graph_cpu_fixture();wp,private=graph_cpu_fixture()
    clock=Interior.DeviceGraphStepIndex(queued.environment.metadata)
    for _ in 1:20
        Interior.device_cellular_step!(wr,reference,.001)
        private.attempted_steps+=1
        Interior.device_graph_enqueue_prototype!(wp,private,.001f0,:projected,nothing,UInt32(private.attempted_steps))
        @test Interior.device_cellular_snapshot(reference)==Interior.device_cellular_snapshot(private)
        queued.attempted_steps+=1
        Interior.device_cellular_enqueue!(wg,queued,.001f0,:projected,nothing,clock)
        Interior.device_cellular_status!(queued)
        @test Interior.device_cellular_snapshot(reference)==Interior.device_cellular_snapshot(queued)
    end
    @test reference.environment.metadata[3]>0
    checkpoint=Interior.device_cellular_checkpoint(wg,queued)
    restored_world=World(wg.grid;seed=911);fill!(restored_world.fields.nutrient,0)
    restored=Interior.device_cellular_restore(restored_world,checkpoint)
    @test Interior.device_cellular_snapshot(restored)==Interior.device_cellular_snapshot(queued)
end

isdefined(Interior,:DeviceCellularGraph) || Base.include(Interior,joinpath(@__DIR__,"..","studies","amd_cellular_graph.jl"))
@testset "Graph replay bookkeeping and failure preflight without HIP" begin
    w,p=graph_cpu_fixture();wr,reference=graph_cpu_fixture();calls=Ref(0)
    clock=Interior.DeviceGraphStepIndex(p.environment.metadata)
    launch=(_)->begin
        calls[]+=1
        Interior.device_cellular_enqueue!(w,p,.001f0,:projected,nothing,clock)
    end
    runtime=(HIP=(;launch),)
    owner=Interior.DeviceCellularGraph(w,p,nothing,nothing,nothing,runtime,.001f0,:projected)
    Interior.device_cellular_graph_replay!(owner,3)
    for _ in 1:3;Interior.device_cellular_step!(wr,reference,.001);end
    @test Interior.device_cellular_snapshot(p)==Interior.device_cellular_snapshot(reference)
    @test calls[]==3 && p.attempted_steps==3 && w.time==0
    checkpoint=Interior.device_cellular_checkpoint(w,p)
    replacement_world=World(w.grid;seed=911);fill!(replacement_world.fields.nutrient,0)
    replacement=Interior.device_cellular_restore(replacement_world,checkpoint)
    original=(state=p.state,environment=p.environment,ownership=p.ownership)
    p.state=replacement.state;p.environment=replacement.environment;p.ownership=replacement.ownership
    @test_nowarn Interior.device_cellular_validate(w,p)
    @test_throws ArgumentError Interior.device_cellular_graph_replay!(owner)
    @test calls[]==3
    p.state=original.state;p.environment=original.environment;p.ownership=original.ownership
    p.attempted_steps+=1
    @test_throws ArgumentError Interior.device_cellular_graph_replay!(owner)
    @test calls[]==3
    p.attempted_steps-=1;p.environment.status[1]=8
    @test_throws ArgumentError Interior.device_cellular_graph_replay!(owner)
    @test calls[]==3
    p.environment.status[1]=0;p.environment.metadata[6]=typemax(Int32);p.attempted_steps=typemax(Int32)
    @test_throws ArgumentError Interior.device_cellular_graph_replay!(owner)
    @test calls[]==3
end
@testset "Extinct device audit has an additive identity" begin
    w,p=graph_cpu_fixture();p.state.active.=0;p.environment.metadata[2]=0
    audit=Interior.device_cellular_audit(Interior.device_cellular_snapshot(p))
    @test audit.population==0
    @test isfinite(audit.material) && isfinite(audit.potential)
end

mutable struct GraphHookOwnershipFixture{O}
    ownership::O
end
@testset "Captured habitat binding rejects replacement before dispatch" begin
    w,p=graph_cpu_fixture();hook=GraphHookOwnershipFixture((marker=Ref(1),));calls=Ref(0)
    runtime=(HIP=(launch=(_)->(calls[]+=1),),)
    owner=Interior.DeviceCellularGraph(w,p,hook,nothing,nothing,runtime,.001f0,:projected)
    hook.ownership=(marker=Ref(1),)
    error=try
        Interior.device_cellular_graph_replay!(owner)
        nothing
    catch exception
        exception
    end
    @test error isa ArgumentError && occursin("habitat",error.msg)
    @test calls[]==0
end
