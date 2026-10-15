using Interior,Test,KernelAbstractions,Random
isdefined(Interior,:DeviceLifecyclePrefixCarrier) || Base.include(Interior,joinpath(@__DIR__,"..","studies","device_lifecycle_prefix.jl"))
isdefined(@__MODULE__,:constructive_device_fixture) || include(joinpath(@__DIR__,"..","scripts","device_habitat_fixture.jl"))
function lifecycle_fixture(C,active,eligible;nextid=100)
    s=(active=Int32.(active),ids=Int32.(1:C),mass=fill(1f0,C),actions=zeros(Float32,C,Interior.CELLULAR_ACTIONS),new_parent=zeros(Int32,C))
    s.mass[eligible].=3f0;s.actions[eligible,10].=1f0
    e=(status=zeros(Int32,1),free=zeros(Int32,C),heap=zeros(Int32,C),heap_keys=zeros(UInt32,C),metadata=Int32[nextid,0,0,0,0,0,0,0,0,0])
    s,e,(division_mass=2f0,seed=UInt32(123)),(capacity=Int32(C),)
end
@testset "Prefix lifecycle preserves all used free slots and birth decisions" begin
    rng=MersenneTwister(643)
    for C in (1,7,64,16384),H in unique((0,min(3,C),C)),trial in 1:8
        active=zeros(Int32,C);H>0 && (active[1:H].=rand(rng,Int32[0,1],H))
        eligible=[i for i in 1:H if active[i]==1 && rand(rng,Bool)]
        s,e,q,g=lifecycle_fixture(C,active,eligible;nextid=C+100);s.ids[H+1:end].=0
        s2=deepcopy(s);e2=deepcopy(e);h=Int32[H]
        Interior.dc_lifecycle!(CPU(),32)(s,e,q,g,UInt32(trial);ndrange=1)
        Interior.dc_lifecycle_prefix!(CPU(),32)(s2,e2,q,g,UInt32(trial),h;ndrange=1)
        KernelAbstractions.synchronize(CPU())
        births=Int(e.metadata[3])
        @test s==s2
        @test e.metadata==e2.metadata && e.status==e2.status && e.heap==e2.heap && e.heap_keys==e2.heap_keys
        @test e.free[1:births]==e2.free[1:births]
        @test h[1]==max(H,births>0 ? maximum(e.free[1:births]) : 0)
    end
end
function lifecycle_canonical(snapshot)
    env=snapshot.cells.environment;names=filter(!=(:free),keys(env))
    merge(snapshot,(cells=merge(snapshot.cells,(environment=NamedTuple{names}(env),)),))
end
if get(ENV,"INTMAT_PREFIX_DIRECT_ONLY","0")!="1"
@testset "Complete coupled prefix trajectory, lifecycle and latch" begin
    w,p=constructive_device_fixture();cp=Interior.device_habitat_checkpoint(w,p)
    w2=World(w.grid;seed=911);fill!(w2.fields.nutrient,0)
    candidate=Interior.device_habitat_restore(w2,cp)
    carrier=Interior.DeviceLifecyclePrefixCarrier(candidate.cells)
    @test only(carrier.highwater)==2
    for step in 1:20
        Interior.device_cellular_step!(w,p.cells,.01;anabolism=p.habitat,lifecycle=:dense)
        Interior.device_lifecycle_prefix_step!(w2,carrier,.01;anabolism=candidate.habitat)
        @test lifecycle_canonical(Interior.device_habitat_snapshot(p))==lifecycle_canonical(Interior.device_habitat_snapshot(candidate))
    end
    @test p.cells.environment.metadata[3]>0
    restored=Interior.DeviceLifecyclePrefixCarrier(candidate.cells)
    @test only(restored.highwater)<=only(carrier.highwater)
    candidate.cells.environment.status[1]=8
    before=Interior.device_habitat_snapshot(candidate;check_status=false);h=copy(carrier.highwater)
    @test_throws ArgumentError Interior.device_lifecycle_prefix_step!(w2,carrier,.01;anabolism=candidate.habitat)
    after=Interior.device_habitat_snapshot(candidate;check_status=false)
    @test before.cells.state==after.cells.state && before.arrays==after.arrays && before.feedstock==after.feedstock
    @test before.cells.environment==after.cells.environment && h==carrier.highwater
end

end # optional complete coupled gate

@testset "Prefix reconstruction and failure boundary" begin
    fake=(state=(ids=Int32[1,0,4,0],active=zeros(Int32,4),new_parent=zeros(Int32,4)),backend=CPU())
    @test Interior.device_lifecycle_host_highwater(fake.state)==3
    fake.state.new_parent[4]=1
    @test Interior.device_lifecycle_host_highwater(fake.state)==4
    fill!(fake.state.ids,0);fill!(fake.state.new_parent,0)
    @test Interior.device_lifecycle_host_highwater(fake.state)==0
    s,e,q,g=lifecycle_fixture(16,Int32[1,1,zeros(Int32,14)...],[1,2];nextid=typemax(Int32)-1)
    h=Int32[2]
    Interior.dc_lifecycle_prefix!(CPU(),32)(s,e,q,g,UInt32(1),h;ndrange=1)
    KernelAbstractions.synchronize(CPU())
    @test e.status==[4] && h==[2] && all(iszero,s.new_parent)
    before=deepcopy((s,e,h))
    Interior.dc_lifecycle_prefix!(CPU(),32)(s,e,q,g,UInt32(2),h;ndrange=1)
    KernelAbstractions.synchronize(CPU())
    @test before==(s,e,h)
end
