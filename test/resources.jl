@testset "finite resource chemistry" begin
    w=World(Grid(8,7;extent=(8,7,1)))
    p=ResourcePool(w;initial=(0,0,0))
    p.fields[1][1,1,1]=1
    mass=resource_mass(p); energy=resource_potential(p)
    step!(w,p,10)
    @test all(a->minimum(a)>=0,p.fields)
    @test resource_mass(p)≈mass atol=2e-6
    @test resource_potential(p)≈energy atol=1e-5
    cell=CartesianIndex(2,2,1)
    p.fields[1][cell]=1
    grants=allocate_resources!(p,[cell,cell],[2.0 1.0; 0 0; 0 0])
    @test grants[1,:]≈[2/3,1/3]
    @test p.fields[1][cell]==0
    @test sum(grants)==p.withdrawn[1]
    @test_throws ArgumentError allocate_resources!(p,[cell],fill(-1.,3,1))
    resource_add!(p,1,cell,1)
    m=resource_mass(p); e=resource_potential(p)
    result=resource_convert!(p,1,3,cell,1;efficiency=.5)
    @test resource_mass(p)≈m atol=1e-6
    @test result.energy+result.heat≈e-resource_potential(p)
    @test result.energy≈1.5
    @test_throws ArgumentError resource_convert!(p,3,1,cell,1)
    @test resource_convert!(p,1,3,cell,1).energy==0
    q=ResourcePool(w;initial=(0,0,0),diffusion=(0,0,0),source_rates=(1,0,0),source_budgets=(2,0,0),sink_rates=(0,1,0),sink_budgets=(0,1,0))
    q.fields[2].=1
    before=resource_mass(q)
    for _ in 1:20; step!(w,q,1); end
    @test q.supplied[1]<=2
    @test q.removed[2]<=1
    @test q.source_remaining[1]>=0
    @test q.sink_remaining[2]>=0
    @test resource_mass(q)≈before+sum(q.supplied)-sum(q.removed) atol=1e-6
    masks=[zeros(Float32,w.grid.size) for _ in 1:3]
    masks[1][1]=1
    localfeed=ResourcePool(w;initial=(0,0,0),diffusion=(0,0,0),source_rates=(1,0,0),source_budgets=(.25,0,0),source_masks=masks)
    step!(w,localfeed,1)
    @test localfeed.fields[1][1]==.25f0
    @test count(!iszero,localfeed.fields[1])==1
    @test localfeed.supplied[1]==.25
    @test Interior.resource_statistics(localfeed).supplied_potential==1
    @test_throws ArgumentError ResourcePool(w;potentials=(Inf,3,1))
    phase=ImplicitPhaseSeparation(w)
    medium=ActiveMedium(w,phase)
    medium.permeability.=0
    blocked=ResourcePool(w;medium,initial=(0,0,0))
    blocked.fields[1][1]=1
    step!(w,blocked,1)
    @test blocked.fields[1][1]==1
    @test resource_mass(blocked)==1
    @test hasproperty(Interior.record_process_fields(blocked),:statistics)
end

@testset "Resource scalar ledgers preserve aliased storage" begin
    w=World(Grid(2,2;extent=(2,2,1)))
    p=ResourcePool(w;initial=(1,1,0),diffusion=(0,0,0),
        source_rates=(.125,.125,0),sink_rates=(.0625,.0625,0),
        source_budgets=(100,100,0),sink_budgets=(100,100,0))
    aliased=deepcopy(p);aliased.removed=aliased.supplied
    @test Interior._resource_ledgers_independent(p)
    @test !Interior._resource_ledgers_independent(aliased)
    advance!(w,p,1f0);advance!(w,aliased,1f0)
    @test isequal(p.fields,aliased.fields)
    @test aliased.supplied===aliased.removed
    @test aliased.supplied==p.supplied+p.removed==[.75,.75,0.]
    @test isequal(p.source_remaining,aliased.source_remaining)
    @test isequal(p.sink_remaining,aliased.sink_remaining)
end
