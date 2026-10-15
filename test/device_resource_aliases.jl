using Interior,Test,Random,KernelAbstractions
isdefined(@__MODULE__,:device_habitat_pilot_fixture) || include(joinpath(@__DIR__,"..","scripts","device_habitat_pilot_fixture.jl"))
function resource_alias_wrapper(p;arrays=p.habitat.arrays,feed_arrays=p.habitat.feedstock.arrays)
    c=p.cells;old=p.habitat;q=old.parameters
    feed=(arrays=feed_arrays,parameters=old.feedstock.parameters,argument=merge(feed_arrays,(status=c.environment.status,)))
    argument=Interior.device_habitat_argument(c,arrays,q)
    ownership=(cells=c,arrays,feedstock=feed,parameters=q,argument,array_shapes=map(size,arrays),feed_shapes=map(size,feed_arrays))
    h=Interior.DeviceCatalyticHabitat(c,arrays,feed,q,argument,ownership)
    Interior.DeviceConstructiveEcology(c,h)
end
@testset "Resource scratch is disjoint from the whole constructive inventory" begin
    w,p=device_habitat_pilot_fixture(CPU();size=8,founders=2,capacity=32,support_strategy=:materialized)
    grants=p.cells.ownership.stages.resource.grants
    @test Interior.validate(w,p)===nothing
    badfeed=resource_alias_wrapper(p;feed_arrays=merge(p.habitat.feedstock.arrays,(fields=view(grants,:,1:3),)))
    @test_throws ArgumentError Interior.validate(w,badfeed)
    badcatalyst=resource_alias_wrapper(p;arrays=merge(p.habitat.arrays,(enzyme_bits=reinterpret(UInt32,view(grants,:,1:2)),)))
    @test_throws ArgumentError Interior.validate(w,badcatalyst)
    badbound=resource_alias_wrapper(p;arrays=merge(p.habitat.arrays,(bound=view(grants,:,1:2),)))
    @test_throws ArgumentError Interior.validate(w,badbound)
    leased=resource_alias_wrapper(p;feed_arrays=merge(p.habitat.feedstock.arrays,(scratch=view(grants,:,1:3),)))
    @test Interior.validate(w,leased)===nothing
    generic_world=World(w.grid;seed=911);fill!(generic_world.fields.nutrient,0)
    before=copy(generic_world.fields.deposit)
    @test_throws ArgumentError step!(generic_world,(Deposition(threshold=0),badfeed),.001)
    @test generic_world.fields.deposit==before && generic_world.time==0
end
