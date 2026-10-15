using Interior,Test,KernelAbstractions
isdefined(@__MODULE__,:constructive_device_fixture) || include(joinpath(@__DIR__,"..","scripts","device_habitat_fixture.jl"))
function lifecycle_habitat_alias(p,role)
    old=p.cells;h=p.habitat
    target=role===:enzyme ? h.arrays.enzyme_bits : role===:bound ? h.arrays.bound : h.feedstock.arrays.fields
    # A valid H=2 and finite nonnegative physical value share these same bits.
    # Bound inventory is given genuine enzyme capacity, so no other physical
    # preflight can accidentally mask this coherent ownership counterexample.
    role===:bound && (h.arrays.enzyme_bits[1,1]=reinterpret(UInt32,1f0))
    target[1]=reinterpret(eltype(target),Int32(2))
    highwater=view(reinterpret(Int32,vec(target)),1:1)
    work=Interior.DeviceLifecycleWorkspace(old.backend,Int(old.geometry.capacity);highwater)
    ownership=Interior.device_cellular_contract(old.state,old.environment,old.parameters,old.geometry;
        resource_workspace=Interior.device_resource_workspace(old),lifecycle_workspace=work)
    c=Interior.DeviceCellular(old.state,old.environment,old.backend,old.parameters,old.geometry,old.attempted_steps,ownership)
    argument=Interior.device_habitat_argument(c,h.arrays,h.parameters)
    contract=(cells=c,arrays=h.arrays,feedstock=h.feedstock,parameters=h.parameters,argument,
        array_shapes=map(size,h.arrays),feed_shapes=map(size,h.feedstock.arrays))
    habitat=Interior.DeviceCatalyticHabitat(c,h.arrays,h.feedstock,h.parameters,argument,contract)
    Interior.DeviceConstructiveEcology(c,habitat)
end
@testset "Lifecycle watermark cannot alias the complete habitat inventory" begin
    for role in (:enzyme,:bound,:feed)
        w,p=constructive_device_fixture();bad=lifecycle_habitat_alias(p,role)
        @test Interior.device_cellular_validate(w,bad.cells)
        @test only(Interior.device_lifecycle_workspace(bad.cells).highwater)==2
        before=Interior.device_habitat_snapshot(bad)
        fill!(w.fields.activity,1)
        fields=deepcopy(w.fields);clock=w.time;attempts=bad.cells.attempted_steps
        @test_throws ArgumentError Interior.validate(w,bad)
        @test_throws ArgumentError step!(w,(Deposition(threshold=0),bad),.001)
        @test fields==w.fields && w.time==clock
        @test bad.cells.attempted_steps==attempts && Interior.device_habitat_snapshot(bad)==before
        @test_throws ArgumentError Interior.device_habitat_checkpoint(w,bad)
        @test bad.cells.attempted_steps==attempts && Interior.device_habitat_snapshot(bad)==before
    end
end
