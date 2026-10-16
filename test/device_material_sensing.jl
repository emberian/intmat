using Interior, Test, KernelAbstractions
isdefined(@__MODULE__,:constructive_device_fixture) || include(joinpath(@__DIR__,"..","scripts","device_habitat_fixture.jl"))
@testset "Explicit material receptor vocabulary" begin
    param=Interior.device_material_sensing_parameters
    @test param((channels=(:raw_a,:bound_total,:catalyst_balance),scales=(2,.02,.004))).material_channels==(1,7,12)
    @test_throws ArgumentError param((channels=(:raw_a,:bound_total,:missing),scales=(2,.02,.004)))
    @test_throws ArgumentError param((channels=(:raw_a,:bound_total,:catalyst_a),scales=(2,0,.004)))
    @test_throws ArgumentError param((channels=(:raw_a,:bound_total,:catalyst_a),scales=(2,Inf,.004)))
    @test_throws ArgumentError param((channels=(:raw_a,),scales=(2,)))
    # Same physical concentrations in two different voxel volumes. Two-edge
    # support average is taken before saturation, and all other slots survive.
    for volume in (1f0,.125f0)
        inputs=fill(9f0,2,24);s=(active=Int32[1,0],inputs=inputs)
        descriptor=(sites=Int32[1 2;1 2],weights=Float32[.25 .75;.25 .75],count=Int32[2,2])
        raw=Float32[2 6;6 2].*volume;bound=Float32[.01 .03;.03 .01].*volume
        enzyme=Float32[.001 .003;.003 .001].*volume
        field=(raw=raw,bound=bound,enzyme_bits=copy(reinterpret(UInt32,enzyme)))
        before=deepcopy(field);status=Int32[0];g=(volume=volume,hidden=Int32(2))
        channels=Int32.(Tuple(1:3))
        Interior.dh_material_sense!(CPU(),64)(s,field,descriptor,status,channels,(5f0,3f0,8f0),g;ndrange=2)
        @test inputs[1,20:22]≈fill(.5f0,3)
        @test inputs[2,:]==fill(9f0,24)
        @test inputs[1,[1:19;23:24]]==fill(9f0,21)
        @test field==before
        Interior.dh_material_sense!(CPU(),64)(s,field,descriptor,status,(Int32(4),Int32(7),Int32(11)),(2f0,.04f0,.004f0),g;ndrange=2)
        @test inputs[1,20:22]≈fill(.5f0,3)
        status[1]=8;fill!(inputs,9)
        Interior.dh_material_sense!(CPU(),64)(s,field,descriptor,status,channels,(1f0,1f0,1f0),g;ndrange=2)
        @test all(==(9),inputs)
    end
end
@testset "Opt-in sensing composes and strictly restarts" begin
    sense=(channels=(:raw_total,:bound_total,:catalyst_total),scales=(.6,.02,.004))
    w,off=constructive_device_fixture(;energy_type=Float64,inventory_type=Float64)
    w2,on=constructive_device_fixture(;energy_type=Float64,inventory_type=Float64,material_sensing=sense)
    @test !hasproperty(off.habitat.parameters,:material_channels)
    @test off.cells.state==on.cells.state
    initial=Interior.device_habitat_audit(Interior.device_habitat_snapshot(on))
    for _ in 1:5;advance!(w,off,.01);advance!(w2,on,.01);end
    h=Int(on.cells.geometry.hidden);active=findall(==(1),on.cells.state.active)
    @test all(==(0),off.cells.state.inputs[active,2h+16:2h+18])
    @test all(>(0),on.cells.state.inputs[active,2h+16])
    @test maximum(abs,on.cells.state.hidden.-off.cells.state.hidden)>1e-5
    final=Interior.device_habitat_audit(Interior.device_habitat_snapshot(on))
    @test abs(final.material-initial.material)<1e-8
    @test abs(final.potential+final.work+final.heat-initial.potential)<1e-8
    cp=Interior.device_habitat_checkpoint(w2,on)
    rw=World(w2.grid;seed=911);fill!(rw.fields.nutrient,0)
    restored=Interior.device_habitat_restore(rw,cp)
    @test Interior.device_habitat_snapshot(restored)==Interior.device_habitat_snapshot(on)
    for _ in 1:3;advance!(w2,on,.01);advance!(rw,restored,.01);end
    @test Interior.device_habitat_snapshot(restored)==Interior.device_habitat_snapshot(on)
    corrupted=merge(cp,(parameters=merge(cp.parameters,(material_scales=(1f0,0f0,1f0),)),))
    @test_throws ArgumentError Interior.device_habitat_restore(rw,corrupted)
end
