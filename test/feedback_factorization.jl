using Interior, Test, KernelAbstractions
include(joinpath(@__DIR__,"..","studies","feedback_factorization.jl"))
include(joinpath(@__DIR__,"..","scripts","device_habitat_fixture.jl"))
@testset "Exact fused feedback complete steps" begin
    for strategy in (:materialized,:axis), channels in ((:raw_total,:bound_total,:catalyst_total),(:raw_balance,:bound_a,:catalyst_balance))
        sensing=(;channels,scales=(.6,.02,.004))
        w,a=constructive_device_fixture(;support_strategy=strategy,inventory_type=Float64,energy_type=Float64,material_sensing=sensing)
        v,b=constructive_device_fixture(;support_strategy=strategy,inventory_type=Float64,energy_type=Float64,material_sensing=sensing)
        for _ in 1:12
            advance!(w,a,.01);feedback_factorized_step!(v,b,.01)
            @test Interior.device_habitat_snapshot(a)==Interior.device_habitat_snapshot(b)
        end
    end
end

@testset "Feedback fusion boundary and disabled guards" begin
    w,a=constructive_device_fixture()
    @test_throws ErrorException feedback_factorized_step!(w,a,.01)
    sense=(channels=(:raw_a,:bound_balance,:catalyst_b),scales=(.6,.02,.004))
    w,a=constructive_device_fixture(;material_sensing=sense)
    v,b=constructive_device_fixture(;material_sensing=sense)
    for p in (a,b)
        p.cells.state.position[1,1]=.001f0
        p.cells.state.position[2,2]=7.999f0
    end
    for _ in 1:3
        advance!(w,a,.01);feedback_factorized_step!(v,b,.01)
        @test Interior.device_habitat_snapshot(a)==Interior.device_habitat_snapshot(b)
    end
end

install_zero_diffusion_factorization!()
@testset "Zero diffusion preserves complete scratch and physical state" begin
    for inventory_type in (Float32,Float64)
        w,a=constructive_device_fixture(;inventory_type,energy_type=inventory_type)
        v,b=constructive_device_fixture(;inventory_type,energy_type=inventory_type)
        for _ in 1:5
            advance!(w,a,.01);feedback_reference_diffusion_step!(v,b,.01)
            @test Interior.device_habitat_snapshot(a)==Interior.device_habitat_snapshot(b)
        end
    end
end

@testset "All-zero diffusion frozen-field oracle" begin
    for T in (Float32,Float64)
        fields=reshape(T.(1:48),16,3);fields[1,1]=-zero(T)
        e=(fields=fields,scratch=similar(fields),permeability=fill(.4f0,16),status=Int32[0])
        f=merge(e,(scratch=similar(fields),))
        g=(dims=(Int32(4),Int32(4),Int32(1)),spacing=(.2f0,.2f0,1f0))
        for diffusion in ((0f0,0f0,0f0),(.01f0,0f0,.03f0))
            q=(;diffusion)
            Interior.dc_diffuse!(CPU(),64)(e,q,g,.01f0;ndrange=(16,3))
            Interior.ff_reference_diffuse!(CPU(),64)(f,q,g,.01f0;ndrange=(16,3))
            @test reinterpret(UInt8,vec(e.scratch))==reinterpret(UInt8,vec(f.scratch))
        end
    end
end
