using Interior,Test
include(joinpath(@__DIR__,"../studies/device_memory_specialization.jl"))
include(joinpath(@__DIR__,"../scripts/device_habitat_fixture.jl"))
using .DeviceMemorySpecialization
@testset "exact specialized coupled storage" begin
    for pages in (false,true)
        w,p=constructive_device_fixture(CPU();support_strategy=:axis)
        w2,p2=deepcopy((w,p));special=DeviceMemorySpecialization.specialize(w2,p2;pages)
        @test DeviceMemorySpecialization.retained_array_bytes(special)<DeviceMemorySpecialization.retained_array_bytes(p)
        for i in 1:12
            step!(w,p,.01);step!(w2,special,.01)
            @test isequal(DeviceMemorySpecialization.persistent_snapshot(p),DeviceMemorySpecialization.persistent_snapshot(special))
        end
    end
end
@testset "storage invariants and portable continuation" begin
    M=DeviceMemorySpecialization
    a=zeros(Float32,7,3);a[:,1].=1:7
    masks=M.mask_columns(a)
    @test isequal(Array(masks),a)
    @test masks.columns[1] isa Vector && masks.columns[2] isa M.UniformColumn
    @test_throws ArgumentError M.zero_tail(ones(Float32,7,3),2)
    z=M.zero_tail(zeros(Float32,7,3),2)
    @test_throws ArgumentError (z[1,3]=1f0)
    @test_throws ArgumentError (z[1,3]=-0f0)
    sparse=M.paged(zeros(Float32,23,2);page_rows=4)
    @test isempty(sparse.pages)
    sparse[2,1]=1f-20;sparse[23,2]=-0f0
    @test sparse[2,1]===1f-20 && isequal(sparse[23,2],-0f0)
    @test length(sparse.pages)==2
    for pages in (false,true)
        w,p=constructive_device_fixture(CPU();support_strategy=:axis)
        p=M.specialize(w,p;pages)
        for _ in 1:5;step!(w,p,.01);end
        checkpoint=Interior.device_habitat_checkpoint(w,p)
        world=World(w.grid,CPU(),NamedTuple(),NamedTuple(),w.time,w.seed)
        restored=Interior.device_habitat_restore(world,checkpoint)
        @test isequal(M.persistent_snapshot(p),M.persistent_snapshot(restored))
        for _ in 1:5
            step!(w,p,.01);step!(world,restored,.01)
            @test isequal(M.persistent_snapshot(p),M.persistent_snapshot(restored))
        end
    end
end
