using Interior, Test, KernelAbstractions
Base.include(Interior,joinpath(@__DIR__,"fixtures","cellular_initialization_reference.jl"))
@testset "Shared cellular founder construction preserves all fields" begin
    for h in (2,8),n in (1,7),seed in (0,911)
        w=World(Grid(8,8;extent=(8,8,1));seed)
        r=ResourcePool(w;potentials=(4,4,0))
        positions=reshape(collect(1.:3n),3,n)
        current=CellularEcology(w,positions,r;hidden=h,seed)
        reference=Interior.cellular_initialization_reference(w,positions,r;hidden=h,seed)
        for field in fieldnames(typeof(current))
            a=getfield(current,field);b=getfield(reference,field)
            if a isa AbstractArray
                @test isequal(a,b)
            elseif a isa Union{Number,Nothing,NamedTuple}
                @test isequal(a,b)
            end
        end
        @test current.resources===reference.resources
        @test current.parameters===Interior.cellular_parameters()
    end
    @test_throws ArgumentError Interior.cellular_parameters(efficiency=2.)
    @test_throws ArgumentError Interior.cellular_parameters(minimum_mass=1.)
    @test_throws ArgumentError Interior.cellular_validate_founder_inventory(Interior.cellular_parameters(),1.,1.,1.)
    @test_throws ArgumentError Interior.cellular_founder_weights(0,1,2)
end
