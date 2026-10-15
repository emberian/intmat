using Interior, Test, Random
include(joinpath(@__DIR__,"..","studies","radial_moment_tiles.jl"))
using .RadialMomentTiles

@testset "Hierarchical radial moments retain the physical support" begin
    rng=MersenneTwister(541)
    for (dims,extent,radius) in (((16,16,1),(4.,4.,1.),.8),
            ((16,12,8),(4.,3.,2.),.9),((9,7,3),(4.,3.,2.),10.),
            ((1,8,1),(4.,4.,1.),.01),((1,1,1),(1.,1.,1.),1e-8))
        grid=Grid(dims...;extent);positions=rand(rng,3,24).*collect(extent)
        positions[:,1].=[.001,.001,0.];positions[:,2].=collect(extent).-1e-5
        field=randn(rng,dims...)
        stencil=CellStencil(grid,positions,radius)
        for width in (2,4)
            pyramid=MomentPyramid(field,grid,radius;leaf_width=width)
            for i in axes(positions,2)
                expected=sum(stencil.weights[e]*field[stencil.sites[e]] for e in stencil.offsets[i]:stencil.offsets[i+1]-1)
                actual=sample(pyramid,view(positions,:,i))
                @test isapprox(actual.value,expected;atol=2e-12,rtol=2e-12)
                @test actual.accepted_voxels==stencil.offsets[i+1]-stencil.offsets[i]
            end
            fill!(field,2.)
            refresh!(pyramid)
            @test sample(pyramid,view(positions,:,3)).value≈2 atol=2e-12
            field.=randn(rng,dims...)
        end
    end
end

@testset "Dense co-location does not omit cells or periodic boundary sites" begin
    grid=Grid(32,32,32;extent=(4.,4.,4.));field=rand(MersenneTwister(66),Float64,grid.size...)
    pyramid=MomentPyramid(field,grid,.9;leaf_width=2)
    positions=repeat(reshape([.01,3.99,.01],3,1),1,32)
    stencil=CellStencil(grid,positions,.9)
    for i in axes(positions,2)
        expected=sum(stencil.weights[e]*field[stencil.sites[e]] for e in stencil.offsets[i]:stencil.offsets[i+1]-1)
        actual=sample(pyramid,view(positions,:,i))
        @test actual.value≈expected atol=2e-12 rtol=2e-12
        @test actual.full_nodes>0 && actual.direct_voxels>0
    end
    @test memory_forecast((512,512,512)).coupled_grid_bytes==35_970_351_104
    @test storage_bytes(pyramid)==sum(sizeof,pyramid.levels)
end

@testset "A held support plan is reusable across fields and refreshes" begin
    rng=MersenneTwister(719);grid=Grid(32,24,16;extent=(4.,3.,2.))
    field=rand(rng,Float32,grid.size...);other=randn(rng,grid.size...)
    p=MomentPyramid(field,grid,.8;leaf_width=2);q=MomentPyramid(other,grid,.8;leaf_width=2)
    for _ in 1:20
        position=rand(rng,3).*collect(grid.extent)
        plan=compile_support(p,position)
        @test sample(p,plan)≈sample(p,position).value atol=2e-12 rtol=2e-12
        @test sample(q,plan)≈sample(q,position).value atol=2e-12 rtol=2e-12
        @test storage_bytes(plan)==16length(plan.blocks)+sum(sizeof,plan.axis_distance2)
    end
    plan=compile_support(p,[.1,.1,.1]);fill!(field,3);refresh!(p)
    @test sample(p,plan)≈3 atol=2e-12
    incompatible=MomentPyramid(field,grid,.9;leaf_width=2)
    @test_throws ArgumentError sample(incompatible,plan)
    # Float32 stored moments are an explicitly different arithmetic option.
    # Check its sensitivity against the Float64 oracle, not bitwise identity.
    approximate=MomentPyramid(other,grid,.8;leaf_width=2,precision=Float32)
    for _ in 1:20
        position=rand(rng,3).*collect(grid.extent)
        @test sample(approximate,position).value≈sample(q,position).value atol=2e-5 rtol=2e-5
    end
end

@testset "Polynomial request scatter is the gather transpose" begin
    rng=MersenneTwister(922);grid=Grid(24,24,24;extent=(4.,4.,4.))
    field=randn(rng,grid.size...);p=MomentPyramid(field,grid,.9;leaf_width=2)
    positions=rand(rng,3,16).*4;positions[:,1].=[.001,3.999,.001]
    positions[:,2:5].=positions[:,1] # arbitrary local crowding keeps all contributors
    plans=[compile_support(p,view(positions,:,i)) for i in axes(positions,2)]
    amounts=rand(rng,size(positions,2));stencil=CellStencil(grid,positions,.9)
    expected=zeros(grid.size)
    for i in axes(positions,2),e in stencil.offsets[i]:stencil.offsets[i+1]-1
        expected[stencil.sites[e]]+=amounts[i]*stencil.weights[e]
    end
    actual=scatter(p,plans,amounts)
    @test all(>=(0),actual.field)
    @test maximum(abs.(actual.field.-expected))<2e-12
    @test sum(actual.field)≈sum(amounts) atol=2e-12 rtol=2e-12
    @test sum(field.*actual.field)≈sum(amounts[i]*sample(p,plans[i]) for i in eachindex(plans)) atol=2e-12 rtol=2e-12
    @test all(iszero,scatter(p,plans,zeros(length(plans))).field)
    @test_throws ArgumentError scatter(p,plans,fill(-1.,length(plans)))
end

@testset "Synchronous request proof funds actual grants from field withdrawal" begin
    rng=MersenneTwister(233);grid=Grid(16,16,16;extent=(4.,4.,4.))
    field=fill(.125f0,grid.size);p=MomentPyramid(field,grid,.9;leaf_width=2)
    positions=rand(rng,3,16).*4;positions[:,2:6].=positions[:,1]
    plans=[compile_support(p,view(positions,:,i)) for i in axes(positions,2)]
    requests=rand(rng,16);initial=copy(field)
    result=allocate_requests(field,p,plans,requests;potential=4.)
    @test field==initial && result.field!==field
    @test all(>=(0),result.field) && all(>=(0),result.grants) && result.inert_residual>=0
    @test all(result.grants.<=requests)
    removed=sum(BigFloat(initial[i])-BigFloat(result.field[i]) for i in eachindex(field))
    received=sum(BigFloat,result.grants)+BigFloat(result.inert_residual)
    @test abs(removed-received)<=64BigFloat(eps(result.withdrawn))
    @test abs(4removed-(4sum(BigFloat,result.grants)+BigFloat(result.chemical_heat)))<=256BigFloat(eps(result.withdrawn))
    empty=allocate_requests(zeros(Float32,grid.size),p,plans,requests)
    @test all(iszero,empty.grants) && empty.withdrawn==empty.inert_residual==0
end

@testset "Lazy moments build only shared queried interior subtrees" begin
    rng=MersenneTwister(532);grid=Grid(48,40,32;extent=(4.,4.,4.))
    field=randn(rng,grid.size...);geometry=MomentGeometry(grid,.55;leaf_width=2)
    positions=rand(rng,3,12).*4;positions[:,2:5].=positions[:,1]
    plans=[compile_support(geometry,view(positions,:,i)) for i in axes(positions,2)]
    dense=MomentPyramid(field,grid,.55;leaf_width=2)
    lazy=LazyMomentPyramid(field,geometry,plans)
    @test storage_bytes(lazy)<storage_bytes(dense)
    @test [sample(lazy,p) for p in plans]==[sample(dense,p) for p in plans]
    repeated=LazyMomentPyramid(field,geometry,repeat(plans,3))
    @test repeated.index==lazy.index && repeated.moments==lazy.moments
    field .*=2;refresh!(dense);refresh!(lazy)
    @test [sample(lazy,p) for p in plans]==[sample(dense,p) for p in plans]
    missing=findfirst(==(0),lazy.index[1])
    @test missing!==nothing
    @test_throws ArgumentError RadialMomentTiles.node_moments(lazy,1,Tuple(missing)...)
    @test_throws ArgumentError LazyMomentPyramid(field,MomentGeometry(grid,.56),plans)
end
