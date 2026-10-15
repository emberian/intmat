using Test,Interior,Random
include(joinpath(@__DIR__,"..","studies","physical_support_prototype.jl"))
using .PhysicalSupportPrototype
@testset "Axis-factored physical supports" begin
    rng=MersenneTwister(41)
    for dims in ((16,16,16),(12,10,1)),radius in (.03,.6,3.)
        grid=Grid(dims...;extent=(4,4,4));positions=4rand(rng,3,12)
        reference=CellStencil(grid,positions,radius);candidate=AxisSupport(grid,positions,radius)
        pool=ResourcePool(World(grid;fields=());initial=(.3,.2,.1))
        for field in pool.fields;field .*= rand(rng,Float32,dims...);end
        concentrations=zeros(3,12);gradients=zeros(3,3,12)
        sample_support!(concentrations,gradients,pool.fields,candidate)
        @test concentrations≈sample_concentrations(pool,reference) rtol=2e-13 atol=1e-12
        @test gradients≈sample_concentration_gradients(pool,reference) rtol=2e-13 atol=1e-12
        for cell in 1:12
            actual=Dict{Int,Float64}()
            foreach_support_edge(candidate,cell) do site,weight,i,j,k;actual[site]=weight;end
            expected=Dict(reference.sites[e]=>reference.weights[e] for e in reference.offsets[cell]:reference.offsets[cell+1]-1)
            @test keys(actual)==keys(expected)
            @test sum(values(actual))≈1 atol=2e-14
            @test all(isapprox(actual[site],weight;rtol=2e-12,atol=1e-15) for (site,weight) in expected)
        end
        # Complete fair withdrawal, not just interpolation on an inert field.
        other=deepcopy(pool);requests=rand(rng,3,12)*.4
        grants=allocate_cell_resources!(pool,reference,requests)
        received=allocate_support!(other,candidate,requests,SupportAllocationScratch(candidate))
        @test grants≈received rtol=2e-7 atol=2e-8
        @test all(maximum(abs.(a.-b))<=4eps(Float32) for (a,b) in zip(pool.fields,other.fields))
        @test sum(received)≈sum(other.withdrawn) rtol=1e-12 atol=1e-12
        @test all(x->x>=0,Iterators.flatten(other.fields))
    end
end
@testset "Float32 descriptor vs existing KA device support" begin
    isdefined(Interior,:DeviceCellular) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_cellular.jl"))
    rng=MersenneTwister(31)
    cases=vcat([(dims,radius) for dims in ((16,16,16),(12,10,1)) for radius in (.03f0,.6f0,3f0)],
               [((n,n,n),.18f0) for n in (64,128,256)])
    for (dims,radius) in cases
        grid=Grid(dims...;extent=(4,4,4));n=11;positions=4f0.*rand(rng,Float32,3,n)
        nd=dims[3]==1 ? 2 : 3
        counts=ntuple(d->d>nd ? 1 : min(dims[d],2ceil(Int,radius/grid.spacing[d])+2),3)
        slots=prod(counts)
        state=(support_count=zeros(Int32,n),active=ones(Int32,n),position=permutedims(positions),
               sites=zeros(Int32,n,slots),support=zeros(Float32,n,slots))
        environment=(status=zeros(Int32,1),)
        parameters=(exchange_radius=radius,)
        geometry=(nd=Int32(nd),dims=Int32.(dims),spacing=grid.spacing,extent=grid.extent,candidates=Int32.(counts))
        Interior.dc_support!(CPU(),128)(state,environment,parameters,geometry;ndrange=n)
        candidate=AxisSupport(grid,positions,radius;precision=Float32)
        for cell in 1:n
            actual=Tuple{Int,Float32}[]
            foreach_support_edge(candidate,cell) do site,weight,i,j,k;push!(actual,(site,weight));end
            expected=[(Int(state.sites[cell,e]),state.support[cell,e]) for e in 1:state.support_count[cell]]
            @test actual==expected
        end
    end
end
@testset "Support input failures precede exchange" begin
    grid=Grid(16,16,16;extent=(4,4,4));x=fill(.3,3,2)
    support=AxisSupport(grid,x,.18);scratch=SupportAllocationScratch(support)
    wrong=ResourcePool(World(Grid(8,8,8);fields=()))
    before=deepcopy(wrong.fields)
    @test_throws DimensionMismatch allocate_support!(wrong,support,zeros(3,2),scratch)
    @test wrong.fields==before && all(iszero,wrong.withdrawn)
    @test_throws ArgumentError AxisSupport(grid,x,1e-100;precision=Float32)
    pool=ResourcePool(World(grid;fields=()))
    @test_throws ArgumentError allocate_support!(pool,support,fill(Inf,3,2),scratch)
end
