isdefined(Interior,:install_cellular_capped_reference!)||Base.include(Interior,joinpath(@__DIR__,"..","studies","cellular_capped_reference.jl"))
@testset "Capped production reference keeps the outer budget horizon" begin
    reference=Interior.install_cellular_capped_reference!().method
    grid=Grid(8,8;extent=(4.,4.,1.));ids=[1,2];radii=fill(.5,2);original=[1. 1.6;1. 1.;0. 0.]
    x=copy(original);y=copy(x);a=CellBonds();b=CellBonds();funds=fill(.001,2);other=copy(funds)
    kwargs=(;adhesion=zeros(2),motility=[.2 -.2;0. 0.;0. 0.])
    default=cellular_mechanics!(grid,ids,x,radii,a,.02;energy=funds,kwargs...)
    uncapped=reference(grid,ids,y,radii,b,.02;energy=other,kwargs...)
    @test isequal(x,y)&&isequal(funds,other)&&isequal(default,uncapped)
    errors=Float64[]
    for h in (.002,.001,.0005)
        y=copy(original)
        r=reference(grid,ids,y,radii,CellBonds(),.02;max_step=h,adhesion=zeros(2))
        @test r.substeps>=round(Int,.02/h)
        @test abs(r.energy_residual)<1e-12
        push!(errors,abs(1-(y[1,2]-y[1,1])-.4exp(-.4)))
    end
    @test errors[1]>errors[2]>errors[3]
    @test errors[2]/errors[3]>1.9
end
