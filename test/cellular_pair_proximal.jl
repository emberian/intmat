isdefined(Interior,:cellular_implicit_mechanics!)||Base.include(Interior,joinpath(@__DIR__,"..","studies","cellular_implicit_mechanics.jl"))
isdefined(Interior,:cellular_pair_relax)||Base.include(Interior,joinpath(@__DIR__,"..","studies","cellular_pair_proximal.jl"))
@testset "Symmetric radial soft-contact flows" begin
    grid=Grid(8,8;extent=(4.,4.,1.));ids=[1,2];radii=fill(.5,2);bonds=CellBonds()
    x=[1. 1.6;1. 1.;0. 0.]
    for h in (.01,.1,1.)
        result=Interior.cellular_pair_relax(grid,ids,x,radii,bonds,h;repulsion=10.)
        @test result.converged
        @test 1-(result.positions[1,2]-result.positions[1,1])≈.4exp(-20h) atol=1e-14
        @test result.potential_after<=result.potential_before
        @test sum(result.positions;dims=2)≈sum(x;dims=2) atol=1e-14
    end
    seam=[.1 3.9;1. 1.;0. 0.]
    result=Interior.cellular_pair_relax(grid,ids,seam,fill(.2,2),bonds,.1)
    @test result.converged
    @test Interior._cellular_length(Interior._cellular_delta(grid,result.positions,1,2))≈.4-.2exp(-2) atol=1e-14
    spring=CellBonds();add_cell_bond!(spring,1,2,.4)
    result=Interior.cellular_pair_relax(grid,ids,x,fill(.1,2),spring,.1;stiffness=3.)
    @test result.converged
    @test result.positions[1,2]-result.positions[1,1]≈.4+.2exp(-.6) atol=1e-14
    # Existing motor funding and global endpoint gate remain mandatory.
    funds=fill(.001,2);before=copy(funds);y=copy(x)
    d=Interior.cellular_implicit_mechanics!(grid,ids,y,radii,bonds,.1;
        relaxer=Interior.cellular_pair_relax,energy=funds,motility=[1. -1.;0. 0.;0. 0.])
    @test funds≈before-d.spent_motor_work atol=1e-15
    @test abs(d.energy_residual)<1e-12
    @test d.relaxation_energy_loss>=0
end
@testset "Radial splitting consistency on a fixed contact graph" begin
    grid=Grid(12,12;extent=(6.,6.,1.));ids=[1,2,3];radii=fill(.5,3)
    original=[1. 1.3 1.6;1. 1. 1.;0. 0. 0.];errors=Float64[]
    # Three collinear overlapping bodies have an exact linear solution while
    # their complete graph remains active over this short physical interval.
    movement=(2/3-.3)*(-expm1(-30*.01))
    exact=copy(original);exact[1,1]-=movement;exact[1,3]+=movement
    for h in (.005,.0025,.00125)
        x=copy(original)
        for _ in 1:round(Int,.01/h)
            result=Interior.cellular_pair_relax(grid,ids,x,radii,CellBonds(),h;repulsion=10.)
            @test result.converged && result.relaxation_energy_loss>=0
            x=result.positions
        end
        push!(errors,sqrt(sum(abs2,x-exact)/length(x)))
    end
    @test errors[1]>errors[2]>errors[3]
    @test errors[2]/errors[3]>3.5
end
@testset "Radial cache preserves complete potential queries" begin
    cache=Interior._CellularRadialCache();ids=[1,2,3];radii=fill(.2,3);bonds=CellBonds()
    grid=Grid(8,8;extent=(4.,4.,1.));x=[.01 3.9 2.;1. 1. 2.;0. 0. 0.]
    for move in (0.,.005,.01,.1,.2)
        y=copy(x);y[1,1]+=move
        pairs=Interior._study_potential_pairs(cache,grid,y,radii)
        index=Dict(id=>i for (i,id) in enumerate(ids))
        cached=Interior._cellular_energy_forces(grid,ids,y,radii,bonds,index,10.,2.,Val(false),pairs).energy
        @test cached==Interior.cellular_mechanical_energy(grid,ids,y,radii,bonds)
    end
    @test cache.builds<5
    before=cache.builds
    Interior._study_potential_pairs(cache,grid,x,radii.*2)
    @test cache.builds==before+1
    before=cache.builds
    Interior._study_potential_pairs(cache,Grid(8,8,8;extent=(4.,4.,1.)),x,radii.*2)
    @test cache.builds==before+1
    funds=ones(3);y=copy(x)
    result=Interior.cellular_implicit_mechanics!(grid,ids,y,radii,CellBonds(),.01;
        relaxer=Interior._CellularRadialCache(),energy=funds,adhesion=zeros(3))
    @test abs(result.energy_residual)<1e-12
    @test result.relaxation_energy_loss>=0
end
