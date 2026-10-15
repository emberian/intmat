isdefined(Interior,:cellular_implicit_mechanics!)||Base.include(Interior,joinpath(@__DIR__,"..","studies","cellular_implicit_mechanics.jl"))
isdefined(Interior,:cellular_pair_relax)||Base.include(Interior,joinpath(@__DIR__,"..","studies","cellular_pair_proximal.jl"))
isdefined(Interior,:cellular_coordinate_relax)||Base.include(Interior,joinpath(@__DIR__,"..","studies","cellular_coordinate_proximal.jl"))
@testset "Local coordinate proximal mechanics" begin
    grid=Grid(8,8;extent=(4.,4.,1.));ids=[1,2];radii=fill(.5,2);bonds=CellBonds()
    x=[1. 1.6;1. 1.;0. 0.]
    for h in (.001,.01,.02)
        result=Interior.cellular_coordinate_relax(grid,ids,x,radii,bonds,h;tolerance=1e-9)
        @test result.converged
        @test result.potential_after<=result.potential_before
        @test result.residual<1e-7
    end
    # A wide local move leaves the conservative chart and must request a
    # smaller macrostep, rather than accept an unconverged clipped position.
    rejected=Interior.cellular_coordinate_relax(grid,ids,x,radii,bonds,.05;tolerance=1e-9)
    @test !rejected.converged && x==[1. 1.6;1. 1.;0. 0.]
    # Opposing tensile and repulsive forces cancel at d=.9 exactly as a model
    # equilibrium, unlike an isolated sequence of nonzero edge flows.
    spring=CellBonds();add_cell_bond!(spring,1,2,.4)
    equilibrium=[1. 1.9;1. 1.;0. 0.]
    result=Interior.cellular_coordinate_relax(grid,ids,equilibrium,radii,spring,.2)
    @test result.converged && result.positions==equilibrium
    errors=Float64[]
    for h in (.01,.005,.0025)
        y=copy(x)
        for _ in 1:round(Int,.02/h)
            result=Interior.cellular_coordinate_relax(grid,ids,y,radii,bonds,h;tolerance=1e-9)
            @test result.converged
            y=result.positions
        end
        push!(errors,abs(1-(y[1,2]-y[1,1])-.4exp(-.4)))
    end
    @test errors[1]>errors[2]>errors[3]
    @test errors[2]/errors[3]>1.7
    funds=fill(.01,2);before=copy(funds);y=copy(x)
    result=Interior.cellular_implicit_mechanics!(grid,ids,y,radii,bonds,.01;
        energy=funds,relaxer=Interior._CellularCoordinateCache(),motility=[.3 -.3;0. 0.;0. 0.])
    @test funds≈before-result.spent_motor_work atol=1e-15
    @test result.relaxation_energy_loss>=0 && abs(result.energy_residual)<1e-12
end
@testset "Symmetric midpoint mechanics and motor quadrature" begin
    grid=Grid(8,8;extent=(4.,4.,1.));ids=[1,2];radii=fill(.5,2);original=[1. 1.6;1. 1.;0. 0.]
    passive_errors=Float64[];active_errors=Float64[];work_errors=Float64[]
    for h in (.01,.005,.0025)
        x=copy(original);cache=Interior._CellularCoordinateCache(midpoint=true)
        for _ in 1:round(Int,.02/h)
            r=Interior.cellular_coordinate_relax(grid,ids,x,radii,CellBonds(),h;cache,tolerance=1e-10)
            @test r.converged && r.relaxation_energy_loss>=0
            x=r.positions
        end
        push!(passive_errors,abs(1-(x[1,2]-x[1,1])-.4exp(-.4)))
        x=copy(original);funds=fill(10.,2);work=0.;bonds=CellBonds()
        for _ in 1:round(Int,.02/h)
            r=Interior.cellular_implicit_mechanics!(grid,ids,x,radii,bonds,h;energy=funds,adhesion=zeros(2),
                motility=[.2 -.2;0. 0.;0. 0.],symmetric_motor=true,
                relaxer=Interior._CellularCoordinateCache(midpoint=true),tolerance=1e-10)
            @test r.relaxation_energy_loss>=0 && abs(r.energy_residual)<1e-12
            work+=sum(r.spent_motor_work)
        end
        overlap=.02+.38exp(-.4)
        exact_potential_work=2*10*.2*(.02*.02+.38*(-expm1(-.4))/20)
        exact_drag=2*.2^2*.02
        push!(active_errors,abs(1-(x[1,2]-x[1,1])-overlap))
        push!(work_errors,abs(work-exact_potential_work-exact_drag))
        @test all(>=(0),funds)
        @test sum(10 .-funds)≈work atol=1e-13
    end
    for errors in (passive_errors,active_errors,work_errors)
        @test errors[1]>errors[2]>errors[3]
        @test errors[2]/errors[3]>3.5
    end
end
@testset "Coordinate endpoint gate shares exact energy arithmetic" begin
    grid=Grid(8,8,8;extent=(4.,4.,4.));ids=[1,2,3]
    x=[.01 3.9 .2;1. 1.1 1.;1. 1. 1.2]
    edges=[Interior._CellularLocalEdge(2,.5,10.,false),Interior._CellularLocalEdge(3,.1,3.,true)]
    for shift in (-.1,0.,.1)
        z=x[:,1].+[shift,0.,0.]
        state=Interior._coordinate_state(grid,ids,1,z,x,x[:,1],edges,.1)
        @test Interior._coordinate_energy(grid,z,x,edges)==state.energy
    end
end
@testset "Symmetric motor funding remains bounded at exhaustion" begin
    grid=Grid(8,8;extent=(4.,4.,1.));ids=[1,2];radii=fill(.5,2)
    for amount in (0.,1e-9,.001,.1)
        x=[1. 1.6;1. 1.;0. 0.];funds=fill(amount,2);before=copy(funds)
        result=Interior.cellular_implicit_mechanics!(grid,ids,x,radii,CellBonds(),.02;
            energy=funds,adhesion=zeros(2),motility=[1. -1.;0. 0.;0. 0.],symmetric_motor=true,
            relaxer=Interior._CellularCoordinateCache(midpoint=true))
        @test all(>=(0),funds)
        @test before-funds≈result.spent_motor_work atol=1e-15
        @test all(result.spent_motor_work.<=before)
        @test result.relaxation_energy_loss>=0 && abs(result.energy_residual)<1e-12
        amount==0&&@test all(==(0),result.active_distance)
    end
end
@testset "Local Hessian and passive center refinement" begin
    grid=Grid(8,8,8;extent=(4.,4.,4.));ids=[1,2,3]
    x=[.01 3.9 .2;1. 1.1 1.;1. 1. 1.2];z=x[:,1];v=[.2,-.1,.3]
    edges=[Interior._CellularLocalEdge(2,.5,10.,false),Interior._CellularLocalEdge(3,.1,3.,true)]
    state=Interior._coordinate_state(grid,ids,1,z,x,z,edges,.1);epsilon=1e-6
    plus=Interior._coordinate_state(grid,ids,1,z+epsilon*v,x,z,edges,.1)
    minus=Interior._coordinate_state(grid,ids,1,z-epsilon*v,x,z,edges,.1)
    @test state.hessian*v≈(plus.gradient-minus.gradient)/(2epsilon) rtol=1e-7 atol=1e-7
    # Coordinate time splitting has a center error even with equal/opposite
    # pair forces. Expose its second-order convergence; do not claim conservation.
    grid2=Grid(8,8;extent=(4.,4.,1.));original=[1. 1.6;1. 1.;0. 0.];drifts=Float64[]
    for h in (.01,.005,.0025)
        y=copy(original);cache=Interior._CellularCoordinateCache(midpoint=true)
        for _ in 1:round(Int,.02/h)
            result=Interior.cellular_coordinate_relax(grid2,[1,2],y,[.5,.5],CellBonds(),h;cache,tolerance=1e-10)
            @test result.converged
            y=result.positions
        end
        push!(drifts,abs(sum(y[1,:]-original[1,:])/2))
    end
    @test drifts[1]>drifts[2]>drifts[3]>0
    @test drifts[2]/drifts[3]>3.5
end
@testset "Fixed outer clock absorbs only floating-point tails" begin
    grid=Grid(8,8;extent=(4.,4.,1.));ids=[1,2]
    for splits in (4,8)
        x=[1. 1.6;1. 1.;0. 0.]
        result=Interior.cellular_implicit_mechanics!(grid,ids,x,[.5,.5],CellBonds(),.1;
            max_step=.1/splits,adhesion=zeros(2),symmetric_motor=true,
            relaxer=Interior._CellularCoordinateCache(midpoint=true),max_evaluations=32)
        @test result.substeps==splits && result.retries==0
        @test result.relaxation_energy_loss>=0 && abs(result.energy_residual)<1e-12
    end
end
