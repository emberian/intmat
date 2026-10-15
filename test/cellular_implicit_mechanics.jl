# This experiment is intentionally not part of the default Interior include tree.
isdefined(Interior,:cellular_proximal_relax) || Base.include(Interior,joinpath(@__DIR__,"..","studies","cellular_implicit_mechanics.jl"))
@testset "Opt-in proximal contact mechanics" begin
    relax=Interior.cellular_proximal_relax
    grid=Grid(16,16;extent=(4.,4.,1.));ids=[1,2];radii=fill(.5,2)
    x=[1. 1.6;1. 1.;0. 0.];bonds=CellBonds()
    # Analytic two-body overlap: o_new=o_old/(1+2kh).
    for h in (.002,.02,.2)
        result=relax(grid,ids,x,radii,bonds,h;repulsion=10.,tolerance=1e-10)
        @test result.converged
        @test 1-(result.positions[1,2]-result.positions[1,1])≈.4/(1+20h) atol=1e-10
        @test result.potential_after<=result.potential_before
        @test result.proximal_objective<=result.potential_before
        @test result.relaxation_energy_loss>=0
    end
    # Passive tensile spring, same analytic backward-Euler contraction.
    spring=CellBonds();add_cell_bond!(spring,1,2,.4)
    result=relax(grid,ids,x,fill(.05,2),spring,.2;stiffness=3.,tolerance=1e-10)
    @test result.converged
    @test result.positions[1,2]-result.positions[1,1]≈.4+.2/(1+1.2) atol=1e-10
    # Force-free translation and isolated bodies remain unchanged.
    separated=[.2 2.;.1 2.;0. 0.]
    result=relax(grid,ids,separated,fill(.1,2),bonds,.5)
    @test result.converged && result.positions==separated
    # Periodic seam uses the local proximal chart, rather than a long wrap move.
    seam=[.1 3.9;1. 1.;0. 0.]
    result=relax(grid,ids,seam,fill(.2,2),bonds,.1;tolerance=1e-10)
    @test result.converged
    @test Interior._cellular_length(Interior._cellular_delta(grid,result.positions,1,2))≈.4-.2/3 atol=1e-10
    # No mutation on failure: caller can retry a smaller macrostep safely.
    original=copy(x);edges=copy(spring.endpoints);rests=copy(spring.rest_lengths)
    result=relax(grid,ids,x,radii,spring,1.;max_iterations=1,tolerance=1e-14)
    @test x==original && spring.endpoints==edges && spring.rest_lengths==rests
    # Timestep refinement converges first-order to continuous two-body relaxation.
    errors=Float64[]
    for h in (.05,.025,.0125)
        current=copy(x)
        for _ in 1:round(Int,.2/h)
            step=relax(grid,ids,current,radii,bonds,h;tolerance=1e-10)
            @test step.converged
            current=step.positions
        end
        push!(errors,abs(1-(current[1,2]-current[1,1])-.4exp(-4)))
    end
    @test errors[1]>errors[2]>errors[3]
    @test errors[2]/errors[3]>1.7
end
@testset "Funded implicit mechanics and topology ledgers" begin
    implicit! = Interior.cellular_implicit_mechanics!
    grid=Grid(12,12;extent=(4.,4.,1.));ids=[1,2];radii=fill(.2,2)
    for funds in (nothing,[0.,0.],[.0001,.0001],[1.,1.])
        x=[1. 1.3;1. 1.;0. 0.];bonds=CellBonds();before=deepcopy(funds)
        motor=[.5 -.5;0. 0.;0. 0.]
        result=implicit!(grid,ids,x,radii,bonds,.1;energy=funds,motility=motor,tolerance=1e-9)
        @test result.relaxation_energy_loss>=0
        @test abs(result.energy_residual)<1e-12
        @test result.max_residual<1e-7
        @test result.motor_drag_loss>=0
        if funds!==nothing
            @test all(>=(0),funds)
            @test before-funds≈result.spent_motor_work atol=1e-15
            @test all(result.spent_motor_work.<=before)
        end
    end
    # Invalid controls and impossible iteration budget do not publish partial state.
    x=[1. 1.3;1. 1.;0. 0.];bonds=CellBonds();funds=ones(2);original=copy(x)
    @test_throws ArgumentError implicit!(grid,ids,x,radii,bonds,.1;energy=funds,max_iterations=0)
    @test x==original && funds==ones(2) && isempty(bonds.endpoints)
    # Consent formation and withdrawal use actual topology energy loss.
    result=implicit!(grid,ids,x,radii,bonds,.1;adhesion=ones(2),stiffness=3.)
    @test result.formed==1 && length(bonds.endpoints)==1
    x[1,2]=x[1,1]+.6
    result=implicit!(grid,ids,x,radii,bonds,0.;adhesion=zeros(2),stiffness=3.)
    @test result.broken==1 && result.topology_energy_loss>0
    @test abs(result.energy_residual)<1e-12
end
@testset "Implicit Hessian and nonlinear solver controls" begin
    grid=Grid(8,8,8;extent=(4.,4.,4.));ids=[1,2,3]
    x=[1. 1.4 1.9;1. 1.2 1.1;1. 1.1 1.3];radii=fill(.35,3)
    bonds=CellBonds();add_cell_bond!(bonds,1,3,.5)
    index=Dict(id=>i for (i,id) in enumerate(ids));h=.05;displacement=zeros(3,3)
    current=Interior._implicit_state(grid,ids,x,displacement,radii,bonds,index,10.,3.,h)
    edges,diagonal=Interior._implicit_edges(grid,ids,x,radii,bonds,index,10.,3.,current.pairs,h,true)
    v=[.2 -.1 .3;-.1 .4 .2;.3 -.2 -.1];hv=similar(v)
    Interior._implicit_hessian!(hv,v,edges,h)
    epsilon=1e-6
    plus=Interior._implicit_state(grid,ids,x,epsilon*v,radii,bonds,index,10.,3.,h)
    minus=Interior._implicit_state(grid,ids,x,-epsilon*v,radii,bonds,index,10.,3.,h)
    @test hv≈(plus.gradient-minus.gradient)/(2epsilon) rtol=1e-7 atol=1e-7
    @test all(>(0),diagonal)
    # Negative curvature is truncated to a descent direction, never a zero step.
    edge=Interior._ImplicitContactEdge(1,2,(1.,0.,0.),1.,-50.)
    gradient=[0. 0.;1. -1.;0. 0.]
    direction,_=Interior._implicit_cg(gradient,ones(3,2),[edge],1.,20)
    @test Interior._implicit_dot(gradient,direction)<0
    for exact_curvature in (false,true)
        result=Interior.cellular_proximal_relax(grid,ids,x,radii,bonds,.02;exact_curvature)
        @test result.converged
        @test result.potential_after<=result.potential_before
    end
end
@testset "Implicit discovery and crowded-state convergence" begin
    grid=Grid(12,12;extent=(6.,6.,1.));ids=[1,2,3];radii=fill(.5,3);bonds=CellBonds()
    x=[1. 1.1 2.15;1. 1. 1.;0. 0. 0.]
    @test !((2,3) in cellular_pairs(grid,x,radii))
    result=Interior.cellular_proximal_relax(grid,ids,x,radii,bonds,.02;tolerance=1e-9)
    @test result.converged
    @test (2,3) in cellular_pairs(grid,result.positions,radii)
    @test result.potential_after≈Interior.cellular_mechanical_energy(grid,ids,result.positions,radii,bonds) atol=1e-14
    # A real nonlinear multi-body state, rather than uncoupled radial springs.
    grid3=Grid(8,8,8;extent=(3.,3.,3.));ids3=collect(1:4);r3=fill(.45,4)
    x3=[1. 1.3 1.2 1.1;1. 1.1 1.4 1.2;1. 1.2 1.1 1.4]
    before=copy(x3);funds=ones(4)
    @test_throws ArgumentError Interior.cellular_implicit_mechanics!(grid3,ids3,x3,r3,bonds,.1;
        energy=funds,max_iterations=1,max_substeps=1,tolerance=1e-12)
    @test x3==before && funds==ones(4) && isempty(bonds.endpoints)
    reference=copy(x3)
    for _ in 1:400
        cellular_mechanics!(grid3,ids3,reference,r3,CellBonds(),.00025;adhesion=zeros(4))
    end
    errors=Float64[]
    for h in (.0125,.00625,.003125)
        y=copy(x3)
        for _ in 1:round(Int,.1/h)
            step=Interior.cellular_proximal_relax(grid3,ids3,y,r3,bonds,h;tolerance=1e-9)
            @test step.converged && step.relaxation_energy_loss>=0
            y=step.positions
        end
        push!(errors,sqrt(sum(abs2,y-reference)/length(y)))
    end
    @test errors[1]>errors[2]>errors[3]
    @test errors[2]/errors[3]>1.6
end
@testset "Implicit objective budget is transactional" begin
    grid=Grid(8,8;extent=(4.,4.,1.));ids=[1,2];x=[1. 1.3;1. 1.;0. 0.];before=copy(x)
    funds=ones(2);bonds=CellBonds()
    @test_throws ArgumentError Interior.cellular_implicit_mechanics!(grid,ids,x,[.2,.2],bonds,.1;energy=funds,max_evaluations=1)
    @test x==before&&funds==ones(2)&&isempty(bonds.endpoints)
end
