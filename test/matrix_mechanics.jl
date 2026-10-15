using Test, Interior, TOML

@testset "Matrix drag preserves the explicit unit path" begin
    g=Grid(16,16;extent=(8.,8.,1.));ids=[3,9,12]
    initial=[2. 2.7 5.;2. 2. 5.;0. 0. 0.];r=fill(.5,3)
    x=copy(initial);y=copy(initial);ex=ones(3);ey=copy(ex)
    motors=[.2 -.1 .15;.1 .2 -.1;0. 0. 0.]
    a=cellular_mechanics!(g,ids,x,r,CellBonds(),.2;motility=motors,energy=ex)
    b=cellular_mechanics!(g,ids,y,r,CellBonds(),.2;motility=motors,energy=ey,drag_multiplier=nothing)
    @test x==y && ex==ey && a==b
    @test_throws ArgumentError cellular_mechanics!(g,ids,copy(initial),r,CellBonds(),.1;drag_multiplier=[1.,.5,1.])
    @test_throws ArgumentError cellular_mechanics!(g,ids,copy(initial),r,CellBonds(),.1;drag_multiplier=[1.,NaN,1.])
end

@testset "Unequal drag changes passive relaxation, not its energy source" begin
    g=Grid(16,16;extent=(8.,8.,1.));ids=[1,2]
    x=[2. 2.6;2. 2.;0. 0.];initial=copy(x);drag=[1.,4.];r=fill(.5,2)
    residual=0.;loss=0.
    for _ in 1:100
        s=cellular_mechanics!(g,ids,x,r,CellBonds(),.001;adhesion=zeros(2),repulsion=2.,drag_multiplier=drag)
        residual+=s.energy_residual;loss+=s.relaxation_energy_loss
    end
    # Isolated overlapping pair: overlap decays at k*(μ1+μ2). Its
    # drag-weighted center is fixed; its ordinary center generally moves.
    expected=.4exp(-2*(1+1/4)*.1)
    @test abs((1-(x[1,2]-x[1,1]))-expected)<1e-4
    @test sum(drag.*x[1,:])≈sum(drag.*initial[1,:]) atol=1e-12
    @test abs(sum(x[1,:])-sum(initial[1,:]))>.01
    @test abs(residual)<1e-12 && loss>0
end

@testset "Motor funding pays the same held drag used in the outer reservation" begin
    g=Grid(16,16;extent=(8.,8.,1.));x=reshape([2.,2.,0.],3,1);energy=[.04]
    s=cellular_mechanics!(g,[1],x,[.5],CellBonds(),.2;motility=reshape([1.,0.,0.],3,1),
        energy,motor_cost=.5,drag_multiplier=[4.])
    @test x[1,1]-2≈sqrt(.04*.2/(.5*4)) atol=1e-12
    @test 0<=only(energy)<1e-12
    @test only(s.spent_motor_work)≈.04 atol=1e-12
    @test s.motor_drag_loss≈.04 atol=1e-12
    @test abs(s.energy_residual)<1e-12
end

@testset "Matrix observation survives death and is held through secretion" begin
    cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","arena","constructed-food.toml"))
    cfg["grid"]["shape"]=[8,8,1];cfg["cellular"]["population"]=3
    cfg["matrix_mechanics"]=Dict("resistance"=>20.)
    a=build_arena(cfg);p=only(filter(x->x isa ConstructiveEcology,a.processes));c=p.cells;h=p.habitat
    h.enzymes[1].=.01;h.enzymes[2].=.02
    context=Interior.CatalyticAnabolism(p)
    @test all(z->isapprox(z,1+20*.03/c.stencil.cell_volume;atol=1e-13),context.drag)
    observed=copy(context.drag);h.enzymes[1].*=2
    @test Interior._cellular_mechanical_drag(a.world,c,context)==observed
    Interior._cellular_select!(c,[3,1])
    @test Interior._cellular_mechanical_drag(a.world,c,context)==observed[[3,1]]
    @test_throws ArgumentError MatrixDrag(-1.)
    @test_throws ArgumentError MatrixDrag(Inf)
end

@testset "Matrix mechanics composes with paid construction, binding and restart" begin
    cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","arena","constructed-food.toml"))
    cfg["grid"]["shape"]=[8,8,1];cfg["cellular"]["population"]=6
    cfg["matrix_mechanics"]=Dict("resistance"=>20.)
    a=build_arena(cfg)
    for _ in 1:40;arena_step!(a);end
    p=only(filter(x->x isa ConstructiveEcology,a.processes));s=arena_snapshot(a)
    @test maximum(p.cells.mechanical.drag_multiplier)>1
    @test all(>=(0),p.cells.energy)
    @test abs(s["metrics"]["material_error"])<1e-4
    @test abs(s["metrics"]["chemical_energy_error"])<1e-4
    mktempdir() do dir
        file=joinpath(dir,"matrix.jls");checkpoint!(a,file);b=restore_arena(file)
        q=only(filter(x->x isa ConstructiveEcology,b.processes))
        @test q.habitat.cells===q.cells && last(b.processes).habitat===q.habitat
        @test q.mechanics.resistance==20.
        for _ in 1:10;arena_step!(a);arena_step!(b);end
        sa=arena_snapshot(a);sb=arena_snapshot(b);delete!(sa,"generated_at_unix");delete!(sb,"generated_at_unix")
        @test sa==sb
    end
end

@testset "Matrix drag identity through a mixed terminal and division step" begin
    cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","arena","constructed-food.toml"))
    cfg["grid"]["shape"]=[8,8,1];cfg["cellular"]["population"]=2
    cfg["matrix_mechanics"]=Dict("resistance"=>20.)
    original=build_arena(cfg);w=original.world
    p=only(filter(x->x isa ConstructiveEcology,original.processes));c=p.cells;h=p.habitat
    c.parameters=merge(c.parameters,(uptake_rate=0.,growth_rate=0.,catabolic_rate=0.,basal_rate=1.,
        neural_cost=0.,starvation_rate=10.,speed=0.))
    c.biomass.=[.26,2.];c.metabolites.=0;c.energy.=[0.,3.]
    c.weights.=0;c.weights[size(c.hidden,1)+10,end,:].=2
    h.parameters=merge(h.parameters,(secretion_rate=0.,))
    h.enzymes[1].=.03;h.enzymes[2].=.01
    expected=Interior.CatalyticAnabolism(p).drag[2]
    a=Arena(w,original.processes;dt=.2);arena_step!(a)
    @test c.ids==[2,3] && c.parents==[0,2] && h.ids==c.ids
    @test c.deaths==1 && c.births==1
    @test c.mechanical.ids==[2] && c.mechanical.drag_multiplier==[expected]
    @test abs(arena_snapshot(a)["metrics"]["material_error"])<1e-5
    @test abs(arena_snapshot(a)["metrics"]["chemical_energy_error"])<1e-5
end

@testset "Matrix drag does not overflow an avoidable catalyst sum" begin
    environment=zeros(4,1);environment[3:4,1].=floatmax(Float64)
    @test Interior._matrix_drag_coefficients(MatrixDrag(0.),environment)==[1.]
    @test only(Interior._matrix_drag_coefficients(MatrixDrag(1e-308),environment))≈1+2*(1e-308*floatmax(Float64))
    @test_throws ArgumentError MatrixDrag(big"1e400")
    @test Interior._matrix_scaled_amount(floatmax(Float64),0.,.5)==0
    @test Interior._matrix_scaled_amount(1e308,1e-308,.5)≈2.
    @test Interior._matrix_scaled_amount(1e-308,1e308,1e308)≈1e-308 rtol=1e-14
end
