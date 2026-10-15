using Test, Random
using Interior: Grid, CellBonds, add_cell_bond!, cellular_pairs, cellular_mechanics!, exchange_cell_metabolites!

@testset "Periodic displacement fast path matches general modulo" begin
    rng=MersenneTwister(413)
    for L in (Float64(Float32(.1)),1.,24.,32.,1e-30,1e30)
        points=vcat((rand(rng,1000).*8 .-4).*L,
            [-0.,0.,prevfloat(L),L,nextfloat(L),-L,prevfloat(-L),nextfloat(-L),2L,-2L])
        @test all(isequal(Interior._cellular_mod_period(x,L),mod(x,L)) for x in points)
    end
end

@testset "Physical periodic cell contacts" begin
    g=Grid(20,10;extent=(10.,4.,1.))
    x=[.05 9.95 5.; 2. 2. 2.; 0. 0. 0.];r=fill(.1,3)
    @test cellular_pairs(g,x,r)==[(1,2)]
    @test cellular_pairs(g,x,r;range=5.)==[(1,2),(1,3),(2,3)]
    # Compare sparse physical bins to a direct oracle, including unequal radii
    # and grids whose periodic neighborhoods alias themselves.
    rng=MersenneTwister(591)
    for extent in ((1.,1.,1.),(10.,2.,3.)), nz in (1,5), scale in (1.,.001)
        grid=Grid(7,8,nz;extent)
        points=rand(rng,3,45).*collect(extent); radii=scale.*(.01 .+ .2rand(rng,45))
        gap=.07scale
        expected=[(i,j) for i in 1:44 for j in i+1:45 if
            sqrt(sum((mod(points[d,j]-points[d,i]+Float64(grid.extent[d])/2,Float64(grid.extent[d]))-Float64(grid.extent[d])/2)^2 for d in 1:(nz==1 ? 2 : 3)))<=radii[i]+radii[j]+gap]
        @test cellular_pairs(grid,points,radii;range=gap)==expected
    end
    ids=[11,42,6];b=CellBonds()
    initial=copy(x)
    result=cellular_mechanics!(g,ids,x,r,b,.1;adhesion=zeros(3))
    @test isempty(b.endpoints)
    @test result.potential_after<result.potential_before
    @test result.relaxation_energy_loss>0
    @test abs(result.energy_residual)<1e-12
    # Seam-local displacement, not arithmetic means of wrapped positions.
    displacement=[mod(x[d,i]-initial[d,i]+Float64(g.extent[d])/2,Float64(g.extent[d]))-Float64(g.extent[d])/2 for d in 1:2,i in 1:3]
    @test maximum(abs,sum(displacement;dims=2))<1e-12
    @test x[:,3]==initial[:,3]
    @test all(x[d,i]>=0 && x[d,i]<g.extent[d] for d in 1:3,i in 1:3)
end

@testset "Persistent cellular identity and consent" begin
    g=Grid(20,20;extent=(10.,10.,1.)); ids=[81,17,56]
    x=[2. 2.35 2.7; 2. 2. 2.; 0. 0. 0.]; r=fill(.2,3)
    b=CellBonds();add_cell_bond!(b,81,17,.4)
    @test !add_cell_bond!(b,17,81,.4)
    @test b.endpoints==[(17,81)]
    perm=[3,1,2]; y=x[:,perm]; c=deepcopy(b)
    d1=cellular_mechanics!(g,ids,x,r,b,.15;max_degree=1)
    d2=cellular_mechanics!(g,ids[perm],y,r[perm],c,.15;max_degree=1)
    @test x≈y[:,invperm(perm)] atol=1e-12
    @test b.endpoints==c.endpoints
    @test d1.potential_after≈d2.potential_after
    @test length(b.endpoints)==1
    # Lost endpoint IDs are removed regardless of the new array's length.
    cellular_mechanics!(g,[81,56],x[:,[1,3]],r[[1,3]],b,0.)
    @test isempty(b.endpoints)
    add_cell_bond!(b,81,56,.1)
    result=cellular_mechanics!(g,[81,56],x[:,[1,3]],r[[1,3]],b,0.;adhesion=[0.,1.])
    @test isempty(b.endpoints)
    @test result.broken==1
    @test result.topology_energy_loss>0
    @test abs(result.energy_residual)<1e-12
    # A tensile seam bond pulls equal and opposite across periodic coordinates.
    x=[.1 9.6; 2. 2.; 0. 0.];b=CellBonds();add_cell_bond!(b,4,9,.3)
    result=cellular_mechanics!(g,[4,9],x,[.1,.1],b,.1;break_strain=2.)
    @test x[1,1]<.1 && x[1,2]>9.6
    @test x[1,1]+x[1,2]≈9.7
    @test result.potential_after<result.potential_before
    # Coincident points receive a reproducible symmetric separating direction.
    x=[2. 2.;2. 2.;0. 0.]
    cellular_mechanics!(g,[1,2],x,[.1,.1],CellBonds(),.1;adhesion=[0.,0.])
    @test x[1,1]<2<x[1,2]
end

@testset "Funded cellular motors" begin
    g=Grid(10,10;extent=(5.,5.,1.))
    x=[2. 2.4;2. 2.;0. 0.];b=CellBonds();add_cell_bond!(b,1,2,.4)
    energy=[.001,.001];before=copy(energy)
    result=cellular_mechanics!(g,[1,2],x,[.1,.1],b,.3;
        motility=[-10. 10.;0. 0.;0. 0.],energy,break_strain=100.)
    @test all(energy.>=0)
    @test before-energy≈result.spent_motor_work atol=1e-14
    @test sum(result.spent_motor_work)>=result.motor_drag_loss
    @test sum(result.spent_motor_work)>=result.motor_potential_work
    @test abs(result.energy_residual)<1e-12
    # Zero budget cannot produce motor motion, even for an isolated cell.
    x=reshape([1.,1.,0.],3,1)
    result=cellular_mechanics!(g,[8],x,[.1],CellBonds(),1.;motility=ones(3,1),energy=[0.])
    @test x==reshape([1.,1.,0.],3,1)
    @test result.spent_motor_work==[0.]
    # Paid translation can cross the seam; no clipping at the domain boundary.
    x=reshape([4.9,1.,0.],3,1);energy=[1.]
    result=cellular_mechanics!(g,[8],x,[.1],CellBonds(),.5;
        motility=reshape([1.,0.,0.],3,1),energy)
    @test x[1]≈.4 atol=1e-12
    @test result.active_distance[1]≈.5
    @test energy[1]≈.5
end

@testset "Bond-mediated conservative intracellular exchange" begin
    b=CellBonds();add_cell_bond!(b,10,30,1.);add_cell_bond!(b,30,50,1.)
    ids=[10,30,50];amounts=[10. 0. 0.;0. 2. 8.];total=sum(amounts;dims=2)
    p=[3,1,2];other=amounts[:,p]
    exchange_cell_metabolites!(amounts,b,ids,[1.,2.,3.],10.;conductance=3.)
    exchange_cell_metabolites!(other,b,ids[p],[1.,2.,3.][p],10.;conductance=3.)
    @test sum(amounts;dims=2)≈total atol=1e-14
    @test all(amounts.>=0)
    @test amounts≈other[:,invperm(p)] atol=1e-14
    before=copy(amounts)
    exchange_cell_metabolites!(amounts,b,ids,ones(3),10.;conductance=3.,openness=zeros(3))
    @test amounts==before
    pair=CellBonds();add_cell_bond!(pair,1,2,1.)
    amounts=[3. 0.]
    exchange_cell_metabolites!(amounts,pair,[1,2],[1.,2.],1.;conductance=1e6)
    @test amounts≈[1. 2.]
end

@testset "Cached broadphase admits newly closing contacts" begin
    grid=Grid(20,20;extent=(4.,4.,1.))
    x=[1. 1.7 2.4; 2. 2. 2.;0. 0. 0.]
    radii=[.1,.18,.12];energy=fill(20.,3)
    result=cellular_mechanics!(grid,[2,4,8],x,radii,CellBonds(),.15;
        motility=[4. 0. -4.;0. 0. 0.;0. 0. 0.],energy,adhesion=zeros(3))
    @test result.substeps>1
    @test result.potential_after>0
    @test abs(result.energy_residual)<1e-11
    @test all(energy.>=0)
end
@testset "Cached contact skin preserves bond formation" begin
    rng=MersenneTwister(384)
    grid=Grid(20,20;extent=(4.,4.,1.))
    ids=randperm(rng,45); radii=.08 .+ .08rand(rng,45)
    old=rand(rng,3,45).*reshape([4.,4.,0.],3,1)
    skin=.4; gap=.07
    candidates=cellular_pairs(grid,old,radii;range=skin+gap)
    moved=copy(old)
    moved[1:2,:].=mod.(moved[1:2,:] .+ .1 .* (rand(rng,2,45) .- .5),4.)
    direct=CellBonds();cached=CellBonds();index=Dict(id=>i for (i,id) in enumerate(ids))
    args=(grid,ids,moved,radii)
    a=Interior._cellular_refresh_bonds!(args...,direct,index,ones(45),.5,.75,3,gap,2.)
    b=Interior._cellular_refresh_bonds!(args...,cached,index,ones(45),.5,.75,3,gap,2.;pairs=candidates)
    @test a==b
    @test direct.endpoints==cached.endpoints
    @test direct.rest_lengths==cached.rest_lengths
end

@testset "Low-precision junction exchange stays nonnegative" begin
    b=CellBonds(); add_cell_bond!(b,1,2,1.)
    amounts=Float32[1 9e-8]; before=sum(Float64,amounts)
    exchange_cell_metabolites!(amounts,[1,2],b,100.;rate=1.,volumes=[1e8,1.])
    @test all(>=(0),amounts)
    @test sum(Float64,amounts)≈before atol=1e-14
end
