using Test, Interior
Base.include(Interior,joinpath(@__DIR__,"..","studies","cellular_implicit_mechanics.jl"))
Base.include(Interior,joinpath(@__DIR__,"..","studies","cellular_pair_proximal.jl"))
Base.include(Interior,joinpath(@__DIR__,"..","studies","cellular_rkc_mechanics.jl"))
@testset "RKC2 stability and order" begin
 for stages in (2,4,8,16)
  coeff=Interior.cellular_rkc_coefficients(stages)
  @test abs(only(Interior.cellular_rkc_step(x->-x,[1.],coeff.stability_interval*.99,coeff)))<=1
  errors=Float64[]
  for steps in (8,16,32)
   x=[1.];h=.4/steps
   for _ in 1:steps;x=Interior.cellular_rkc_step(y->-3y,x,h,coeff);end
   push!(errors,abs(x[1]-exp(-1.2)))
  end
  @test errors[1]/errors[2]>3.7
  @test errors[2]/errors[3]>3.7
 end
end
@testset "passive RKC soft-contact flow" begin
 grid=Grid(16,16,1;extent=(16,16,1));ids=UInt64[1,2];radii=[.5,.5]
 x=[4. 4.8;4. 4.;0. 0.];b=CellBonds()
 error=Float64[]
 for parts in (1,2,4)
  y=copy(x);cache=Interior._CellularRKCCache()
  for _ in 1:parts
   result=Interior.cellular_rkc_relax(grid,ids,y,radii,b,.05/parts;repulsion=8.,cache)
   @test result.converged
   @test result.potential_after<=result.potential_before
   y=result.positions
  end
  exact_distance=1-.2exp(-16*.05)
  push!(error,abs(y[1,2]-y[1,1]-exact_distance))
  @test sum(y[1,:])≈sum(x[1,:]) atol=1e-14
 end
 @test error[1]/error[2]>3.5
 @test error[2]/error[3]>3.5
end
@testset "funded stabilized flow and work convergence" begin
 grid=Grid(8,8;extent=(4.,4.,1.));ids=[1,2];radii=fill(.5,2);original=[1. 1.6;1. 1.;0. 0.]
 position_errors=Float64[];work_errors=Float64[]
 for h in (.01,.005,.0025)
  x=copy(original);funds=fill(10.,2);work=0.;cache=Interior._CellularRKCCache()
  for _ in 1:round(Int,.02/h)
   result=Interior.cellular_implicit_mechanics!(grid,ids,x,radii,CellBonds(),h;energy=funds,
    adhesion=zeros(2),motility=[.2 -.2;0. 0.;0. 0.],symmetric_motor=true,relaxer=cache)
   @test result.relaxation_energy_loss>=0&&abs(result.energy_residual)<1e-12
   work+=sum(result.spent_motor_work)
  end
  overlap=.02+.38exp(-.4)
  expected_work=2*10*.2*(.02*.02+.38*(-expm1(-.4))/20)+2*.2^2*.02
  push!(position_errors,abs(1-(x[1,2]-x[1,1])-overlap))
  push!(work_errors,abs(work-expected_work))
  @test all(>=(0),funds)
  @test sum(10 .-funds)≈work atol=1e-13
 end
 for errors in (position_errors,work_errors)
  @test errors[1]>errors[2]>errors[3]
  @test errors[2]/errors[3]>3.5
 end
 for amount in (0.,1e-9,.001,.1)
  x=copy(original);funds=fill(amount,2)
  result=Interior.cellular_implicit_mechanics!(grid,ids,x,radii,CellBonds(),.02;energy=funds,
   adhesion=zeros(2),motility=[1. -1.;0. 0.;0. 0.],symmetric_motor=true,relaxer=Interior._CellularRKCCache())
  @test all(>=(0),funds)
  @test fill(amount,2)-funds≈result.spent_motor_work atol=1e-15
 end
end

Base.include(Interior,joinpath(@__DIR__,"..","studies","cellular_rkc_funding.jl"))
@testset "cached funding schedule matches original RKC" begin
 grid=Grid(16,16;extent=(8.,8.,1.));ids=[5,2,9];radii=fill(.5,3)
 original=[1. 1.7 2.3;1. 1.1 1.;0. 0. 0.]
 for amount in (0.,.001,1.),splits in (1,4),adhesion in (zeros(3),ones(3))
  x=copy(original);y=copy(original);a=CellBonds();b=CellBonds();ea=fill(amount,3);eb=copy(ea)
  kwargs=(motility=[.2 -.1 -.1;0. .02 -.02;0. 0. 0.],adhesion=adhesion,max_step=.02/splits,max_evaluations=256)
  left=Interior.cellular_implicit_mechanics!(grid,ids,x,radii,a,.02;energy=ea,relaxer=Interior._CellularRKCCache(),symmetric_motor=true,kwargs...)
  right=Interior.cellular_rkc_funded_mechanics!(grid,ids,y,radii,b,.02;energy=eb,kwargs...)
  @test isequal(x,y)
  @test isequal(ea,eb)
  @test a.endpoints==b.endpoints&&a.rest_lengths==b.rest_lengths
  @test isequal(left.spent_motor_work,right.spent_motor_work)
 end
end
