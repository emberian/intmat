using Test, Interior
isdefined(Interior,:DielectricGrowth) || Base.include(Interior,joinpath(@__DIR__,"..","src","electrical.jl"))
@testset "Dielectric growth" begin
 w=World(Grid(15,15,15);seed=29);initial=Float64(material_mass(w))
 p=Interior.DielectricGrowth(w;eta=3.,tolerance=1e-5)
 @test p.last_residual<=p.tolerance
 @test all(==(0),p.potential[1,:,:])
 @test p.potential[p.occupied]==[1f0]
 @test minimum(p.potential)>=-1e-5 && maximum(p.potential)<=1+1e-5
 initial_potential=copy(p.potential)
 for _ in 1:15;step!(w,p,.05f0);end
 @test p.events==15
 @test p.last_residual<=p.tolerance
 @test count(p.occupied)==16
 @test p.potential != initial_potential
 @test sum(w.fields.deposit)≈16p.deposit_per_cell atol=2e-5
 @test sum(Float64,w.fields.nutrient)+sum(Float64,w.fields.deposit)≈initial atol=2e-5
 # Every grown site remains face-connected to the conductor seed.
 seen=falses(size(p.occupied));stack=[CartesianIndex(8,8,8)];seen[stack[1]]=true
 while !isempty(stack)
  I=pop!(stack)
  for J in Interior.electrical_neighbors(I)
   if checkbounds(Bool,seen,J) && p.occupied[J] && !seen[J]
    seen[J]=true;push!(stack,J)
   end
  end
 end
 @test seen==p.occupied
 w2=World(Grid(15,15,15);seed=29);p2=Interior.DielectricGrowth(w2;eta=3.,tolerance=1e-5)
 for _ in 1:75;step!(w2,p2,.01f0);end
 @test p2.events==p.events
 @test p2.occupied==p.occupied
 @test w2.fields.deposit==w.fields.deposit
 @test p2.potential==p.potential
 @test_throws ArgumentError Interior.DielectricGrowth(World(Grid(4,4,4)))
 @test_throws ArgumentError Interior.DielectricGrowth(World(Grid(8,8,8));event_interval=0)
 @test_throws ArgumentError Interior.DielectricGrowth(World(Grid(8,8,8));root=(1,3,3))
end
@testset "Electrical harmonic benchmark and controls" begin
 w=World(Grid(5,5,5));p=Interior.DielectricGrowth(w;tolerance=2e-6)
 # Independently solve the free-cell linear system around the held central seed.
 free=[I for I in CartesianIndices(p.potential) if !Interior.electrical_boundary(I,size(p.potential)) && !p.occupied[I]]
 lookup=Dict(I=>i for (i,I) in enumerate(free));A=zeros(Float64,length(free),length(free));b=zeros(Float64,length(free))
 for (i,I) in enumerate(free)
  A[i,i]=6
  for J in Interior.electrical_neighbors(I)
   if p.occupied[J];b[i]+=1
   elseif haskey(lookup,J);A[i,lookup[J]]-=1
   end
  end
 end
 exact=A\b
 @test maximum(abs.([p.potential[I] for I in free].-exact))<2e-5
 olddeposit=copy(w.fields.deposit);oldclock=w.time;p.event_interval=NaN
 @test_throws ArgumentError step!(w,(WaveForcing(),p),.02f0)
 @test w.fields.deposit==olddeposit && w.time==oldclock
 @test all(iszero,w.fields.pressure)
 p.event_interval=.05;p.frontier=falses(6,6,6)
 @test_throws DimensionMismatch step!(w,p,.02f0)
 # Failed initialization leaves the world untouched.
 badw=World(Grid(15,15,15));before=copy(badw.fields.nutrient)
 @test_throws ErrorException Interior.DielectricGrowth(badw;max_iterations=1,tolerance=1e-12)
 @test badw.fields.nutrient==before && all(iszero,badw.fields.deposit)
end
