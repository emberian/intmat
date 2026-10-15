using Interior,Test,Random
Base.include(Interior,joinpath(@__DIR__,"..","studies","device_contact_candidates.jl"))
function candidate_fixture(;n=16,capacity=n+4,dense=false,density=128.,bin_width=1.)
 w=World(Grid(8,8,8;extent=(4,4,4)));fill!(w.fields.nutrient,0)
 r=ResourcePool(w;potentials=(4,4,0),initial=(.2,.3,.4),diffusion=(0,0,0))
 positions=rand(MersenneTwister(41),3,n).*4
 dense&&(positions.=1 .+ 0.04 .* positions)
 c=CellularEcology(w,positions,r;capacity,density,speed=.03)
 p=Interior.DeviceCellular(w,c;mechanical_substeps=16,bin_width,energy_type=Float64)
 p.state.actions[:,1].=.5f0;p.state.actions[:,2].=-.3f0
 w,p
end
function candidate_equal(a,b)
 before=Interior.device_cellular_snapshot(a);after=Interior.device_cellular_snapshot(b)
 all(isequal(getproperty(before.state,k),getproperty(after.state,k)) for k in keys(before.state))&&
 all(isequal(getproperty(before.environment,k),getproperty(after.environment,k)) for k in keys(before.environment))
end
@testset "Order-certified device contact candidates" begin
 for dense in (false,true),K in (1,64)
  _,a=candidate_fixture(;dense);b=deepcopy(a)
  cache=Interior.device_contact_candidate_cache(b;capacity=K,skin=.025f0)
  Interior.device_cellular_mechanics!(a,.0001f0,Val(true))
  Interior.device_candidate_mechanics!(b,.0001f0,cache)
  @test a.environment.status[1]==b.environment.status[1]==0
  @test candidate_equal(a,b)
  @test dense&&K==1 ? sum(cache.fallbacks)>0 : sum(cache.hits)>0
  @test dense&&K==1 ? any(>(0),cache.overflow) : true
 end
end
@testset "Coverage cache and explicit conservative invalidation" begin
 _,p=candidate_fixture();q=deepcopy(p)
 cache=Interior.device_contact_candidate_cache(q;capacity=64,skin=.025f0,strict_order=false)
 Interior.device_cellular_mechanics!(p,.0001f0,Val(true));Interior.device_candidate_mechanics!(q,.0001f0,cache)
 @test candidate_equal(p,q) # sparse fixture has no rounding-sensitive multi-contact sum
 @test cache.invalid[1]==0
 @test sum(cache.hits)>0
 # Check every cached/omitted pair against the full cube certificate.
 for i in eachindex(q.state.active)
  q.state.active[i]==1||continue
  selected=Set(cache.indices[i,1:cache.count[i]])
  for j in eachindex(q.state.active)
   q.state.active[j]==1&&i!=j||continue
   contact=Interior.dc_distance(q.state,q.geometry,i,j)<q.state.radii[i]+q.state.radii[j]
   @test !contact||j in selected
  end
 end
 q.state.position[1,1]+=cache.skin
 Interior.dc_contact_validate_at!(q.state,q.environment,q.geometry,cache,1)
 @test cache.invalid[1]==1
 cache.invalid[1]=0;q.state.radii[1]*=2
 Interior.dc_contact_validate_at!(q.state,q.environment,q.geometry,cache,1)
 @test cache.invalid[1]==1
 @test_throws ArgumentError Interior.device_contact_candidate_cache(q;skin=.2)
 @test_throws ArgumentError Interior.device_contact_candidate_cache(q;capacity=0)
end
@testset "Candidate overflow and exhausted skin have complete fallback" begin
 for dense in (false,true)
  _,a=candidate_fixture(;n=8,dense);b=deepcopy(a)
  # Tiny skin expires in the moving sparse case; capacity1 overflows in dense.
  cache=Interior.device_contact_candidate_cache(b;capacity=dense ? 1 : 64,skin=.0002f0,strict_order=false)
  a.state.actions[:,1].=1;b.state.actions[:,1].=1
  dt=dense ? .0001f0 : .01f0
  Interior.device_cellular_mechanics!(a,dt,Val(true));Interior.device_candidate_mechanics!(b,dt,cache)
  @test a.environment.status[1]==b.environment.status[1]==0
  @test sum(cache.fallbacks)>0
  @test maximum(abs.(a.state.position-b.state.position))<1e-5
  @test maximum(abs.(a.state.energy-b.state.energy))<1e-8
  @test a.state.bonds==b.state.bonds
  @test dense ? any(>(0),cache.overflow) : cache.invalid[1]==1
 end
end
@testset "Candidate mechanics composes with the complete cellular schedule" begin
 w,a=candidate_fixture(;n=8);w2=deepcopy(w);b=deepcopy(a)
 cache=Interior.device_contact_candidate_cache(b;strict_order=true)
 carrier=Interior.DeviceContactCandidateCarrier(b,cache)
 for _ in 1:3
  Interior.device_cellular_step!(w,a,.001)
  Interior.device_candidate_step!(w2,carrier,.001)
  @test candidate_equal(a,b)
  @test a.attempted_steps==b.attempted_steps
 end
end
struct CandidateScheduleTrace{P}
 process::P
 launches::Vector{Symbol}
end
Base.getproperty(p::CandidateScheduleTrace,key::Symbol)=key in (:process,:launches) ? getfield(p,key) : getproperty(getfield(p,:process),key)
Interior.device_cellular_launch!(p::CandidateScheduleTrace,k,args...;ndrange)=push!(p.launches,nameof(k))
@testset "Whole mechanical schedule removes launches only in coverage mode" begin
 _,p=candidate_fixture(;n=2)
 trace=CandidateScheduleTrace(p,Symbol[])
 Interior.device_cellular_mechanics!(trace,.001f0,Val(true));@test length(trace.launches)==82
 for strict_order in (false,true)
  empty!(trace.launches);cache=Interior.device_contact_candidate_cache(p;strict_order)
  Interior.device_candidate_mechanics!(trace,.001f0,cache)
  @test length(trace.launches)==(strict_order ? 83 : 53)
  @test count(==(:dc_cached_forces!),trace.launches)==16
 end
end
@testset "Skin query covers approaching pairs beyond old bins and periodic seam" begin
 for seam in (false,true)
  _,p=candidate_fixture(;n=2,density=30000.,bin_width=.1)
  s=p.state;g=p.geometry;e=p.environment
  s.position[:,:].=.5f0;s.actions[:,1:3].=0
  s.position[1,1]=seam ? .03f0 : .5f0;s.position[2,1]=seam ? 3.72f0 : .81f0
  cache=Interior.device_contact_candidate_cache(p;skin=.4f0,capacity=8,strict_order=false)
  Interior.device_candidate_mechanics!(p,.000001f0,cache)
  key=Interior.dc_bin(s,g,1);range=Interior.dc_search_range(s,e,g,1,1)
  oldbins=Set(Interior.dc_neighbor_bin(key[1],dx,g.bins[1],range) for dx in 0:Interior.dc_search_count(g.bins[1],range)-1)
  @test Interior.dc_bin(s,g,2)[1]∉oldbins
  @test 2 in cache.indices[1,1:cache.count[1]]
  shift=seam ? -.14f0 : .14f0
  s.position[1,1]=Interior.dc_wrap(s.position[1,1]+shift,g.extent[1])
  s.position[2,1]=Interior.dc_wrap(s.position[2,1]-shift,g.extent[1])
  for i in 1:2;Interior.dc_contact_validate_at!(s,e,g,cache,i);end
  @test cache.invalid[1]==0
  @test Interior.dc_distance(s,g,1,2)<s.radii[1]+s.radii[2]
  oracle=deepcopy(p)
  Interior.device_cellular_launch!(oracle,Interior.dc_clear_heads!,oracle.environment;ndrange=length(e.heads))
  Interior.device_cellular_launch!(oracle,Interior.dc_geometry!,oracle.state,oracle.environment,oracle.parameters,g;ndrange=length(s.active))
  Interior.device_cellular_launch!(oracle,Interior.dc_forces!,oracle.state,oracle.environment,oracle.parameters,g,.000001f0;ndrange=length(s.active))
  Interior.device_cellular_launch!(p,Interior.dc_cached_forces!,s,e,p.parameters,g,.000001f0,cache;ndrange=length(s.active))
  @test isequal(s.force,oracle.state.force)
  @test isequal(s.rows,oracle.state.rows)
 end
end
