using Interior, Test, FFTW, Serialization
include(joinpath(@__DIR__,"..","studies","compact_brine_runtime.jl"))
FFTW.set_num_threads(2)
@testset "Compact brine exact physical trajectory and workspace lifetime" begin
 c=TOML.parsefile(joinpath(@__DIR__,"..","experiments/art/brine-porous-384-f64.toml"))
 for key in ("nutrient_diffusion","activity_relaxation","interface_activity","mineral_pinning","bulk")
  @test_throws ArgumentError compact_brine_fixture(merge(c,Dict(key=>.3)),4)
 end
 @test_throws ArgumentError compact_brine_fixture(merge(c,Dict("chemical_storage"=>"Float32")),4)
 for n in (9,16,24)
  ref=phase_runtime_fixture(c,n;chemical_type=Float64)
  x=compact_brine_fixture(c,n)
  dt=.0025f0
  @test x.p.phi==ref.p.phi && x.p.interface==ref.p.interface
  @test x.scratch===x.p.real_buffer
  @test propertynames(x.w.fields)==(:deposit,)
  for _ in 1:12
   coupled_runtime_step!(ref,c,dt,true);compact_brine_step!(x,c,dt)
   @test x.p.phi==ref.p.phi
   @test x.p.interface==ref.p.interface
   @test x.p.chemical_forcing==ref.p.chemical_forcing
   @test x.a==ref.a && x.b==ref.b && x.w.fields.deposit==ref.w.fields.deposit
   @test x.permeability==ref.permeability
   @test x.w.time==ref.w.time
  end
  checkpoint=compact_brine_checkpoint(x,c)
  io=IOBuffer();serialize(io,checkpoint);seekstart(io)
  y=compact_brine_restore(deserialize(io),c)
  @test y.scratch===y.p.real_buffer
  @test y.p.real_buffer!==x.p.real_buffer
  @test y.a==x.a && y.p.phi==x.p.phi && y.w.time==x.w.time
  for _ in 1:3;compact_brine_step!(x,c,dt);compact_brine_step!(y,c,dt);end
  @test y.a==x.a && y.b==x.b && y.w.fields.deposit==x.w.fields.deposit
  @test y.p.phi==x.p.phi && y.p.interface==x.p.interface && y.w.time==x.w.time
  @test_throws ArgumentError compact_brine_restore(merge(checkpoint,(version=2,)),c)
  @test_throws ArgumentError compact_brine_restore(merge(checkpoint,(source=Dict(),)),c)
  badstate=merge(checkpoint.state,(scratch=copy(checkpoint.state.scratch),))
  @test_throws ArgumentError compact_brine_restore(merge(checkpoint,(state=badstate,)),c)
  @test_throws ArgumentError compact_brine_restore(checkpoint,merge(c,Dict("reaction_rate"=>0.)))
  @test_throws ArgumentError compact_brine_step!(x,c,NaN)
  # Count owned arrays exactly once; runtime-native workareas are excluded.
  arrays=(x.p.phi,x.p.interface,x.p.chemical_forcing,x.p.real_buffer,x.p.spectrum,x.p.nonlinear_spectrum,x.p.eigen_axes...,x.a,x.b,x.w.fields.deposit,x.permeability)
  @test sum(sizeof,arrays)==compact_brine_memory(n).owned_bytes
 end
end
