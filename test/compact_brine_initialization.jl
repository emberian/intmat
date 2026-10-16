using Interior, Test, FFTW, Serialization
isdefined(@__MODULE__,:compact_brine_seed)||include(joinpath(@__DIR__,"..","studies","compact_brine_runtime.jl"))
FFTW.set_num_threads(2)
@testset "Data-only brine seed owns no FFT resources" begin
 c=TOML.parsefile(joinpath(@__DIR__,"..","experiments/art/brine-porous-384-f64.toml"))
 for n in (9,16,24)
  seed=compact_brine_seed(c,n)
  @test propertynames(seed.p)==(:phi,:interface)
  @test propertynames(seed)==(:w,:p,:a,:b,:recipe)
  @test sum(sizeof,(seed.p.phi,seed.p.interface,seed.a,seed.b,seed.w.fields.deposit))==compact_brine_seed_memory(n)
  x=compact_brine_seed_upload(seed,c,Interior.CPU())
  ref=compact_brine_fixture(c,n)
  @test x.p.phi==ref.p.phi && x.p.interface==ref.p.interface
  @test x.a==ref.a && x.b==ref.b && x.w.fields.deposit==ref.w.fields.deposit
  @test x.p.chemical_forcing==ref.p.chemical_forcing && x.permeability==ref.permeability
  @test x.scratch===x.p.real_buffer
  @test x.p.phi!==seed.p.phi && x.a!==seed.a && x.w.fields.deposit!==seed.w.fields.deposit
  for _ in 1:5;compact_brine_step!(x,c,.0025f0);compact_brine_step!(ref,c,.0025f0);end
  @test x.p.phi==ref.p.phi && x.p.interface==ref.p.interface
  @test x.a==ref.a && x.b==ref.b && x.w.fields.deposit==ref.w.fields.deposit && x.w.time==ref.w.time
  @test_throws ArgumentError compact_brine_seed_upload(seed,merge(c,Dict("reaction_rate"=>0.)),Interior.CPU())
  @test_throws ArgumentError compact_brine_seed_upload(merge(seed,(a=Float32.(seed.a),)),c,Interior.CPU())
  badtime=deepcopy(seed);badtime.w.time=1
  @test_throws ArgumentError compact_brine_seed_upload(badtime,c,Interior.CPU())
  io=IOBuffer();serialize(io,seed);seekstart(io)
  y=compact_brine_seed_upload(deserialize(io),c,Interior.CPU())
  @test y.p.phi==seed.p.phi && y.a==seed.a && y.scratch===y.p.real_buffer
 end
 @test compact_brine_seed_memory(512)==4*1024^3
end
