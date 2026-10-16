using Test
include(joinpath(@__DIR__,"..","studies","brine_initialization_performance.jl"))
install_brine_initializer_reference!()
@testset "Specialized brine initialization preserves physical fields" begin
 for chemical_type in (Float32,Float64), phase in ("spinodal","lamellar"), reservoir in ("lens","sheet")
  c=TOML.parsefile(joinpath(@__DIR__,"..","experiments/art/brine-porous-384.toml"))
  c["phase_initialization"]=phase;c["reservoir_shape"]=reservoir
  phase=="lamellar"&&delete!(c,"phase_perturbation_rms")
  a=phase_runtime_fixture(c,9;chemical_type);b=phase_runtime_fixture(c,9;chemical_type)
  Base.invokelatest(brine_initialize_reference!,b.w,b.p,b.a,b.b,c)
  @test brine_initialization_equal(a,b)
 end
end
