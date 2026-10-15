using Interior,Test,Random
include(joinpath(@__DIR__,"..","studies","lattice_split_bank.jl"))
mods=install_lattice_split_variants(mktempdir())
Base.invokelatest() do
 @testset "Split finite genotype coefficients and ownership" begin
  for loci in (1,6,12)
   audit=lattice_split_coefficient_audit(mods.split;loci)
   @test audit.violations==0
   @test audit.maximum_delta<=audit.conservative_delta_bound
   @test audit.new_error<2e-6
  end
  bank=mods.split.lattice_splitbank()
  @test sizeof(bank)==2*64*14*11*4
  @test sizeof(bank)<sizeof(LatticeEcologyStudy.lattice_genotypes())/30
  @test_throws BoundsError bank[0,1,1]
  @test_throws BoundsError bank[4097,1,1]
  w=LatticeEcologyStudy.lattice_world(Grid(8,8;extent=(8,8,1)))
  p=mods.split.LatticeEcology(w)
  @test Interior.validate(w,p)===nothing
  @test !Base.mightalias(p.weights,p.cells)
  @test_throws ArgumentError mods.split.lattice_checkpoint(w,p,joinpath(mktempdir(),"unsupported.jls"))
 end
end
Base.invokelatest() do
 @testset "Independent bounded controller dot error" begin
  audit=lattice_split_controller_audit(mods.split)
  @test audit.violations==0
  @test audit.maximum_dot_difference<=audit.conservative_dot_bound
  @test audit.hidden_error<2e-6
  @test audit.action_error<2e-6
  @test audit.neurons==256*11
 end
 @testset "Zero-strength split and packed coefficients" begin
  old=LatticeEcologyStudy.lattice_genotypes(;loci=6,strength=0)
  bank=mods.split.lattice_splitbank(;loci=6,strength=0)
  @test isequal(old,Array(bank))
  w=LatticeEcologyStudy.lattice_world(Grid(8,8;extent=(8,8,1)))
  a=LatticeEcologyStudy.LatticeEcology(w;loci=6)
  b=mods.split.LatticeEcology(deepcopy(w);loci=6)
  initiala=LatticeEcologyStudy.lattice_statistics(w,a);initialb=mods.split.lattice_statistics(w,b)
  @test initiala.material==initialb.material
  @test initiala.energy_accounted==initialb.energy_accounted
  for _ in 1:12;Interior.step!(w,a,.1);end
  w2=LatticeEcologyStudy.lattice_world(Grid(8,8;extent=(8,8,1)))
  for _ in 1:12;Interior.step!(w2,b,.1);end
  after=mods.split.lattice_statistics(w2,b)
  @test after.finite
  @test abs(after.material-initialb.material)<1e-4
  @test abs(after.energy_accounted-initialb.energy_accounted)<1e-4
 end
end
