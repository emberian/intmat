using Interior, Test, FFTW
include(joinpath(@__DIR__,"..","studies","implicit_phase_runtime.jl"))
@testset "External implicit phase runtime preserves physical continuation" begin
 w=World(Grid(9,7,5;extent=(9,14,7)));p=ImplicitPhaseSeparation(w;nutrient_diffusion=.01,activity_relaxation=.2)
 v,q=deepcopy((w,p))
 r=Interior.ImplicitPhaseRuntime(p)
 @test r.forward isa FFTW.FFTWPlan
 @test r!==Interior.ImplicitPhaseRuntime(q)
 for i in 1:10
  @. p.chemical_forcing=.1f0*p.phi
  @. q.chemical_forcing=.1f0*q.phi
  step!(w,ReusableImplicitPhase(p,r),.1);step!(v,UncachedImplicitPhase(q),.1)
 end
 @test p.phi==q.phi
 @test p.interface==q.interface
 @test p.chemical_potential==q.chemical_potential
 @test w.fields==v.fields
 io=IOBuffer();Interior.Serialization.serialize(io,(w,p,p));seekstart(io)
 rw,rp,alias=Interior.Serialization.deserialize(io)
 @test rp===alias
 @test Interior.ImplicitPhaseRuntime(rp)!==r
 step!(w,p,.1);step!(rw,rp,.1)
 @test p.phi==rp.phi
 @test w.fields==rw.fields

end

install_phase_coupling_reference!()
@testset "Zero diffusion preserves complete phase coupling" begin
 for D in (0.,.01), relaxation in (0.,.2)
  w=World(Grid(11,9,7));p=ImplicitPhaseSeparation(w;nutrient_diffusion=D,activity_relaxation=relaxation)
  v,q=deepcopy((w,p))
  for i in 1:4
   step!(w,p,.001);step!(v,ReferenceImplicitPhase(q,Interior.ImplicitPhaseRuntime(q)),.001)
  end
  @test p.phi==q.phi
  @test p.interface==q.interface
  @test w.fields==v.fields
  @test w.scratch==v.scratch
 end
 c=TOML.parsefile(joinpath(@__DIR__,"..","experiments/art/brine-lace-96.toml"))
 a=phase_runtime_fixture(c,16);b=phase_runtime_fixture(c,16)
 for i in 1:5
  coupled_runtime_step!(a,c,.01f0,true)
  coupled_runtime_step!(b,c,.01f0,true;reference=true)
 end
 @test a.p.phi==b.p.phi
 @test a.w.fields==b.w.fields
 @test a.a==b.a&&a.b==b.b
 @test a.p.interface==b.p.interface
end
