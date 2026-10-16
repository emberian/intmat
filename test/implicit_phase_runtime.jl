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

@testset "Implicit runtime layout contract" begin
 w=World(Grid(9,7,5));p=ImplicitPhaseSeparation(w)
 v=World(Grid(8,7,5));q=ImplicitPhaseSeparation(v)
 r=Interior.ImplicitPhaseRuntime(p)
 @test_throws ArgumentError Interior.advance!(v,q,.001,r)
 @test_throws ArgumentError Interior.advance!(w,p,-1)
 @test fieldnames(typeof(p))==(:reference,:chemical_forcing,:stabilization,:eigenvalues,:real_buffer,:spectrum,:nonlinear_spectrum)
 @test eltype(p.phi)===Float32&&eltype(p.real_buffer)===Float64&&eltype(p.spectrum)===ComplexF64
end

struct NoImplicitFloat64Backend <: Interior.KernelAbstractions.Backend end
Interior.KernelAbstractions.supports_float64(::NoImplicitFloat64Backend)=false
@testset "No silent implicit phase precision downgrade" begin
 cpu=World(Grid(3,3))
 unsupported=World(cpu.grid,NoImplicitFloat64Backend(),cpu.fields,cpu.scratch,cpu.time,cpu.seed)
 @test_throws ArgumentError ImplicitPhaseSeparation(unsupported)
end

@testset "Alias-preserving physical implicit phase checkpoints" begin
 w=World(Grid(7,5,3));p=ImplicitPhaseSeparation(w;noise=0,mean=.3,nutrient_diffusion=0,activity_relaxation=0)
 p=ImplicitPhaseSeparation(p.reference,p.phi,p.stabilization,p.eigenvalues,p.real_buffer,p.spectrum,p.nonlinear_spectrum)
 w=World(w.grid,w.backend,w.fields,(nutrient=w.fields.nutrient,),w.time,w.seed)
 cp=Interior.implicit_phase_checkpoint(w,p)
 @test cp.phase.chemical_forcing===cp.phase.phi
 @test cp.world.scratch.nutrient===cp.world.fields.nutrient
 @test cp.phase.phi!==p.phi
 io=IOBuffer();Interior.Serialization.serialize(io,(cp,cp.phase));seekstart(io)
 payload,alias=Interior.Serialization.deserialize(io)
 @test payload.phase===alias
 v,q=Interior.implicit_phase_restore(payload)
 @test q.chemical_forcing===q.phi
 @test v.scratch.nutrient===v.fields.nutrient
 @test q.real_buffer!==p.real_buffer
 q.chemical_forcing[1]=.4f0;p.chemical_forcing[1]=.4f0
 @test q.phi[1]==.4f0
 step!(w,p,.01);step!(v,q,.01)
 @test p.phi==q.phi&&w.fields==v.fields
 @test_throws ArgumentError Interior.implicit_phase_restore(merge(payload,(version=2,)))
end

@testset "Forcing shape is checked before transfer" begin
 w=World(Grid(9,7,5))
 @test_throws DimensionMismatch ImplicitPhaseSeparation(w;chemical_forcing=zeros(Float32,9,35))
 @test_throws DimensionMismatch ImplicitPhaseSeparation(w;chemical_forcing=zeros(Float32,15,7,3))
end

@testset "Explicit chemical storage reuses complete brine law" begin
 c=TOML.parsefile(joinpath(@__DIR__,"..","experiments/art/brine-lace-96.toml"))
 a=phase_runtime_fixture(c,8);b=phase_runtime_fixture(c,8;chemical_type=Float32)
 @test a.w.fields==b.w.fields&&a.a==b.a&&a.b==b.b
 @test_throws ArgumentError phase_runtime_fixture(c,8;chemical_type=Float16)
 d=phase_runtime_fixture(c,8;chemical_type=Float64)
 @test eltype(d.a)===eltype(d.b)===eltype(d.scratch)===eltype(d.w.fields.deposit)===Float64
 @test eltype(d.p.phi)===eltype(d.p.chemical_forcing)===eltype(d.permeability)===Float32
 @test d.a==Float64.(a.a)&&d.b==Float64.(a.b)&&d.p.phi==a.p.phi
 initial=sum(d.a)+sum(d.b)
 for _ in 1:20;coupled_runtime_step!(d,c,.01f0,true);end
 @test abs((sum(d.a)+sum(d.b)+sum(d.w.fields.deposit)-initial)/initial)<1e-12
 @test minimum(d.a)>=0&&minimum(d.b)>=0
end
