# Experimental brine-only representation. No physical law is replaced here.
# Include after implicit_phase_runtime.jl; runtime resources are never serialized.
isdefined(@__MODULE__, :phase_runtime_fixture) || include("implicit_phase_runtime.jl")
function compact_brine_contract(c)
 for key in ("nutrient_diffusion","activity_relaxation","interface_activity","mineral_pinning")
  get(c,key,0)==0 || throw(ArgumentError("compact brine requires $key=0; use the ordinary phase runtime for this coupling"))
 end
 get(c,"bulk",1)==1 || throw(ArgumentError("compact brine retains the original recipe bulk=1"))
 get(c,"chemical_storage","Float64")=="Float64" || throw(ArgumentError("compact brine requires Float64 chemical storage"))
 nothing
end
struct CompactBrinePhase{A,R,C,E}
 phi::A
 interface::A
 chemical_forcing::A
 real_buffer::R
 spectrum::C
 nonlinear_spectrum::C
 eigen_axes::E
 bulk::Float32
 kappa::Float32
 mobility::Float32
 stabilization::Float64
end
function compact_brine_phase(grid,backend,c)
 z(T,shape)=Interior.KernelAbstractions.zeros(backend,T,shape...)
 dims=grid.size; hs=(dims[1]÷2+1,dims[2],dims[3])
 axes=ntuple(3) do d
  h=[dims[d]>1 ? 4sinpi((i-1)/dims[d])^2/grid.spacing[d]^2 : 0.0 for i in 1:hs[d]]
  a=z(Float64,(length(h),));copyto!(a,h);a
 end
 all(v->isfinite(v)&&v>=0,(get(c,"bulk",1),c["kappa"],c["mobility"],c["stabilization"])) && get(c,"bulk",1)>0 && c["kappa"]>0 || throw(ArgumentError("invalid phase coefficients"))
 CompactBrinePhase(z(Float32,dims),z(Float32,dims),z(Float32,dims),z(Float64,dims),z(ComplexF64,hs),z(ComplexF64,hs),axes,
  Float32(get(c,"bulk",1)),Float32(c["kappa"]),Float32(c["mobility"]),Float64(c["stabilization"]))
end
@kernel function compact_brine_interface_kernel!(interface,phi,h)
 i,j,k=@index(Global,NTuple)
 interface[i,j,k]=Interior.phase_interface_value(phi,(i,j,k),h)
end
Interior.phase_observe!(w,p::CompactBrinePhase)=Interior.launch!(w,compact_brine_interface_kernel!,p.interface,p.phi,w.grid.spacing)
@kernel function compact_brine_modes_kernel!(spectrum,nonlinear,axes,dt,m,s,kappa)
 i,j,k=@index(Global,NTuple)
 q=Float64(Float32(axes[1][i]+axes[2][j]+axes[3][k]))
 spectrum[i,j,k]=Interior.implicit_phase_mode(spectrum[i,j,k],nonlinear[i,j,k],q,dt,m,s,kappa)
end
function compact_brine_fixture(c,n;backend=Interior.CPU())
 compact_brine_contract(c)
 backend isa Interior.CPU || error("Initialize on CPU, then compact_brine_transfer; native initializer is a separate gate")
 grid=Grid(n,n,n;extent=ntuple(_->c["extent"],3))
 mineral=zeros(Float64,grid.size)
 w=World(grid,backend,(deposit=mineral,),NamedTuple(),0.0,17)
 p=compact_brine_phase(grid,backend,c)
 a=zeros(Float64,grid.size);b=similar(a);permeability=zeros(Float32,grid.size)
 brine_initialize!(w,p,a,b,c)
 (;w,p,a,b,scratch=p.real_buffer,permeability,runtime=Interior.implicit_phase_runtime(p))
end
function compact_brine_step!(x,c,dt)
 compact_brine_contract(c)
 (;w,p,a,b,scratch,permeability,runtime)=x
 isfinite(dt)&&dt>0&&isfinite(Float32(dt))&&Float32(dt)>0 || throw(ArgumentError("positive finite Float32 timestep required"))
 dt=Float32(dt)
 # Same explicit solute positivity constraint as the original brine producer.
 dt*6max(Float32(c["diffusion_a"]),Float32(c["diffusion_b"]))/minimum(w.grid.spacing)^2<.9 || throw(ArgumentError("brine diffusion positivity bound"))
 runtime.layout==Interior.implicit_phase_layout(p) || throw(ArgumentError("compact runtime layout differs"))
 scratch===p.real_buffer || throw(ArgumentError("compact transport scratch must share dead FFT real workspace"))
 feedback=Float32(c["mineral_phase_feedback"])
 @. p.chemical_forcing=feedback*w.fields.deposit*p.phi
 all(isfinite,p.phi) && all(isfinite,p.chemical_forcing) || throw(ArgumentError("phase composition and forcing must be finite"))
 Interior.implicit_phase_real!(w,p,runtime)
 mul!(p.spectrum,runtime.forward,p.real_buffer)
 Interior.implicit_phase_nonlinear!(w,p,runtime)
 mul!(p.nonlinear_spectrum,runtime.forward,p.real_buffer)
 compact_brine_modes_kernel!(w.backend,256)(p.spectrum,p.nonlinear_spectrum,p.eigen_axes,dt,Float64(p.mobility),p.stabilization,p.kappa;ndrange=size(p.spectrum))
 Interior.implicit_phase_inverse!(w,p,runtime)
 Interior.phase_observe!(w,p)
 w.time+=Float64(dt)
 phase_barrier=Float32(c["phase_transport_barrier"]);mineral_barrier=Float32(c["mineral_transport_barrier"])
 @. permeability=1f0/(1f0+phase_barrier*(1f0+tanh(2f0*p.phi))/2f0+mineral_barrier*w.fields.deposit)
 Interior.medium_diffuse!(scratch,a,permeability,w.grid,Float32(c["diffusion_a"]),dt)
 Interior.medium_diffuse!(scratch,b,permeability,w.grid,Float32(c["diffusion_b"]),dt)
 Interior.launch!(w,brine_react!,a,b,w.fields.deposit,p.interface,dt,Float32(c["reaction_rate"]),Float32(c["basal_reaction"]),Float32(c["nucleation_gradient"]),Float32(c["saturation_gradient"]))
 x
end
# One memo keeps aliases across state and phase; FFT resources rebuilt afterward.
function compact_brine_transfer(x,backend)
 memo=IdDict()
 transfer(a)=get!(memo,a) do
  b=Interior.KernelAbstractions.zeros(backend,eltype(a),size(a)...);copyto!(b,a);b
 end
 w=World(x.w.grid,backend,map(transfer,x.w.fields),NamedTuple(),x.w.time,x.w.seed)
 p=x.p
 q=CompactBrinePhase(transfer(p.phi),transfer(p.interface),transfer(p.chemical_forcing),transfer(p.real_buffer),transfer(p.spectrum),transfer(p.nonlinear_spectrum),map(transfer,p.eigen_axes),p.bulk,p.kappa,p.mobility,p.stabilization)
 (;w,p=q,a=transfer(x.a),b=transfer(x.b),scratch=transfer(x.scratch),permeability=transfer(x.permeability),runtime=Interior.implicit_phase_runtime(q))
end
# Existing default remains unchanged. This forecast excludes driver baseline and
# unmeasured size-specific FFT workareas; it is never permission to allocate.
function compact_brine_memory(n)
 N=n^3;H=(n÷2+1)*n^2;axes=8*(n÷2+1+2n)
 (;n,owned_bytes=48N+32H+axes,native_active_bytes=48N+64H+axes,full_four_field_bytes=28N)
end
const COMPACT_BRINE_SOURCE=Dict(f=>bytes2hex(SHA.sha256(read(joinpath(@__DIR__,"..",f)))) for f in
 ("src/core.jl","src/phase.jl","src/implicit_phase.jl","src/compartments.jl","studies/brine_memory.jl","studies/implicit_phase_runtime.jl","studies/compact_brine_runtime.jl"))
function compact_brine_checkpoint(x,c)
 host=compact_brine_transfer(x,Interior.CPU())
 # Runtime plans have deliberately no route into the payload.
 state=(;host.w,host.p,host.a,host.b,host.scratch,host.permeability)
 (;format="Interior.CompactBrine",version=1,julia_version=VERSION,source=copy(COMPACT_BRINE_SOURCE),recipe=deepcopy(c),state)
end
function compact_brine_restore(payload,c;backend=Interior.CPU())
 compact_brine_contract(c)
 payload.format=="Interior.CompactBrine" && payload.version==1 && payload.julia_version==VERSION && payload.source==COMPACT_BRINE_SOURCE && payload.recipe==c || throw(ArgumentError("compact brine checkpoint source/schema/recipe differs"))
 payload.state.scratch===payload.state.p.real_buffer || throw(ArgumentError("compact checkpoint scratch alias differs"))
 compact_brine_transfer(payload.state,backend)
end

# Initialization-only host representation: no spectral buffers or FFT plans.
# This is not a resumable compact state; initialization happens at time zero.
struct CompactBrineSeed{A}
 phi::A
 interface::A
end
Interior.phase_observe!(w,p::CompactBrineSeed)=Interior.launch!(w,compact_brine_interface_kernel!,p.interface,p.phi,w.grid.spacing)
function compact_brine_seed(c,n)
 compact_brine_contract(c)
 grid=Grid(n,n,n;extent=ntuple(_->c["extent"],3));backend=Interior.CPU()
 w=World(grid,backend,(deposit=zeros(Float64,grid.size),),NamedTuple(),0.0,17)
 p=CompactBrineSeed(zeros(Float32,grid.size),zeros(Float32,grid.size))
 a=zeros(Float64,grid.size);b=similar(a)
 brine_initialize!(w,p,a,b,c)
 (;w,p,a,b,recipe=deepcopy(c))
end
function compact_brine_seed_upload(seed,c,backend)
 compact_brine_contract(c)
 seed.recipe==c || throw(ArgumentError("compact brine seed recipe differs"))
 seed.w.time==0 || throw(ArgumentError("compact brine seed is initialization-only"))
 Interior.KernelAbstractions.supports_float64(backend) || throw(ArgumentError("compact brine requires Float64 arithmetic and FFT support"))
 dims=seed.w.grid.size
 all(size(a)==dims for a in (seed.p.phi,seed.p.interface,seed.a,seed.b,seed.w.fields.deposit)) || throw(DimensionMismatch("compact seed grid differs"))
 upload(a)=begin b=Interior.KernelAbstractions.zeros(backend,eltype(a),size(a)...);copyto!(b,a);b end
 w=World(seed.w.grid,backend,(deposit=upload(seed.w.fields.deposit),),NamedTuple(),0.0,seed.w.seed)
 p=compact_brine_phase(w.grid,backend,c)
 copyto!(p.phi,seed.p.phi);copyto!(p.interface,seed.p.interface)
 a=upload(seed.a);b=upload(seed.b)
 permeability=Interior.KernelAbstractions.zeros(backend,Float32,dims...)
 (;w,p,a,b,scratch=p.real_buffer,permeability,runtime=Interior.implicit_phase_runtime(p))
end
compact_brine_seed_memory(n)=32n^3
