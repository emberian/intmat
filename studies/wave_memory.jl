# GPU-native wave-specific finite material memory. Shared AcousticField stays unchanged.
using Interior, SHA, Statistics
using Interior: JSON3, KernelAbstractions
using Interior.KernelAbstractions: @kernel, @index, CPU
import Interior: advance!,validate,stability_limit
struct WaveMemory{A}
 host::A
 dose::A
 envelope::A
 gain::Float32
 threshold::Float32
 relaxation::Float32
end
function WaveMemory(w,p;host,gain=12,threshold=.0004,relaxation=0)
 h=KernelAbstractions.zeros(w.backend,Float32,w.grid.size...);copyto!(h,host)
 WaveMemory(h,KernelAbstractions.zeros(w.backend,Float32,w.grid.size...),p.envelope,Float32(gain),Float32(threshold),Float32(relaxation))
end
function validate(w,p::WaveMemory)
 all(isfinite,(p.gain,p.threshold,p.relaxation))&&min(p.gain,p.threshold,p.relaxation)>=0||error("invalid memory parameters")
 size(p.host)==size(w.fields.deposit)==size(p.envelope)||error("memory shape mismatch")
end
@inline function wave_transfer(m,host,intensity,dt,gain,threshold,relaxation)
 # Both transfers are finite fractions of their respective donor inventories.
 put=(host-m)*(-expm1(-dt*gain*max(0f0,intensity-threshold)))
 take=m*(-expm1(-dt*relaxation))
 m+put-take
end
@kernel function wave_memory_kernel!(mineral,nutrient,dose,host,envelope,dt,gain,threshold,relaxation)
 i=@index(Global,Linear)
 intensity=envelope[i];dose[i]+=dt*intensity
 mineral[i]=wave_transfer(mineral[i],host[i],intensity,dt,gain,threshold,relaxation)
 nutrient[i]=host[i]-mineral[i]
end
function advance!(w::World,p::WaveMemory,dt::Real)
 Interior.launch!(w,wave_memory_kernel!,w.fields.deposit,w.fields.nutrient,p.dose,p.host,p.envelope,Float32(dt),p.gain,p.threshold,p.relaxation)
end
function memory_fixture(backend,n;gain=12,threshold=.0004,relaxation=0)
 w=World(Grid(n,n,n);backend,fields=(:nutrient,:deposit,:activity,:pressure),scratch_fields=())
 p=AcousticField(w;speed=.44,material_feedback=12,damping=.10,source_amplitude=30,frequency=2.2,source_radius=.08,source_center=(-.42,0,0),activity_gain=0,envelope_time=.12)
 host=[Float32(clamp((1-sqrt(((2(i-.5)/n-1)/.86)^2+((2(j-.5)/n-1)/.63)^2+((2(k-.5)/n-1)/.83)^2))/.065,0,1)) for i=1:n,j=1:n,k=1:n]
 m=WaveMemory(w,p;host,gain,threshold,relaxation);copyto!(w.fields.deposit,.015f0.*host);copyto!(w.fields.nutrient,.985f0.*host)
 w,p,m
end
function scalar_memory!(w,m,dt)
 dt=Float32(dt)
 for i in eachindex(m.host)
  m.dose[i]+=dt*m.envelope[i]
  w.fields.deposit[i]=wave_transfer(w.fields.deposit[i],m.host[i],m.envelope[i],dt,m.gain,m.threshold,m.relaxation)
  w.fields.nutrient[i]=m.host[i]-w.fields.deposit[i]
 end
end
function wave_memory_step!(w,p,m,dt;scalar=false)
 Interior.advance!(w,p,dt)
 scalar ? scalar_memory!(w,m,dt) : Interior.advance!(w,m,dt)
 w.time+=Float64(Float32(dt))
end
function wave_memory_fields(w,p,m)
 (;deposit=Array(w.fields.deposit),nutrient=Array(w.fields.nutrient),pressure=Array(p.pressure),previous=Array(p.previous),envelope=Array(p.envelope),dose=Array(m.dose))
end

"""Authored broad coherent packets; zero-velocity setup uses production initializer."""
function packet_memory_fixture(backend,n;relaxation=.25)
 w,p,m=memory_fixture(backend,n;gain=5,threshold=.018,relaxation)
 p.source_amplitude=0;p.damping=.065f0;p.envelope_time=.18f0
 pressure=zeros(Float32,n,n,n)
 for k=1:n,j=1:n,i=1:n
  x,y,z=ntuple(d->2*((i,j,k)[d]-.5)/n-1,3);q=0.
  for (u,v,t,c,phase) in ((1.,.32,.22,-.34,0.),(-.27,1.,.41,.29,1.1),(.35,-.24,1.,-.18,2.2))
   d=(u*x+v*y+t*z)/sqrt(u*u+v*v+t*t)
   q+=.7cos(34d+phase)*exp(-((d-c)/.19)^2)*exp(-((x/.86)^8+(y/.86)^8+(z/.86)^8))
  end
  pressure[i,j,k]=q
 end
 Interior.seed_acoustic_field!(w,p,pressure;dt=.006f0)
 w,p,m
end
