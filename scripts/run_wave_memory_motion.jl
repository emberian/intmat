include(joinpath(@__DIR__,"../studies/wave_memory.jl"))
kind=get(ENV,"WAVE_BACKEND","cpu")
if kind=="metal";using Metal;backend=Metal.MetalBackend()
elseif kind=="amd";using AMDGPU;backend=AMDGPU.ROCBackend()
else;backend=CPU();end
function produce(out;backend,n=128,steps=600,cadence=15)
 ispath(out)&&error("fresh output required");mkpath(out)
 dt=.006f0;w,p,m=packet_memory_fixture(backend,n)
 hashes=Dict(f=>bytes2hex(sha256(read(joinpath(@__DIR__,"..",f)))) for f in ("studies/wave_memory.jl","src/acoustics.jl","scripts/run_wave_memory_motion.jl"))
 frames=Any[]
 function capture(k)
  path=joinpath(out,"frames",lpad(k,5,'0'));d=Array(w.fields.deposit);host=Array(m.host);pr=Array(p.pressure)
  export_volume((deposit=d,density=host,activity=abs.(pr),envelope=Array(p.envelope)),path)
  push!(frames,(step=k,time=w.time,directory=relpath(path,out),mineral_extrema=extrema(d),mineral_sum=sum(Float64,d),pressure_extrema=extrema(pr),inventory_error=maximum(abs.(d.+Array(w.fields.nutrient).-host))))
  println(last(frames));flush(stdout)
 end
 capture(0);elapsed=@elapsed for k=1:steps
  wave_memory_step!(w,p,m,dt)
  (k%cadence==0||k==steps)&&capture(k)
 end
 meta=(title="Where the note goes",shape=(n,n,n),steps=steps,dt=dt,elapsed_including_capture_seconds=elapsed,backend=kind,frames=frames,source_sha256=hashes,forcing="Three authored broad coherent packets, prepared at zero velocity with production face operator; no ongoing source. Camera does not move.",model="Shared AcousticField plus finite exponential precursor/mineral conversion and independent.25-rate dissolution; both donors bounded. Phenomenological conversion, no latent or acoustic energy debit.",actual_states=true,temporal_interpolation=false)
 open(io->JSON3.pretty(io,meta),joinpath(out,"run.json"),"w")
end
produce(ARGS[1];backend,n=parse(Int,get(ENV,"WAVE_N","128")),steps=parse(Int,get(ENV,"WAVE_STEPS","600")),cadence=parse(Int,get(ENV,"WAVE_CADENCE","15")))
