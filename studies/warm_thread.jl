# Continuous, mutually inventory-limited thermal knots inside a changing acoustic medium.
using Interior, Metal, SHA, Statistics
using Interior: JSON3
struct KnotWriter{W}
 writer::W
 kind::Int
end
function knot_position(p::KnotWriter,t)
 q=p.writer;phase=2pi*(t-q.start)/q.period
 if p.kind==1
  r=.38+.125cos(7phase);return Float32.((r*cos(3phase),r*sin(3phase),.47sin(7phase)))
 end
 r=.32+.12cos(5phase)
 Float32.((.46sin(5phase),r*cos(2phase),r*sin(2phase)))
end
Interior.validate(w,p::KnotWriter)=Interior.validate(w,p.writer)
Interior.stability_limit(w,p::KnotWriter)=Interior.stability_limit(w,p.writer)
function Interior.advance!(w::World,p::KnotWriter,dt::Real)
 q=p.writer;t=w.time+Float64(dt)/2;power=q.start<=t<q.stop ? q.heat : 0f0
 Interior.launch!(w,Interior.write_heat_kernel!,q.scratch,q.temperature,q.dose,w.fields.nutrient,w.fields.deposit,w.grid.spacing,w.grid.extent,Float32(dt),knot_position(p,t),q.width,power,q.diffusion,q.cooling,q.threshold,q.conversion)
 copyto!(q.temperature,q.scratch);w
end
function warm_thread(out;n=384,steps=2500,dt=.0032)
 estimate=4*28*n^3;estimate<=parse(Int,get(ENV,"WARM_MEMORY_BUDGET_MB","8192"))*1024^2 || error("Estimated resident plus capture memory exceeds budget")
 ispath(out)&&error("Fresh output required");mkpath(out)
 sources=("studies/warm_thread.jl","src/writing.jl","src/acoustics.jl")
 for f in sources;cp(f,joinpath(out,basename(f)));end
 w=World(Grid(n,n,n);seed=73,backend=Metal.MetalBackend(),fields=(:nutrient,:deposit,:activity,:pressure),scratch_fields=())
 writers=(KnotWriter(ThermalWriter(w;period=7,width=.026,heat=140,diffusion=.00012,cooling=8,threshold=.035,conversion=14,stop=7),1),KnotWriter(ThermalWriter(w;period=6.8,width=.026,heat=125,diffusion=.00012,cooling=8,threshold=.035,conversion=12,start=.45,stop=7.25),2))
 acoustic=AcousticField(w;speed=.37,material_feedback=14,damping=.13,source_amplitude=30,frequency=1.7,source_radius=.07,source_center=(-.34,-.1,-.15),activity_gain=0)
 processes=(writers...,acoustic);@assert dt<=minimum(stability_limit(w,p) for p in processes)
 maxjumps=[maximum(sqrt(sum((knot_position(p,t+dt)[d]-knot_position(p,t)[d])^2 for d in 1:3)) for t in 0:dt:steps*dt) for p in writers]
 maximum(maxjumps)<.5*.026 || error("Toolpath sampling coarser than half a source sigma")
 original=n^3;frames=Any[];times=Float64[]
 for step in 1:steps
  elapsed=@elapsed step!(w,processes,dt);step>1&&push!(times,elapsed)
  if step in (round(Int,steps/3),round(Int,2steps/3),steps)
   directory=joinpath(out,"frames",lpad(step,5,'0'));mkpath(directory)
   temp=Array(writers[1].writer.temperature).+Array(writers[2].writer.temperature);dose=Array(writers[1].writer.dose).+Array(writers[2].writer.dose)
   solid=Array(w.fields.deposit);food=Array(w.fields.nutrient);pressure=Array(acoustic.pressure)
   d1=-expm1.(-writers[1].writer.conversion.*Array(writers[1].writer.dose));d2=-expm1.(-writers[2].writer.conversion.*Array(writers[2].writer.dose))
   overlap=sum(Float64,d1.*d2);composition_error=maximum(abs.(solid.-(d1.+d2.-d1.*d2)))
   @assert all(isfinite,temp)&&all(isfinite,pressure)&&minimum(food)>=0
   residual=(sum(Float64,solid)+sum(Float64,food)-original)/original;@assert abs(residual)<1e-6&&composition_error<2e-5
   for (name,field) in (("temperature",temp),("dose",dose),("solid",solid),("pressure",pressure));write(joinpath(directory,name*".f32"),field);end
   push!(frames,(step=step,time=w.time,directory=relpath(directory,out),mineral_max=maximum(solid),mineral_sum=sum(Float64,solid),temperature_max=maximum(temp),pressure_max=maximum(abs,pressure),shared_inventory_suppression=overlap,shared_composition_max_error=composition_error,material_relative_error=residual))
   println(JSON3.write(last(frames)));flush(stdout)
  elseif step%100==0
   println("STEP ",step," median=",median(times));flush(stdout)
  end
 end
 report=(id="warm-thread",title="The warm thread remembers",description="Two continuous thermal knot histories braid mineral into an amber medium; sound crosses the material they leave behind.",shape=(n,n,n),extent=(2,2,2),dt=dt,steps=steps,time=w.time,frames=frames,median_step_seconds=median(times),evolution_seconds=sum(times),spacing=2/n,width_cells=.026*n/2,maximum_source_step_distance=maxjumps,thermal_stability_limits=[stability_limit(w,p) for p in writers],acoustic_stability_limit=stability_limit(w,acoustic),source_sha256=Dict(f=>bytes2hex(sha256(read(joinpath(out,basename(f))))) for f in sources),simulation="Existing Julia write_heat_kernel!/thermal_exposure ×2 conservatively convert shared nutrient into deposit using accumulated supra-threshold exposure. New source-position functions prescribe continuous3:7 and rotated2:5 thermal knots, not emergent motion. AcousticField propagates a damped scalar wave with deposit-dependent face wave speed. Paths are designed, mineral and pressure are solved; model units, not a calibrated physical specimen.",writers=[(kind=p.kind,start=p.writer.start,stop=p.writer.stop,period=p.writer.period,width=p.writer.width,heat=p.writer.heat,diffusion=p.writer.diffusion,cooling=p.writer.cooling,threshold=p.writer.threshold,conversion=p.writer.conversion) for p in writers])
 open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w")
end
if abspath(PROGRAM_FILE)==@__FILE__
 warm_thread(ARGS[1];n=parse(Int,get(ENV,"WARM_N","384")),steps=parse(Int,get(ENV,"WARM_STEPS","2500")))
end
