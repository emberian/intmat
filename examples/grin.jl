using Interior,JSON3,SHA
out=isempty(ARGS) ? joinpath(@__DIR__,"..","outputs","grin") : abspath(ARGS[1]);mkpath(out)
N=64;g=Grid(N,N,N);w=World(g)
writer=ThermalWriter(w;path=:lissajous,width=.12,heat=5,cooling=.6,diffusion=.002,stop=1.4,period=2,conversion=1)
dt=.01
for _ in 1:140;step!(w,writer,dt);end
T=Array(writer.temperature);peak=maximum(T)
index=ones(Float32,N,N,N)
for k in 1:N,j in 1:N,i in 1:N
    x,y,z=ntuple(d->((i,j,k)[d]-.5)*2/N-1,3)
    taper=prod(max(0.,1-(v/.95)^2)^2 for v in (x,y,z))
    # Deliberately exaggerated positive thermo-optic response, not calibrated glass.
    index[i,j,k]=1+.2*T[i,j,k]/peak*taper
end
f=GradedIndex(index);constant=GradedIndex(ones(Float32,N,N,N))
Interior.export_volume((temperature=T,refractive_index=index,deposit=Array(w.fields.deposit)),out)
settings=(width=320,height=320,ds=.012,background_z=1.7,span=2.4)
plain=render_grin(constant;settings...);curved=render_grin(f;settings...)
Interior.write_ppm(joinpath(out,"constant-index.ppm"),plain.rgb;exposure=1)
Interior.write_ppm(joinpath(out,"thermal-index.ppm"),curved.rgb;exposure=1)
rays=[trace_grin(f,(x,0.,-1.),(0.,0.,1.);ds=.008) for x in -.8:.1:.8]
refined=render_grin(f;width=64,height=64,ds=.006,background_z=1.7,span=2.4)
coarse=render_grin(f;width=64,height=64,ds=.012,background_z=1.7,span=2.4)
finer=render_grin(f;width=64,height=64,ds=.003,background_z=1.7,span=2.4)
maxdifference(a,b)=maximum(abs.(a.hits.-b.hits))
record=(state=(source="ThermalWriter lissajous, 64^3, 140 steps dt=.01",time=w.time,temperature_peak=peak,index_range=extrema(index),index_law="1 + .2 temperature/peak * boundary taper; invented exaggerated response",index_sha256=bytes2hex(sha256(reinterpret(UInt8,vec(index))))),observer=(camera="orthographic toward +z",settings=settings,light="fixed luminous background chart",model="midpoint unit-direction eikonal rays; no scattering/interfaces"),checks=(unresolved=curved.unresolved,coarse_refined_max_hit_difference=maxdifference(coarse,refined),refined_finer_max_hit_difference=maxdifference(refined,finer),constant_variable_max_hit_difference=maxdifference(plain,curved)),rays=rays)
open(joinpath(out,"observation.json"),"w") do io;JSON3.write(io,record);end
println(JSON3.write(record.checks))
