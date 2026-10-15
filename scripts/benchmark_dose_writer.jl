# Paired warm kernel timings, no rendering; original and dose-only writer.
using Interior,Statistics,Metal
using Interior: JSON3
sizes=isempty(ARGS) ? (96,384) : parse.(Int,split(ARGS[1],",")); all(n->n<=384,sizes)||error("bounded benchmark only")
rows=Any[]
for n in sizes
    backend=Metal.MetalBackend();g=Grid(n,n,n)
    nutrient=Interior.KernelAbstractions.zeros(backend,Float32,g.size...);fill!(nutrient,1f0)
    deposit=similar(nutrient);fill!(deposit,0f0)
    a=World(g,backend,(;nutrient,deposit),NamedTuple(),0.,91);b=World(g,backend,NamedTuple(),NamedTuple(),0.,91)
    old=ThermalWriter(a;width=.018,heat=30,threshold=.05,conversion=18);new=DoseThermalWriter(b;width=.018,heat=30,threshold=.05,conversion=18)
    dt=min(.003125,Float64(stability_limit(a,old))*.9)
    for _ in 1:3;step!(a,old,dt);step!(b,new,dt);end
    x=Float64[];y=Float64[]
    # Alternate order to reduce slow temperature/device scheduling bias.
    for i in 1:30
        if isodd(i);push!(x,@elapsed step!(a,old,dt));push!(y,@elapsed step!(b,new,dt))
        else;push!(y,@elapsed step!(b,new,dt));push!(x,@elapsed step!(a,old,dt));end
    end
    difference=maximum(abs.(Array(old.temperature).-Array(new.writer.temperature)))
    @assert difference<=2f-5
    push!(rows,(;n,dt,original_seconds=median(x),dose_pingpong_seconds=median(y),ratio=median(x)/median(y),temperature_error=difference,original_bytes=20*n^3,dose_bytes=12*n^3))
end
result=(;schema="dose-writer-metal-warm",shared_load=true,rows,sources=Interior.arena_sources())
println(JSON3.write(result));length(ARGS)>1&&open(io->JSON3.pretty(io,result),ARGS[2],"w")
