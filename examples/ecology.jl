using Interior, Random, JSON3
# A finite nutrient inventory; no replenishment.
base=joinpath(@__DIR__,"..","outputs","ecology");mkpath(base)
out=isempty(ARGS) ? mktempdir(base;prefix="run-",cleanup=false) : abspath(ARGS[1])
ispath(joinpath(out,"state-relaxed-final")) && error("Choose a fresh output directory; this one already contains a state recording")
mkpath(out)
rng=MersenneTwister(47)
w=World(Grid(64,64;extent=(64,64,1));seed=47)
for y in 1:64,x in 1:64
    w.fields.nutrient[x,y,1]=Float32(0.03+1.5*exp(-((x-32)^2+(y-32)^2)/600))
end
p=Ecology(vcat(64rand(rng,Float32,2,240),zeros(Float32,1,240));seed=47,capacity=1200,lifespan=50,speed=0.8)
initial=ecology_mass(w,p)
stats=Any[]
for tick in 0:1200
    if tick%20==0
        push!(stats,merge((tick=tick,time=w.time,mass=ecology_mass(w,p),mass_error=ecology_mass(w,p)-initial),ecology_statistics(p)))
    end
    tick<1200 && step!(w,p,0.2f0)
end
open(joinpath(out,"statistics.json"),"w") do io; JSON3.pretty(io,stats); end
# Extrude deposit into a shallow volume; agent storage is separately recorded.
volume=repeat(w.fields.deposit,1,1,8)
open(joinpath(out,"deposit-64x64x8.f32"),"w") do io; write(io,volume); end
open(joinpath(out,"population.json"),"w") do io
    JSON3.pretty(io,(ids=p.ids,parents=p.parents,lineages=p.lineages,generations=p.generations,positions=[collect(p.positions[:,i]) for i in eachindex(p.ids)],deposit_shape=size(volume),initial_material=initial,statistics=ecology_statistics(p)))
end
record_state(w,(p,),joinpath(out,"state-relaxed-final");step=1200)
println(last(stats))

println("Saved ecology to ",out)
