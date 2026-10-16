include(joinpath(@__DIR__,"../studies/wave_memory.jl"))
using Test
kind=get(ENV,"WAVE_BACKEND","cpu")
if kind=="metal";using Metal;backend=Metal.MetalBackend();else;backend=CPU();end
function verify(backend)
 a=packet_memory_fixture(CPU(),16);b=packet_memory_fixture(backend,16)
 for _=1:100;wave_memory_step!(a...,.006);wave_memory_step!(b...,.006);end
 af=wave_memory_fields(a...);bf=wave_memory_fields(b...)
 errors=Dict(string(k)=>maximum(abs.(getfield(af,k).-getfield(bf,k))) for k in keys(af));@test maximum(values(errors))<1e-4
 host=Array(b[3].host);@test minimum(bf.deposit)>=0;@test minimum(bf.nutrient)>=0
 inventory=maximum(abs.(bf.deposit.+bf.nutrient.-host));@test inventory<=1.2e-7
 # With no acoustic pressure, only independently specified dissolution remains.
 w,p,m=memory_fixture(CPU(),16;relaxation=.25);p.source_amplitude=0;initial=copy(w.fields.deposit)
 for _=1:100;wave_memory_step!(w,p,m,.006);end
 expected=initial.*exp(-100*.006f0*.25f0);decay=maximum(abs.(w.fields.deposit.-expected));@test decay<2e-7;@test all(iszero,p.pressure)
 (backend=kind,steps=100,errors=errors,local_inventory_error=inventory,zero_wave_decay_error=decay,nonnegative_donors=true,qualification="16³ portability and finite-donor fixture; not artwork-resolution convergence.")
end
result=verify(backend);open(io->JSON3.pretty(io,result),ARGS[1],"w");println(JSON3.write(result))
