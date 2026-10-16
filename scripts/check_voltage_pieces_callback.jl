include(joinpath(@__DIR__,"../studies/voltage_between_pieces.jl"))
using Test
root=ARGS[1];ispath(root)&&error("fresh check directory required");mkpath(root)
produce(joinpath(root,"plain");dims=(5,4,3),steps=60)
calls=Ref(0)
function observe(s)
 calls[]+=1
 @test all(isfinite,s.potential);@test s.solver.residual<2e-6
 @test s.step>=0;@test s.time>=0
 Interior.mechanical_energy(s.process) # Same read-only diagnostic used by film writer.
end
produce(joinpath(root,"observed");dims=(5,4,3),steps=60,snapshot=observe)
a=read(joinpath(root,"plain","state.json"));b=read(joinpath(root,"observed","state.json"));@test a==b;@test calls[]==4
result=(snapshot_calls=calls[],read_only_callback_final_state_byte_exact=true,state_sha256=bytes2hex(sha256(a)),qualification="Small complete coupled trajectory at unchanged default cadence, callback disabled versus read-only observation. No artwork rerun.")
open(io->JSON3.pretty(io,result),joinpath(root,"result.json"),"w");println(JSON3.write(result))
