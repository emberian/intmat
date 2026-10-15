# Invoke with --project pointing at the checkpoint's frozen source/environment.
using Interior, Serialization, JSON3
include(joinpath(dirname(pathof(Interior)),"..","studies","constructive_world_batches.jl"))
function resume_constructive_world_batch(checkpoint_path,out,steps;dt=.01)
    ispath(out)&&throw(ArgumentError("fresh output directory required"))
    steps>=0||throw(ArgumentError("nonnegative additional steps required"))
    payload=open(deserialize,checkpoint_path)
    worlds=map(payload.members) do member
        g=member.cells.snapshot.geometry
        grid=Grid(Tuple(Int.(g.dims))...;extent=g.extent)
        World(grid,CPU(),NamedTuple(),NamedTuple(),0.,member.cells.world_seed)
    end
    batch=constructive_world_batch_restore(worlds,payload)
    initial=[Interior.device_habitat_audit(constructive_world_batch_snapshot(batch,i)) for i in eachindex(worlds)]
    for _ in 1:steps;constructive_world_batch_step!(worlds,batch,dt);end
    final=[Interior.device_habitat_audit(constructive_world_batch_snapshot(batch,i)) for i in eachindex(worlds)]
    material=[b.material-a.material-b.supplied+a.supplied+b.removed-a.removed for (a,b) in zip(initial,final)]
    potential=[b.potential-a.potential+b.work-a.work+b.heat-a.heat-b.supplied_potential+a.supplied_potential+b.removed_potential-a.removed_potential for (a,b) in zip(initial,final)]
    result=(source=CONSTRUCTIVE_BATCH_SOURCE,member_indices=payload.member_indices,
        additional_steps=steps,dt=dt,times=[w.time for w in worlds],initial=initial,final=final,
        material_errors=material,potential_work_heat_errors=potential)
    mkdir(out)
    open(joinpath(out,"result.json"),"w") do io;JSON3.pretty(io,result);println(io);end
    open(joinpath(out,"final.jls"),"w") do io;serialize(io,constructive_world_batch_checkpoint(worlds,batch));end
    result
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==3||error("usage: CHECKPOINT.jls FRESH_OUT_DIR ADDITIONAL_STEPS")
    result=resume_constructive_world_batch(ARGS[1],ARGS[2],parse(Int,ARGS[3]))
    println((members=length(result.final),steps=result.additional_steps,times=result.times,
        material=maximum(abs,result.material_errors),potential=maximum(abs,result.potential_work_heat_errors)))
end
