t0=time_ns()
include(joinpath(@__DIR__,"..","scripts/device_cellular_setup.jl"))
include(joinpath(@__DIR__,"..","studies/learning_world_batches.jl"))
include(joinpath(@__DIR__,"..","studies/native_joint_construction_credit.jl"))
const J=NativeJointConstructionCredit
register_joint_construction_world_batches!(J)
profile=joinpath(@__DIR__,"..","experiments/native/constructed-garden-workbench-16.toml")
f=J.C.fixture(profile;reserve_fraction=.2)
loading_seconds=(time_ns()-t0)/1e9
function cohort(n)
 w=[World(f.grid,backend,NamedTuple(),NamedTuple(),0.,f.seed) for _ in 1:n]
 p=[Interior.device_habitat_restore(x,f.checkpoint) for x in w]
 h=[J.sidecar(x;dt=.02,seed=isodd(i) ? 6201 : 6277,memory=isodd(i)) for (i,x) in enumerate(p)]
 w,p,h
end
sw,sp,sh=cohort(8)
bw,bp,bh=cohort(8)
b=learning_world_batch(bp,bh;identities=["joint-$i" for i in 1:8])
println("joint batch owned bytes ",learning_world_batch_memory(b));flush(stdout)
origin=[Interior.device_habitat_audit(x.physical) for x in learning_world_batch_snapshot(b)]
errors=Float64[];timings=[]
# Identical observation schedule and asynchronous serial baseline, one barrier per block.
for (trial,steps) in enumerate((1,99,100,100))
 s=@elapsed begin
  learning_serial_advance!(sw,sp,sh,steps,.02)
  global serial=learning_serial_snapshot(sp,sh)
  map(learning_snapshot_observation,serial)
 end
 q=@elapsed begin
  learning_world_batch_advance!(bw,b,steps,.02)
  global batch=learning_world_batch_snapshot(b)
  map(learning_snapshot_observation,batch)
 end
 push!(errors,learning_batch_parity(serial,batch;tolerance=2e-5))
 push!(timings,(;trial,steps,serial_seconds=s,batch_seconds=q,speedup=s/q))
 println("joint trial ",trial," steps ",steps," serial ",s," batch ",q," error ",last(errors));flush(stdout)
end
cp=learning_world_batch_checkpoint(bw,b)
rw=[World(f.grid,backend,NamedTuple(),NamedTuple(),0.,f.seed) for _ in 1:8]
r=learning_world_batch_restore(rw,cp)
learning_world_batch_snapshot(r)==learning_world_batch_snapshot(b)||error("strict restart differs")
residuals=map(zip(origin,learning_world_batch_snapshot(b))) do (a,x)
 z=Interior.device_habitat_audit(x.physical)
 material=z.material-a.material-z.supplied+a.supplied+z.removed-a.removed
 potential=z.potential-a.potential+z.work-a.work+z.heat-a.heat-z.supplied_potential+a.supplied_potential+z.removed_potential-a.removed_potential
 abs(material)<1e-8&&abs(potential)<1e-8&&sum(x.state.paid)>0||error("paid accounting failure")
 (;material,potential,paid=sum(x.state.paid))
end
open(io->JSON3.pretty(io,(;loading_seconds,timings,errors,residuals,memory_bytes=learning_world_batch_memory(b),source=LEARNING_BATCH_SOURCE,policy_source=b.contract.policy.source,total_steps=300,dt=.02)),ARGS[1],"w")
println("joint native batch gate complete");flush(stdout)
