include(joinpath(@__DIR__,"device_cellular_setup.jl"))
include(joinpath(@__DIR__,"..","studies","cohort_construction_world_batches.jl"))
include(joinpath(@__DIR__,"..","studies","native_cohort_construction_credit.jl"))
const Q=NativeCohortConstructionCredit
register_cohort_construction_world_batches!(Q)
f=Q.C.fixture(joinpath(@__DIR__,"..","experiments/native/constructed-garden-workbench-16.toml");reserve_fraction=.2)
cw=[World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed) for _ in 1:4]
cp=[Interior.device_habitat_restore(w,f.checkpoint) for w in cw]
ch=[Q.sidecar(p;dt=.1,horizon=2.,act_time=4.,seed=6201,memory=isodd(i),learning=i<=2,critic=i>2,rate=i<=2 ? (i==1 ? .003 : -.003) : 0.) for (i,p) in enumerate(cp)]
for i in 1:4;ch[i].options=merge(ch[i].options,(critic_gain=isodd(i) ? .01f0 : .03f0,));end
b=learning_world_batch(cp,ch;backend,identities=["actor-positive","actor-reversed","critic-gain01","critic-gain03"])
bw=[World(f.grid,backend,NamedTuple(),NamedTuple(),0.,f.seed) for _ in 1:4]
origin=[Interior.device_habitat_audit(x.physical) for x in learning_world_batch_snapshot(b)]
errors=[];checks=0;restored=nothing;rw=nothing
for n in (1,18,20,21)
 learning_serial_advance!(cw,cp,ch,n,.1)
 learning_world_batch_advance!(bw,b,n,.1)
 x=learning_world_batch_snapshot(b)
 push!(errors,learning_batch_parity(learning_serial_snapshot(cp,ch),x;tolerance=2e-5));global checks+=1
 if restored!==nothing
  learning_world_batch_advance!(rw,restored,n,.1)
  learning_world_batch_snapshot(restored)==x||error("continued native partial-queue restart differs")
  global checks+=1
 end
 if b.learner.steps==19
  payload=learning_world_batch_checkpoint(bw,b)
  global rw=[World(f.grid,backend,NamedTuple(),NamedTuple(),0.,f.seed) for _ in 1:4]
  global restored=learning_world_batch_restore(rw,payload)
  learning_world_batch_snapshot(restored)==x||error("immediate native partial-queue restart differs")
  global checks+=1
 end
 println("cohort batch steps ",b.learner.steps," error ",last(errors));flush(stdout)
end
rows=map(zip(origin,learning_world_batch_snapshot(b))) do (a,x)
 z=Interior.device_habitat_audit(x.physical)
 material=z.material-a.material-z.supplied+a.supplied+z.removed-a.removed
 potential=z.potential-a.potential+z.work-a.work+z.heat-a.heat-z.supplied_potential+a.supplied_potential+z.removed_potential-a.removed_potential
 paid=sum(x.state.paid);matured=sum(x.state.matured)
 abs(material)<1e-8&&abs(potential)<1e-8&&paid>0&&matured>0||error("native cohort paid accounting/maturity failed")
 (;material,potential,paid,matured,critic_gain=x.options.critic_gain,memory=x.properties.memory,actor_updates=sum(x.state.actor_updates),critic_updates=sum(x.state.critic_updates))
end
checks+=4
open(io->JSON3.pretty(io,(;checks,errors,rows,time=getproperty.(bw,:time),source=LEARNING_BATCH_SOURCE,policy=b.contract.policy.source)),ARGS[1],"w")
println("native cohort batch gate passed ",checks," checks");flush(stdout)
