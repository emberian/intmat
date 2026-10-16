using Interior, Test, KernelAbstractions
include(joinpath(@__DIR__,"..","studies","cohort_construction_world_batches.jl"))
include(joinpath(@__DIR__,"..","studies","native_cohort_construction_credit.jl"))
const Q=NativeCohortConstructionCredit
register_cohort_construction_world_batches!(Q)
f=Q.C.fixture(joinpath(@__DIR__,"..","experiments/native/constructed-garden-workbench-16.toml");reserve_fraction=.2)
w=[World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed) for _ in 1:2]
p=[Interior.device_habitat_restore(x,f.checkpoint) for x in w]
h=[Q.sidecar(x;dt=.1,horizon=2.,act_time=4.,seed=6201,memory=i==1,learning=i==1,critic=i==2,rate=i==1 ? .003 : 0.) for (i,x) in enumerate(p)]
b=learning_world_batch(p,h)
bw=[World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed) for _ in 1:2]
@testset "Decision-aligned cohort queues share original hooks" begin
 for n in (1,38,21)
  learning_serial_advance!(w,p,h,n,.1)
  learning_world_batch_advance!(bw,b,n,.1)
  a=learning_world_batch_snapshot(b)
  @test learning_batch_parity(learning_serial_snapshot(p,h),a)==0
  for i in 1:2
   @test a[i].state.queue_ready==h[i].state.queue_ready
   @test a[i].state.queue_value==h[i].state.queue_value
   @test a[i].cache.original_masks==h[i].original_masks
  end
  cp=learning_world_batch_checkpoint(bw,b)
  rw=[World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed) for _ in 1:2]
  restored=learning_world_batch_restore(rw,cp)
  @test learning_world_batch_snapshot(restored)==a
 end
 a=learning_world_batch_snapshot(b)
 @test sum(a[1].state.matured)>0
 @test sum(a[2].state.matured)>0
 @test sum(a[1].state.paid)>0
 @test sum(a[2].state.paid)>0
end
