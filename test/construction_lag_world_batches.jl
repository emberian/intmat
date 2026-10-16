include(joinpath(@__DIR__,"..","studies","construction_lag_world_batches.jl"))

using Interior, Test, KernelAbstractions
include(joinpath(@__DIR__,"..","studies","native_construction_lag.jl"))
const P=NativeConstructionLag
register_construction_lag_world_batches!(P)
f=P.C.fixture(joinpath(@__DIR__,"..","experiments","native","constructed-garden-workbench-16.toml");reserve_fraction=.2)
w=[World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed) for _ in 1:2]
p=[Interior.device_habitat_restore(x,f.checkpoint) for x in w]
h=[P.pulse(P.J.sidecar(x;dt=.1,seed=6201,memory=i==1);axis=i,amplitude=isodd(i) ? .15 : -.15,producer=i,duration=.2) for (i,x) in enumerate(p)]
b=learning_world_batch(p,h)
bw=[World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed) for _ in 1:2]
@testset "Pulse policy private commands and strict partial restart" begin
 for n in (1,1,18)
  learning_serial_advance!(w,p,h,n,.1)
  learning_world_batch_advance!(bw,b,n,.1)
  a=learning_world_batch_snapshot(b)
  @test learning_batch_parity(learning_serial_snapshot(p,h),a)==0
  for i in 1:2
   @test a[i].cache.saved_drive==h[i].saved_drive
   @test a[i].properties.axis==i
   @test a[i].properties.producer==i
  end
  cp=learning_world_batch_checkpoint(bw,b)
  rw=[World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed) for _ in 1:2]
  restored=learning_world_batch_restore(rw,cp)
  @test learning_world_batch_snapshot(restored)==a
 end
 other=P.PulseHabitat(h[2].base,h[2].axis,h[2].amplitude,h[2].producer,h[2].saved_drive,1,h[2].duration_steps)
 @test_throws ArgumentError learning_world_batch(p,[h[1],other])
end
