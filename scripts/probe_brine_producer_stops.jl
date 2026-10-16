# Failure-path test: diagnostics must survive stops and a failed full export.
include(joinpath(@__DIR__,"..","studies","native_brine_producer.jl"))
using Test
root=ARGS[1];mkpath(root)
c=TOML.parsefile(joinpath(@__DIR__,"..","experiments/art/brine-porous-384-f64.toml"))
c["resolution"]=16;c["steps"]=1;c["replay_resolution"]=8;c["history_stride"]=1;c["replay_stride"]=1
for mode in (:complete,:wall,:export_failure)
 c["wall_seconds"]=mode===:wall ? 0. : 1000.
 recipe=joinpath(root,string(mode)*".toml");open(io->TOML.print(io,c),recipe,"w")
 capture=function(x,out,step,m,frames;full=false)
  full&&mode===:export_failure&&error("injected diagnostic export failure")
  producer_capture(x,out,step,m,frames;full)
 end
 out=joinpath(root,string(mode));report=native_brine_producer(recipe,out;capture)
 @test report.success==(mode===:complete)
 @test !report.resumable
 @test isfile(joinpath(out,"status.json"))
 if mode===:export_failure
  @test occursin("injected",report.export_failure)
 else
  frame=last(report.frames);@test frame.kind=="final"
  @test filesize(joinpath(out,frame.directory,"a.f64"))==16^3*8
  @test filesize(joinpath(out,frame.directory,"phi.f32"))==16^3*4
 end
end
println("producer stop/export tests passed")
