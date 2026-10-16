using Test
# Evaluate only the actual CLI error-handler block with a failing mock producer;
# this tests output ownership without importing or launching a GPU backend.
source=Meta.parseall(read(joinpath(@__DIR__,"..","studies","native_brine_producer.jl"),String))
entry=last(filter(ex->ex isa Expr&&ex.head===:if,source.args)).args[2]
@testset "Producer failure metadata never modifies historical outputs" begin
 mktempdir() do root
  old=joinpath(root,"old");mkpath(old);write(joinpath(old,"marker"),"preserved")
  for (path,fresh) in ((old,false),(joinpath(root,"fresh"),true))
   m=Module(gensym(:ProducerEntryMock))
   Core.eval(m,:(using Interior: JSON3))
   Core.eval(m,:(const ARGS=["recipe",$path]))
   Core.eval(m,:(function native_brine_producer(recipe,out;preflight_only=false)
    ispath(out)&&error("fresh output directory required")
    mkpath(out);error("injected setup failure")
   end))
   @test_throws ErrorException Core.eval(m,entry)
   @test isfile(joinpath(path,"failure.json"))==fresh
  end
  @test readdir(old)==["marker"]
  @test read(joinpath(old,"marker"),String)=="preserved"
 end
end
