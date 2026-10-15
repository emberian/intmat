# Run standalone: this test injects an arena_step! failure in its Julia process.
# It must not be included in the ordinary in-process physics test harness.
using Interior, Test, TOML, JSON3
include(joinpath(@__DIR__,"..","scripts","cellular_soak.jl"))
@eval Interior function arena_step!(a::Arena)
    a.world.fields.nutrient[1]=777f0
    error("injected mid-step failure after field mutation")
end
@testset "Soak failure boundaries and resume profile contract" begin
    mktempdir() do dir
        cfg=Dict("model"=>"cellular","seed"=>81,"dt"=>.1,"grid"=>Dict("shape"=>[8,8,1],"extent"=>[8.,8.,1.]),"cellular"=>Dict("population"=>2))
        profile=joinpath(dir,"profile.toml");open(io->TOML.print(io,cfg),profile,"w")
        output=joinpath(dir,"run")
        @test_throws ErrorException soak_main([profile,output,"--seconds","300","--budget-bytes","32000000"])
        status=JSON3.read(read(joinpath(output,"status.json"),String))
        @test status.state=="failed"
        @test status.unsafe_in_memory_state
        @test status.checkpoint_preserved
        @test status.last_completed_step==0
        pointer=JSON3.read(read(joinpath(output,"checkpoint.json"),String))
        checkpoint=joinpath(output,String(pointer.file))
        saved=restore_arena(checkpoint)
        @test saved.steps==0
        @test saved.world.fields.nutrient[1]==0f0
        @test length(filter(x->startswith(x,"checkpoint-")&&endswith(x,".jls"),readdir(output)))==1
        @test !ispath(joinpath(output,".arena-writer"))
        before=read(joinpath(output,"status.json"))
        cfg["seed"]=82;open(io->TOML.print(io,cfg),profile,"w")
        @test_throws ErrorException soak_main([profile,output,"--resume",checkpoint,"--seconds","1"])
        @test read(joinpath(output,"status.json"))==before
    end
end
