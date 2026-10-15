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

@testset "Explicit source migration forks and throughput provenance" begin
    mktempdir() do dir
        cfg=Dict("model"=>"cellular","seed"=>83,"dt"=>.1,"grid"=>Dict("shape"=>[8,8,1],"extent"=>[8.,8.,1.]),"cellular"=>Dict("population"=>2),"sleep_seconds"=>1.0)
        profile=joinpath(dir,"profile.toml");open(io->TOML.print(io,cfg),profile,"w")
        parent=joinpath(dir,"parent.jls");checkpoint!(build_arena(cfg),parent)
        e=open(Serialization.deserialize,parent);old=copy(e.sources);old["cellular.jl"]="previous implementation"
        open(io->Serialization.serialize(io,merge(e,(sources=old,))),parent,"w")
        strict=joinpath(dir,"strict")
        @test_throws ArgumentError soak_main([profile,strict,"--resume",parent,"--seconds","1"])
        @test !ispath(strict)
        @test_throws ErrorException soak_main([profile,strict,"--allow-source-change"])
        existing=joinpath(dir,"existing");mkdir(existing);write(joinpath(existing,"keep"),"user evidence")
        @test_throws ErrorException soak_main([profile,existing,"--resume",parent,"--allow-source-change"])
        @test read(joinpath(existing,"keep"),String)=="user evidence"
        output=joinpath(dir,"fork")
        @test_throws ErrorException soak_main([profile,output,"--resume",parent,"--allow-source-change","--throughput","--seconds","300","--budget-bytes","32000000"])
        @test read(joinpath(output,"parent-checkpoint.jls"))==read(parent)
        sessionfile=only(filter(x->startswith(x,"session-"),readdir(output)))
        session=JSON3.read(read(joinpath(output,sessionfile),String))
        @test session.throughput
        @test session.batch_sleep_seconds==0
        @test session.source_migration.initial_step==0
        @test session.source_migration.parent_checkpoint_sha256==bytes2hex(sha256(read(parent)))
        @test session.source_migration.parent_sources["cellular.jl"]=="previous implementation"
        @test session.source_migration.loaded_sources["cellular.jl"]==Interior.arena_sources()["cellular.jl"]
        pointer=JSON3.read(read(joinpath(output,"checkpoint.json"),String))
        @test restore_arena(joinpath(output,String(pointer.file))).world.fields.nutrient[1]==0f0
    end
end
