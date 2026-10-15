# Runner validation does not require a GPU package.
include(joinpath(@__DIR__,"..","scripts","run_recurrent_material.jl"))
@testset "bounded portable recurrent runner" begin
    o=parse_recurrent_args(["8","10000000","unused","--backend","cpu","--profile","all"])
    @test length(o.selected)==3
    @test length(o.frames)==12 && length(o.checkpoints)==1
    @test maximum(o.frames)==o.steps
    @test o.capture_mode=="slices"
    @test o.estimated_bytes<=o.budget_bytes
    @test_throws ArgumentError parse_recurrent_args(["8","10","unused","--backend","unknown"])
    @test_throws ArgumentError parse_recurrent_args(["8","10","unused","--captures","-1"])
    @test_throws ArgumentError parse_recurrent_args(["256","1000","unused","--capture-mode","volumes","--budget-mb","1"])
    @test_throws ArgumentError parse_recurrent_args(["512","1000","unused","--resident-mb","1024"])
    @test recurrent_capture_steps(3,20)==Set(1:3)
    @test isempty(recurrent_capture_steps(3,0))
    mktempdir() do dir
        root=joinpath(dir,"cpu")
        o=parse_recurrent_args(["3","9",root,"--backend","cpu","--profile","braided","--captures","3","--states","2"])
        run_recurrent_material(o)
        @test count(p->endswith(p,".f32"),readdir(joinpath(root,"braided")))==18
        @test isdir(joinpath(root,"braided","state-4")) && isdir(joinpath(root,"braided","state-9"))
        @test isfile(joinpath(root,"run-config.json"))
        @test_throws ArgumentError run_recurrent_material(o)
        allroot=joinpath(dir,"all-profiles")
        alloptions=parse_recurrent_args(["3","3",allroot,"--backend","cpu","--profile","all","--captures","1","--states","0"])
        run_recurrent_material(alloptions)
        @test all(name->isfile(joinpath(allroot,name,"metrics.json")),first.(recurrent_profiles))
        @test count(p->endswith(p,".f32"),readdir(joinpath(allroot,"slow-memory")))==6
    end
end
