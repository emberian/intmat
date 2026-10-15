# Runner validation does not require a GPU package.
include(joinpath(@__DIR__,"..","scripts","run_recurrent_material.jl"))
@testset "bounded portable recurrent runner" begin
    o=parse_recurrent_args(["8","10000000","unused","--backend","cpu","--profile","all"])
    @test length(o.selected)==3
    @test length(o.frames)==12 && length(o.checkpoints)==1
    @test maximum(o.frames)==o.steps
    @test o.capture_mode=="slices"
    @test o.estimated_bytes<=o.budget_bytes
    @test o.resident_array_bytes==8^3*84
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
        @test read(joinpath(root,"loaded-core-source.jl"),String)==loaded_core_source
        @test count(p->endswith(p,".f32"),readdir(joinpath(root,"braided","state-9")))==4
        @test_throws ArgumentError run_recurrent_material(o)
        allroot=joinpath(dir,"all-profiles")
        alloptions=parse_recurrent_args(["3","3",allroot,"--backend","cpu","--profile","all","--captures","1","--states","0"])
        run_recurrent_material(alloptions)
        @test all(name->isfile(joinpath(allroot,name,"metrics.json")),first.(recurrent_profiles))
        @test count(p->endswith(p,".f32"),readdir(joinpath(allroot,"slow-memory")))==6
    end
end

@testset "compact world preserves recurrent trajectory" begin
    grid=Grid(5,6,4;extent=(5,6,4))
    full=World(grid;seed=31)
    compact=recurrent_world(grid,Interior.CPU())
    @test propertynames(compact.fields)==(:nutrient,:deposit,:activity)
    @test propertynames(compact.scratch)==(:nutrient,:deposit)
    buffers=(values(compact.fields)...,values(compact.scratch)...)
    @test all(!Base.mightalias(buffers[i],buffers[j]) for i in eachindex(buffers) for j in 1:i-1)
    p=RecurrentMaterial(full);q=RecurrentMaterial(compact)
    old_initial=clamp.(.12f0.*randn(MersenneTwister(31),Float32,grid.size...,8),-1f0,1f0)
    @test p.state==q.state==old_initial
    for _ in 1:101
        step!(full,p,.05);step!(compact,q,.05)
    end
    @test full.time==compact.time
    @test p.state==q.state && p.scratch==q.scratch
    @test all(getproperty(full.fields,n)==getproperty(compact.fields,n) for n in propertynames(compact.fields))
    @test sum(sizeof,buffers)+sizeof(q.state)+sizeof(q.scratch)==prod(grid.size)*84
    @test_throws ArgumentError World(grid;fields=(:nutrient,:nutrient))
    @test_throws ArgumentError World(grid;fields=(:unknown,))
    @test_throws ArgumentError World(grid;fields=(:activity,),scratch_fields=(:nutrient,))
    @test_throws ArgumentError World(grid;fields=(:activity,),scratch_fields=(:activity,:activity))
end
