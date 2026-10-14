@testset "Arena budgets, checkpoints and control" begin
    w=World(Grid(8,8;extent=(8,8,1)))
    a=Arena(w,(Chemostat(source_rate=.03,sink_rate=.02),);dt=.1)
    for _ in 1:10; arena_step!(a); end
    snap=arena_snapshot(a;max_side=4,max_agents=0)
    @test abs(snap["metrics"]["material_error"])<1e-10
    @test snap["fields"]["nutrient"]["shape"]==[4,4]
    mktempdir() do dir
        checkpoint!(a,joinpath(dir,"saved.jls")); b=restore_arena(joinpath(dir,"saved.jls"))
        for _ in 1:10; arena_step!(a); arena_step!(b); end
        @test a.world.fields.nutrient==b.world.fields.nutrient
        @test a.world.time==b.world.time
        @test a.processes[1].supplied==b.processes[1].supplied
        run_arena!(b,joinpath(dir,"output");steps=2,checkpoint_seconds=100)
        @test b.steps==22
        @test isfile(joinpath(dir,"output","checkpoint.json"))
    end
    Interior.arena_command!(a,Dict("id"=>"pause-1","action"=>"pause")); n=a.steps; arena_step!(a); @test a.steps==n
    Interior.arena_command!(a,Dict("id"=>"resume-1","action"=>"resume")); arena_step!(a); @test a.steps==n+1
    @test_throws ArgumentError Interior.arena_command!(a,Dict("id"=>"invalid","action"=>"set","values"=>Dict("dt"=>-1)))
end
@testset "Coupled arena deterministic restart" begin
    cfg=Dict("grid"=>Dict("shape"=>[8,8,1],"extent"=>[8,8,1]),"ecology"=>Dict("population"=>8,"capacity"=>20))
    a=build_arena(cfg); b=build_arena(cfg)
    for _ in 1:5; arena_step!(a); arena_step!(b); end
    @test a.world.fields==b.world.fields
    @test a.processes[2].weights==b.processes[2].weights
    mktempdir() do dir
        path=joinpath(dir,"restart.jls"); checkpoint!(a,path); c=restore_arena(path)
        @test c.processes[3].phase===c.processes[4]
        @test c.processes[2].medium===c.processes[3]
        @test c.processes[3].membrane===c.world.fields.deposit
        for _ in 1:5; arena_step!(a); arena_step!(c); end
        @test a.world.fields==c.world.fields
        @test a.processes[2].weights==c.processes[2].weights
        @test a.processes[2].plastic_marks==c.processes[2].plastic_marks
        @test a.processes[4].phi==c.processes[4].phi
    end
end
@testset "Arena output ownership and writer lock" begin
    mktempdir() do dir
        unowned=joinpath(dir,"unowned");mkpath(unowned)
        write(joinpath(unowned,"latest.json"),"user content")
        a=Arena(World(Grid(4,4)),(Chemostat(),);dt=.1)
        @test_throws ArgumentError run_arena!(a,unowned;steps=0)
        @test read(joinpath(unowned,"latest.json"),String)=="user content"
        owned=joinpath(dir,"owned")
        lock=Interior.arena_lock(owned)
        @test_throws ArgumentError Interior.arena_lock(owned;resume=true)
        Interior.arena_unlock(lock)
        @test_throws ArgumentError run_arena!(a,owned;steps=0)
        run_arena!(a,owned;steps=1,resume=true)
        c=restore_arena(joinpath(owned,"checkpoint-1.jls"))
        run_arena!(c,owned;steps=1,resume=true)
        @test c.steps==2
        @test !ispath(joinpath(owned,".arena-writer"))
        stale=joinpath(owned,".arena-writer");mkdir(stale)
        write(joinpath(stale,"pid.json"),"{\"schema\":\"interior-arena-lock\",\"pid\":2147483647}")
        run_arena!(c,owned;steps=0,resume=true)
        @test !ispath(stale)
    end
end

@testset "Loaded implementation provenance is isolated" begin
    original=Interior.arena_sources()
    edited=Interior.arena_sources()
    edited["arena.jl"]="modified caller copy"
    @test Interior.arena_sources()==original
    @test edited!=original
end
