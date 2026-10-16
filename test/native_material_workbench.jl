using Test, Interior, TOML, Serialization
include(joinpath(@__DIR__,"..","scripts/native_workbench_worker.jl"))
const MW=NativeWorkbench
@testset "Fork-only physical material controls" begin
    cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments/native/living-patch-recruitment-16.toml"))
    cfg["grid"]["shape"]=[8,8,8];cfg["cellular"]["founders"]=4;cfg["cellular"]["capacity"]=256
    mktempdir() do temp
        s=MW.create_session(cfg,joinpath(temp,"worker");side=4,living=true,limits=MW.Limits(wall_seconds=3600.,max_commands=50))
        cmd(id,action,branch;kwargs...)=MW.process_command!(s,merge(Dict{String,Any}("id"=>id,"action"=>action,"branch"=>branch),Dict(string(k)=>v for (k,v) in kwargs)))
        rates=Dict{String,Any}("decay_rate"=>2.,"desorption_rate"=>1.,"raw_diffusion"=>.03)
        original=Interior.device_habitat_snapshot(s.branches["sensing"].process)
        for branch in ("both","off","sensing")
            ack=cmd("original-"*branch,"material_rates",branch;rates)
            @test !ack.ok && !ack.accepted
        end
        @test Interior.device_habitat_snapshot(s.branches["sensing"].process)==original
        @test cmd("grow","step","sensing";steps=4).ok
        @test cmd("fork","fork","sensing").ok
        parent=s.branches["sensing"];child=s.branches["fork1"]
        parent_before=Interior.device_habitat_snapshot(parent.process);child_before=Interior.device_habitat_snapshot(child.process)
        policy_before=MW.living_snapshot(child.policy);parent_checkpoint=parent.checkpoint;parent_hash=MW.hashfile(joinpath(s.root,parent_checkpoint.file))
        for (i,bad) in enumerate((merge(rates,Dict("decay_rate"=>true)),merge(rates,Dict("raw_diffusion"=>-1.)),merge(rates,Dict("decay_rate"=>1e-100)),merge(rates,Dict("decay_rate"=>11.)),merge(rates,Dict("raw_diffusion"=>1.)),Dict("decay_rate"=>1.),merge(rates,Dict("capacity_per_enzyme"=>1.))))
            ack=cmd("invalid-"*string(i),"material_rates","fork1";rates=bad)
            @test !ack.ok && !ack.accepted
        end
        @test_throws ArgumentError MW.material_rates(merge(rates,Dict("desorption_rate"=>Inf)))
        @test_throws ArgumentError MW.material_rates(merge(rates,Dict("decay_rate"=>NaN)))
        @test !cmd("wrong-field","observe","fork1";rates).ok
        @test Interior.device_habitat_snapshot(child.process)==child_before
        apply=Dict{String,Any}("id"=>"apply","action"=>"material_rates","branch"=>"fork1","rates"=>rates)
        ack=MW.process_command!(s,apply);@test ack.ok
        changed=Interior.device_habitat_snapshot(child.process)
        @test changed.cells==child_before.cells && changed.arrays==child_before.arrays && changed.feedstock==child_before.feedstock
        @test MW.living_snapshot(child.policy)==policy_before
        @test child.policy.habitat===child.process.habitat
        @test typeof(child.process)===typeof(parent.process)
        @test changed.parameters.decay_rate==2f0 && changed.parameters.desorption_rate==1f0
        @test changed.feed_parameters.diffusion==(.03f0,.03f0,0f0)
        @test Interior.device_habitat_snapshot(parent.process)==parent_before
        @test MW.hashfile(joinpath(s.root,parent_checkpoint.file))==parent_hash
        events=child.history[end-1:end]
        metas=[Interior.JSON3.read(read(joinpath(s.root,"branches/fork1",e.file),String)) for e in events]
        @test metas[1].time==metas[2].time
        @test metas[1].binary_sha256==metas[2].binary_sha256
        @test metas[2].feed_parameters.diffusion[1]≈.03
        @test metas[2].intervention.action=="material_rates"
        @test MW.process_command!(s,apply)==ack
        @test length(child.interventions)==1
        payload=deserialize(joinpath(s.root,child.checkpoint.file))
        restored=Interior.device_habitat_restore(World(child.world.grid,CPU(),NamedTuple(),NamedTuple(),child.world.time,child.world.seed),payload)
        @test Interior.device_habitat_snapshot(restored)==changed
        @test payload.living==policy_before
        @test cmd("fork-edited","fork","fork1").ok
        @test Interior.device_habitat_snapshot(s.branches["fork2"].process)==changed
        @test s.branches["fork2"].interventions==child.interventions
        @test cmd("child-future","step","fork1";steps=8).ok
        @test Interior.device_habitat_snapshot(parent.process)==parent_before
        @test cmd("parent-future","step","sensing";steps=8).ok
        @test child.world.time==parent.world.time
        @test child.observation.material_summary.catalyst!=parent.observation.material_summary.catalyst
        @test child.observation.material_summary.bound_food!=parent.observation.material_summary.bound_food
        @test MW.accounting(child).audit.finite && MW.accounting(parent).audit.finite
        @test s.dt==Float32(cfg["native"]["dt"])
        @test cmd("zero-rates","material_rates","fork2";rates=Dict(k=>0. for k in keys(rates))).ok
        @test s.branches["fork2"].process.habitat.feedstock.parameters.diffusion==(0f0,0f0,0f0)
        @test cmd("stop","stop","both").ok
    end
end
