using Test, Interior, TOML, Serialization
include(joinpath(@__DIR__,"..","scripts","native_workbench_worker.jl"))
using .NativeWorkbench
@testset "Bounded native workbench commands and actual observations" begin
    cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","native","constructed-garden-material-sensing-64.toml"))
    cfg["grid"]["shape"]=[8,8,8];cfg["grid"]["extent"]=[1.,1.,1.]
    cfg["cellular"]["founders"]=4;cfg["cellular"]["capacity"]=256;cfg["native"]["dt"]=.01
    for pool in ("resources","feedstock"),key in ("source_rates","source_budgets","sink_rates","sink_budgets")
        cfg[pool][key]=cfg[pool][key]./64
    end
    mktempdir() do temp
        root=joinpath(temp,"worker")
        token=NativeWorkbenchStartup.begin_startup(root,joinpath(@__DIR__,"..","scripts","native_workbench_worker.jl"))
        initial=Interior.JSON3.read(read(joinpath(root,"worker.json"),String))
        @test initial.status=="starting"
        @test_throws ErrorException NativeWorkbenchStartup.begin_startup(root,@__FILE__)
        @test_throws ErrorException NativeWorkbench.create_session(cfg,root;side=4,startup_token="wrong")
        s=NativeWorkbench.create_session(cfg,root;side=4,startup_token=token,
            limits=NativeWorkbench.Limits(max_steps=10,max_steps_per_command=4,max_model_time=1.,wall_seconds=3600.,idle_seconds=60.,max_commands=40))
        @test all(b->b.world.time==0,values(s.branches))
        metadata=Interior.JSON3.read(read(joinpath(s.root,"worker.json"),String))
        @test isconst(NativeWorkbench,:LOADED_CODE_PINS)
        @test metadata.driver_sha256==NativeWorkbench.LOADED_CODE_PINS.driver
        @test metadata.factory_sha256==NativeWorkbench.LOADED_CODE_PINS.factory
        @test metadata.capture_helper_sha256==NativeWorkbench.LOADED_CODE_PINS.capture
        @test metadata.startup_helper_sha256==NativeWorkbench.LOADED_CODE_PINS.startup
        @test s.started==initial.started_unix
        @test_throws ErrorException NativeWorkbench.create_session(cfg,root;side=4,startup_token=token)
        a=NativeWorkbench.process_command!(s,Dict("id"=>"one","action"=>"step","branch"=>"both","steps"=>3))
        @test a.ok && a.accepted && a.status=="completed"
        @test all(b->b.process.cells.attempted_steps==3,values(s.branches))
        duplicate=NativeWorkbench.process_command!(s,Dict("id"=>"one","action"=>"step","branch"=>"both","steps"=>3))
        @test duplicate==a && s.sequence==1
        altered=NativeWorkbench.process_command!(s,Dict("id"=>"one","action"=>"step","steps"=>4))
        @test !altered.ok && !altered.accepted
        for (id,extra) in (("badcount",Dict("steps"=>5)),("boolcount",Dict("steps"=>true)),("badbranch",Dict("steps"=>1,"branch"=>"../off")),("eval",Dict("steps"=>1,"code"=>"1+1")))
            ack=NativeWorkbench.process_command!(s,merge(Dict{String,Any}("id"=>id,"action"=>"step"),extra))
            @test !ack.ok && !ack.accepted
        end
        @test all(b->b.process.cells.attempted_steps==3,values(s.branches))
        paused=NativeWorkbench.process_command!(s,Dict("id"=>"pause","action"=>"pause"))
        @test paused.ok && all(b->b.world.time≈3Float64(s.dt),values(s.branches))
        checkpoint=NativeWorkbench.process_command!(s,Dict("id"=>"save","action"=>"checkpoint","branch"=>"sensing"))
        @test checkpoint.ok
        branch=s.branches["sensing"];payload=deserialize(joinpath(s.root,branch.checkpoint.file))
        w=World(branch.world.grid,CPU(),NamedTuple(),NamedTuple(),0.,cfg["seed"])
        restored=Interior.device_habitat_restore(w,payload)
        @test Interior.device_habitat_snapshot(restored)==Interior.device_habitat_snapshot(branch.process)
        meta=Interior.JSON3.read(read(joinpath(s.root,branch.observation.file),String))
        @test meta.step==3 && meta.shape==[4,4,4] && meta.native_shape==[8,8,8]
        @test meta.arrays.fields.bytes==4^3*9*4
        @test meta.arrays.hidden.shape==[4,8] && meta.arrays.material_inputs.shape==[4,3]
        binary=joinpath(dirname(joinpath(s.root,branch.observation.file)),meta.binary)
        @test NativeWorkbench.hashfile(binary)==meta.binary_sha256
        index=Interior.JSON3.read(read(joinpath(s.root,"branches","sensing","observatory.json"),String))
        @test length(index.frames)==2 && index.frames[end].step==3
        # A paused intervention changes the concrete receptor parameters only.
        branch=s.branches["sensing"];oldprocess=branch.process
        oldsnapshot=Interior.device_habitat_snapshot(oldprocess)
        for (name,values) in (("zero",[0.,1.,1.]),("underflow",[1e-100,1.,1.]),("overflow",[1e100,1.,1.]),("bool",Any[true,1.,1.]),("short",[1.,1.]),("nan",[NaN,1.,1.]))
            @test_throws ArgumentError NativeWorkbench.validate_command(Dict("id"=>name,"action"=>"receptor_scales","branch"=>"sensing","scales"=>values),s)
        end
        @test_throws ArgumentError NativeWorkbench.validate_command(Dict("id"=>"offscale","action"=>"receptor_scales","branch"=>"off","scales"=>[1.,1.,1.]),s)
        intervention=Dict("id"=>"receptors","action"=>"receptor_scales","branch"=>"sensing","scales"=>[.5,.03,.006])
        changed=NativeWorkbench.process_command!(s,intervention)
        @test changed.ok && changed.accepted
        @test branch.process!==oldprocess && branch.process.cells===oldprocess.cells
        @test branch.process.habitat.arrays===oldprocess.habitat.arrays
        @test branch.process.habitat.feedstock===oldprocess.habitat.feedstock
        now=Interior.device_habitat_snapshot(branch.process)
        @test now.cells==oldsnapshot.cells && now.arrays==oldsnapshot.arrays && now.feedstock==oldsnapshot.feedstock
        @test now.parameters.material_scales==(.5f0,.03f0,.006f0)
        @test now.parameters.material_channels==oldsnapshot.parameters.material_channels
        @test branch.world.time==3Float64(s.dt) && branch.process.cells.attempted_steps==3
        @test length(branch.interventions)==1 && branch.interventions[1].step==3
        @test NativeWorkbench.process_command!(s,intervention)==changed && length(branch.interventions)==1
        post=deserialize(joinpath(s.root,branch.checkpoint.file))
        @test post.parameters.material_scales==now.parameters.material_scales && post.workbench.intervention.phase=="after"
        @test post.workbench.interventions[1].id=="receptors"
        @test Interior.device_habitat_snapshot(Interior.device_habitat_restore(w,post))==now
        events=Interior.JSON3.read(read(joinpath(s.root,"branches","sensing","observatory.json"),String)).frames[end-1:end]
        @test events[1].time==events[2].time && events[1].source_sha256!=events[2].source_sha256
        beforemeta=Interior.JSON3.read(read(joinpath(s.root,"branches","sensing",events[1].file),String))
        aftermeta=Interior.JSON3.read(read(joinpath(s.root,"branches","sensing",events[2].file),String))
        @test beforemeta.binary_sha256==aftermeta.binary_sha256
        @test beforemeta.intervention.phase=="before" && aftermeta.intervention.phase=="after"
        @test aftermeta.material_parameters.material_scales≈[.5,.03,.006]
        @test_throws ArgumentError NativeWorkbench.reserve!(s,s.limits.artifact_bytes)
        more=NativeWorkbench.process_command!(s,Dict("id"=>"more","action"=>"step","steps"=>4))
        @test more.ok && all(b->b.process.cells.attempted_steps==7,values(s.branches))
        @test maximum(abs,Array(branch.process.cells.state.inputs).-oldsnapshot.cells.state.inputs)>1e-4
        excess=NativeWorkbench.process_command!(s,Dict("id"=>"excess","action"=>"step","steps"=>4))
        @test !excess.ok && !excess.accepted && all(b->b.process.cells.attempted_steps==7,values(s.branches))
        # Fork only a current complete state; parent and child have independent arrays.
        parent=s.branches["sensing"];parentsnapshot=Interior.device_habitat_snapshot(parent.process)
        fork=NativeWorkbench.process_command!(s,Dict("id"=>"fork-one","action"=>"fork","branch"=>"sensing"))
        @test fork.ok && length(s.branches)==3
        child=s.branches["fork1"]
        @test child.world.time==parent.world.time && child.world.seed==parent.world.seed
        @test child.initial==parent.initial && child.parent.branch=="sensing"
        @test child.process.cells.state.position!==parent.process.cells.state.position
        @test child.process.habitat.arrays.bound!==parent.process.habitat.arrays.bound
        for (owned,original) in zip((child.process.cells.state,child.process.cells.environment,child.process.habitat.arrays,child.process.habitat.feedstock.arrays),
                                     (parent.process.cells.state,parent.process.cells.environment,parent.process.habitat.arrays,parent.process.habitat.feedstock.arrays))
            @test all(k->getproperty(owned,k)!==getproperty(original,k),keys(owned))
        end
        @test child.buffer!==parent.buffer
        @test Interior.device_habitat_snapshot(child.process)==parentsnapshot
        @test Interior.device_habitat_snapshot(parent.process)==parentsnapshot
        checkpointpath=joinpath(s.root,child.parent.checkpoint.file)
        immutable_hash=NativeWorkbench.hashfile(checkpointpath)
        @test immutable_hash==child.parent.checkpoint.sha256
        @test NativeWorkbench.process_command!(s,Dict("id"=>"replay-parent","action"=>"step","branch"=>"sensing","steps"=>1)).ok
        @test NativeWorkbench.process_command!(s,Dict("id"=>"replay-child","action"=>"step","branch"=>"fork1","steps"=>1)).ok
        @test Interior.device_habitat_snapshot(child.process)==Interior.device_habitat_snapshot(parent.process)
        @test NativeWorkbench.process_command!(s,Dict("id"=>"child-scales","action"=>"receptor_scales","branch"=>"fork1","scales"=>[2.,.04,.008])).ok
        @test child.process.habitat.parameters.material_scales!=(parent.process.habitat.parameters.material_scales)
        @test Interior.device_habitat_snapshot(parent.process).parameters==now.parameters
        parentbefore=Interior.device_habitat_snapshot(parent.process)
        @test NativeWorkbench.process_command!(s,Dict("id"=>"changed-child-future","action"=>"step","branch"=>"fork1","steps"=>1)).ok
        @test Interior.device_habitat_snapshot(parent.process)==parentbefore
        @test NativeWorkbench.process_command!(s,Dict("id"=>"parent-future","action"=>"step","branch"=>"sensing","steps"=>1)).ok
        @test child.world.time==parent.world.time
        @test maximum(abs,Array(child.process.cells.state.inputs).-Array(parent.process.cells.state.inputs))>1e-3
        @test !NativeWorkbench.process_command!(s,Dict("id"=>"child-absolute-cap","action"=>"step","branch"=>"fork1","steps"=>2)).accepted
        @test child.process.cells.attempted_steps==9
        @test NativeWorkbench.hashfile(checkpointpath)==immutable_hash
        @test NativeWorkbench.process_command!(s,Dict("id"=>"fork-two","action"=>"fork","branch"=>"off")).ok
        @test length(s.branches)==4 && s.branches["fork2"].parent.branch=="off"
        @test !NativeWorkbench.process_command!(s,Dict("id"=>"too-many-forks","action"=>"fork","branch"=>"sensing")).accepted
        @test_throws ArgumentError NativeWorkbench.validate_command(Dict("id"=>"off-child-scales","action"=>"receptor_scales","branch"=>"fork2","scales"=>[1.,1.,1.]),s)
        @test_throws ArgumentError NativeWorkbench.validate_command(Dict("id"=>"both-fork","action"=>"fork","branch"=>"both"),s)
        @test_throws ArgumentError NativeWorkbench.validate_command(Dict("id"=>"path-fork","action"=>"fork","branch"=>"../old"),s)
        before=filesize(binary)
        stop=NativeWorkbench.process_command!(s,Dict("id"=>"stop","action"=>"stop"))
        @test stop.ok && s.stop_reason=="stop command"
        @test !NativeWorkbench.process_command!(s,Dict("id"=>"after","action"=>"step","steps"=>1)).ok
        @test filesize(binary)==before
        NativeWorkbench.serve!(s)
        finalindex=Interior.JSON3.read(read(joinpath(s.root,"branches","sensing","observatory.json"),String))
        @test finalindex.run_complete
    end
end
@testset "Wide-inventory block means have a uniform accumulator" begin
    g=(dims=(Int32(4),Int32(4),Int32(4)),volume=.125f0)
    for T in (Float32,Float64)
        fields=reshape(T.(1:192)./1000,64,3);feed=(fields=fields.*2,)
        e=(fields=fields,permeability=Float32.(1:64)./64)
        enzymes=reshape(Float32.(1:128)./2000,64,2)
        h=(enzyme_bits=copy(reinterpret(UInt32,enzymes)),bound=enzymes.*3)
        dst=zeros(Float32,8,9)
        NativeWorkbench.habitat_coarse!(CPU(),128)(dst,e,h,feed,g,Int32(2);ndrange=8)
        for coarse in CartesianIndices((2,2,2)),k in 1:9
            sites=[LinearIndices((4,4,4))[(coarse[1]-1)*2+x,(coarse[2]-1)*2+y,(coarse[3]-1)*2+z] for x in 1:2,y in 1:2,z in 1:2]
            at(site)=k<=2 ? Float64(e.fields[site,k])/g.volume : k<=4 ? Float64(feed.fields[site,k-2])/g.volume :
                k<=6 ? Float64(enzymes[site,k-4])/g.volume : k<=8 ? Float64(h.bound[site,k-6])/g.volume : Float64(e.permeability[site])
            expected=Float32(sum(at,sites)/8)
            @test dst[LinearIndices((2,2,2))[coarse],k]≈expected rtol=2eps(Float32)
        end
    end
end

@testset "Partial branches cannot fork" begin
    fake=(branches=Dict("sensing"=>(failed=true,)),limits=NativeWorkbench.Limits())
    @test_throws ErrorException NativeWorkbench.fork_current!(fake,"sensing")
end
