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
            limits=NativeWorkbench.Limits(max_steps=8,max_steps_per_command=4,max_model_time=1.,wall_seconds=3600.,idle_seconds=60.,max_commands=20))
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
        @test_throws ArgumentError NativeWorkbench.reserve!(s,s.limits.artifact_bytes)
        more=NativeWorkbench.process_command!(s,Dict("id"=>"more","action"=>"step","steps"=>4))
        @test more.ok && all(b->b.process.cells.attempted_steps==7,values(s.branches))
        excess=NativeWorkbench.process_command!(s,Dict("id"=>"excess","action"=>"step","steps"=>2))
        @test !excess.ok && !excess.accepted && all(b->b.process.cells.attempted_steps==7,values(s.branches))
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
