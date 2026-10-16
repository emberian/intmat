# Standalone bounded CPU composition gate (two Julia threads are sufficient).
using Test, Interior, TOML, Serialization
include(joinpath(@__DIR__,"..","scripts","native_workbench_worker.jl"))
const NW=NativeWorkbench
@testset "Living resident complete forks" begin
    cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments/native/living-patch-recruitment-16.toml"))
    cfg["cellular"]["founders"]=4;cfg["cellular"]["capacity"]=256
    mktempdir() do dir
        s=NW.create_session(cfg,joinpath(dir,"live");side=4,living=true,
            limits=NW.Limits(max_steps=256,max_model_time=8.,wall_seconds=3600.))
        b=s.branches["sensing"]
        @test b.policy!==nothing
        @test b.policy.inheritance_fraction==.5f0
        @test s.config["living_policy"]["horizon"]==2.
        @test s.config["living_policy"]["actor_rate"]==.003
        @test s.config["living_policy"]["noise_seed"]==6201
        @test b.policy.options.seed==UInt32(6201)
        @test b.policy.options.sigma==.3f0
        # Reach pending delayed cohorts through real native steps.
        for _ in 1:12;NW.advance_one!(b,s.dt);end
        @test any(!iszero,Array(b.policy.state.queue_ids))
        # Seed existing learned private state to test preservation even before the
        # first two-second paid maturity; this is a checkpoint fixture, not a
        # claim that learning has occurred during these twelve steps.
        b.policy.state.learned[1,1,1]=.125f0
        # Exercise the actual funded birth-copy transaction on a declared
        # two-resident birth-transition fixture. Native division geometry is not
        # being asserted by this checkpoint composition test.
        oldh=b.policy
        b.policy=typeof(oldh)(oldh.base,oldh.events,.5f0,oldh.birth_transfer_adjusted)
        h=b.policy;c=b.process.cells;st=c.state;q=c.parameters
        st.new_parent[2]=1
        view=(active=st.active,ids=st.ids,new_parent=st.new_parent,energy=st.energy,metabolites=st.metabolites,ledger=st.ledger,ledger_low=st.ledger_low)
        chemistry=(efficiency=q.efficiency,potentials=q.potentials,kappa=h.options.kappa)
        NW.LivingLaw.L.policy_launch!(h,NW.LivingLaw.birth!,view,c.environment.status,chemistry,h.state,h.events,b.world.time,h.inheritance_fraction,h.birth_transfer_adjusted;ndrange=Int(c.geometry.capacity))
        st.new_parent[2]=0
        @test h.events.copy_work[Int(st.ids[2])]>0
        @test h.state.paid[2]>0
        @test h.state.learned[2,1,1]==.0625f0
        s.sequence=1
        NW.fork_current!(s,"sensing")
        f=s.branches["fork1"]
        @test NW.living_snapshot(b.policy)==NW.living_snapshot(f.policy)
        for key in keys(b.policy.state)
            @test getproperty(b.policy.state,key)!==getproperty(f.policy.state,key)
        end
        for key in keys(b.policy.events)
            @test getproperty(b.policy.events,key)!==getproperty(f.policy.events,key)
        end
        @test b.policy.original_masks!==f.policy.original_masks
        @test f.policy.habitat===f.process.habitat
        @test b.policy.steps==12
        before=Interior.device_habitat_snapshot(b.process)
        NW.advance_one!(f,s.dt)
        @test Interior.device_habitat_snapshot(b.process)==before
        NW.advance_one!(b,s.dt)
        @test Interior.device_habitat_snapshot(b.process)==Interior.device_habitat_snapshot(f.process)
        @test NW.living_snapshot(b.policy)==NW.living_snapshot(f.policy)
        old=NW.living_snapshot(f.policy)
        s.sequence=2
        NW.intervene_receptors!(s,"fork1","receptors",(0.5f0,1f0,0.5f0))
        @test NW.living_snapshot(f.policy)==old
        @test f.policy.habitat===f.process.habitat
        @test b.process.habitat.parameters.material_scales!=f.process.habitat.parameters.material_scales
        payload=deserialize(joinpath(s.root,f.checkpoint.file))
        @test payload.living==old
        bad=merge(old,(source=Dict("bad"=>"hash"),))
        @test_throws ErrorException NW.living_restore(f.process,bad)
    end
end
