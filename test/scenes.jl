include(joinpath(@__DIR__,"..","studies","verify.jl"))

@testset "composed scenes" begin
    for dimension in (2,3)
        config=SceneVerification.configuration(;dimension,resolution=8,particles=12)
        a=build_scene(config);b=build_scene(config)
        @test a.world.fields.activity==b.world.fields.activity
        @test a.particles.positions==b.particles.positions
        @test a.particles.hidden==b.particles.hidden
        @test size(a.world.fields.nutrient)==(8,8,dimension==2 ? 1 : 8)
        initial=SceneVerification.inventory(a)
        for _ in 1:50
            step!(a.world,a.processes,0.02)
        end
        @test SceneVerification.inventory(a) ≈ initial rtol=1e-5
        @test all(v->all(isfinite,v),values(a.world.fields))
        @test minimum(a.world.fields.nutrient)>=0
        @test minimum(a.world.fields.deposit)>=0
        @test sum(a.world.fields.deposit)>0
        @test total_particle_nutrient(a.particles)>0
        @test a.particles.positions!=b.particles.positions
        @test a.particles.hidden!=b.particles.hidden

        no_particles=deepcopy(config);no_particles["particles"]["count"]=0
        removed=build_scene(no_particles)
        @test removed.particles===nothing
        @test removed.world.fields.activity==b.world.fields.activity
        @test removed.world.fields.nutrient==b.world.fields.nutrient
        @test length(removed.processes)==length(b.processes)-1
        before_hidden=copy(b.particles.hidden)
        SceneVerification.ablate!(b,"no_neural_secretion")
        @test b.particles.hidden==before_hidden
        @test all(iszero,b.particles.rule.secretion)
    end
end
