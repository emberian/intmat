@testset "grid and conservative composition" begin
    @test_throws ArgumentError Grid(0,8)
    for dims in ((12,10,1),(8,7,6))
        w=World(Grid(dims...))
        seed_activity!(w)
        initial=material_mass(w)
        for _ in 1:20
            step!(w,(Diffusion(:nutrient,0.01),Excitable(),Deposition(),WaveForcing()),0.005)
        end
        @test w.time ≈ 0.1f0 atol=1f-6
        @test all(isfinite,w.fields.activity)
        @test minimum(w.fields.nutrient)>=0
        @test minimum(w.fields.deposit)>=0
        @test material_mass(w) ≈ initial rtol=2f-6
        @test sum(w.fields.deposit)>0
    end
    a=World(Grid(12,12)); b=World(Grid(12,12))
    seed_activity!(a); seed_activity!(b)
    @test a.fields.activity == b.fields.activity
    @test_throws ArgumentError step!(a,Excitable(),1)
    @test_throws ArgumentError step!(a,Excitable(),-0.1)
end
@testset "diffusion transports without producing matter" begin
    w=World(Grid(12,12,3)); fill!(w.fields.nutrient,0)
    w.fields.nutrient[6,6,2]=1
    advance!(w,Diffusion(:nutrient,0.01),0.01)
    @test sum(w.fields.nutrient) ≈ 1f0 atol=1f-6
    @test 0<w.fields.nutrient[6,6,2]<1
    @test w.fields.nutrient[5,6,2]>0
    @test minimum(w.fields.nutrient)>=0
end
@testset "mineral alters conservative transport" begin
    w=World(Grid(10,9,4)); reference=World(w.grid)
    fill!(w.fields.nutrient,0); fill!(reference.fields.nutrient,0)
    w.fields.nutrient[5,5,2]=1; reference.fields.nutrient[5,5,2]=1
    w.fields.deposit[5,5,2]=4
    advance!(w,Diffusion(:nutrient,0.01;material_feedback=2),0.01)
    advance!(reference,Diffusion(:nutrient,0.01),0.01)
    @test sum(w.fields.nutrient) ≈ 1f0 atol=1f-6
    @test minimum(w.fields.nutrient)>=0
    @test w.fields.nutrient[5,5,2]>reference.fields.nutrient[5,5,2]
    before=copy(w.fields.deposit)
    @test_throws ArgumentError step!(w,(Deposition(),Excitable()),1)
    @test w.fields.deposit==before
    @test_throws ArgumentError step!(w,(Deposition(rate=NaN32),),0.01)
end
