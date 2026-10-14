@testset "focused thermal writing" begin
    w=World(Grid(12,12,12))
    p=ThermalWriter(w;path=:stationary,width=.1,heat=10,threshold=0)
    step!(w,p,.01)
    @test p.temperature≈reverse(p.temperature;dims=1)
    @test p.temperature≈reverse(p.temperature;dims=2)
    @test p.temperature≈reverse(p.temperature;dims=3)
    @test minimum(p.temperature)>=0
    @test material_mass(w)≈prod(w.grid.size)
    @test minimum(w.fields.nutrient)>=0
    @test maximum(w.fields.deposit)<=1
    before=copy(p.dose)
    for _ in 1:30;step!(w,p,.01);end
    @test all(p.dose.>=before)
    @test sum(w.fields.deposit)>0
    @test material_mass(w)≈prod(w.grid.size)
    @test_throws ArgumentError step!(w,p,2stability_limit(w,p))
    q=ThermalWriter(w;path=:lissajous)
    @test beam_position(q,0)!=beam_position(q,.2)
    @test collect(beam_position(q,0))≈collect(beam_position(q,q.period))
    old=copy(w.fields.deposit)
    mismatch=ThermalWriter(World(Grid(3,3)))
    @test_throws ArgumentError step!(w,(p,mismatch),.01)
    @test w.fields.deposit==old
    @test_throws ArgumentError ThermalWriter(w;width=0)
    # Closed shutter still cools the specimen and preserves its memory.
    closed=ThermalWriter(w;stop=0,threshold=100)
    fill!(closed.temperature,1)
    fill!(closed.dose,2)
    step!(w,closed,.01)
    @test maximum(closed.temperature)<1
    @test all(==(2),closed.dose)
end
