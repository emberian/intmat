@testset "long-run clock preserves solver time" begin
    w=World(Grid(1,1)); dt=1f-6; steps=1_000_000
    for _ in 1:steps; step!(w,(),dt); end
    @test w.time isa Float64
    @test w.time==steps*Float64(dt)
    @test eltype(w.fields.activity)==Float32
    # Small solver steps still move a clock far beyond Float32's resolution.
    w.time=65536.0; initial=w.time
    for _ in 1:1000; step!(w,(),dt); end
    @test w.time>initial
    @test w.time≈initial+1000Float64(dt) atol=1e-8 rtol=0
    # The public Real argument is rounded to solver precision before accumulation.
    w=World(Grid(2,2)); requested=0.001234567890123
    step!(w,(WaveForcing(),Deposition()),requested)
    @test w.time==Float64(Float32(requested))
    @test all(isfinite,w.fields.pressure)
    a=World(Grid(4,4)); b=World(a.grid); a.time=b.time=.25
    step!(a,WaveForcing(frequency=.7f0),.001f0)
    step!(b,WaveForcing(frequency=-.7f0),.001f0)
    @test a.fields.pressure≈b.fields.pressure
    @test a.fields.vx≈-b.fields.vx
    @test maximum(abs,a.fields.vx)>0
end
