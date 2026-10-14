@testset "acoustic propagation" begin
    w=World(Grid(40,40)); p=AcousticField(w;source_amplitude=0,activity_gain=0,damping=0)
    seed_acoustic!(w,p;radius=.08,dt=.02)
    initial=copy(p.pressure)
    for _ in 1:50; step!(w,p,.02); end
    @test all(isfinite,p.pressure)
    @test maximum(abs.(p.pressure))<1.1
    @test sum(abs,p.pressure[27:31,18:23,1])>10sum(abs,initial[27:31,18:23,1])
    @test maximum(abs.(p.pressure.-initial))>.1
    @test_throws ArgumentError step!(w,p,1)
    @test_throws ArgumentError step!(w,p,.01)
    @test_throws ArgumentError AcousticField(w;speed=-1)
end
@testset "material and acoustic feedback" begin
    a=World(Grid(20,20,12)); b=World(a.grid)
    pa=AcousticField(a;source_amplitude=0,activity_gain=.2,damping=.1)
    pb=AcousticField(b;source_amplitude=0,activity_gain=.2,damping=.1)
    fill!(b.fields.deposit,3f0)
    seed_acoustic!(a,pa;dt=.01); seed_acoustic!(b,pb;dt=.01)
    mass=material_mass(a)
    for _ in 1:40
        step!(a,(pa,Deposition(threshold=.0001)),.01)
        step!(b,pb,.01)
    end
    @test sum(a.fields.deposit)>0
    @test material_mass(a)≈mass rtol=2f-6
    @test maximum(abs.(pa.pressure.-pb.pressure))>.01
    @test sum(a.fields.activity)>0
    @test minimum(pa.envelope)>=0
end
@testset "damped wave loses mechanical energy" begin
    w=World(Grid(24,24)); p=AcousticField(w;source_amplitude=0,activity_gain=0,damping=1)
    seed_acoustic!(w,p;radius=.12,dt=.01)
    function wave_energy()
        kinetic=sum(abs2,(p.pressure.-p.previous)./.01f0)/2
        strain=sum(sum(abs2,(circshift(p.pressure,axis==1 ? (1,0,0) : (0,1,0)).-p.pressure)./w.grid.spacing[axis]) for axis in 1:2)*p.speed^2/2
        kinetic+strain
    end
    initial=wave_energy()
    for _ in 1:200; step!(w,p,.01); end
    @test wave_energy()<.25initial
end
