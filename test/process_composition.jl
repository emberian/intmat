@testset "written scaffold changes neural deposition" begin
    function written_trial(feedback;dt=.01f0)
        w=World(Grid(12,12)); initial=Float64(material_mass(w))
        writer=Interior.ThermalWriter(w;path=:stationary,width=.18,heat=12,
            threshold=.02,conversion=3,start=0,stop=.2)
        neural=NeuralField(w;initial_noise=0,material_feedback=feedback,
            nutrient_feedback=0,diffusion=0,response_time=.2)
        for _ in 1:round(Int,.4/dt)
            step!(w,(writer,neural,Deposition(threshold=0,rate=.2)),dt)
        end
        (w=w,writer=writer,neural=neural,initial=initial)
    end
    a=written_trial(2); b=written_trial(0)
    @test sum(a.writer.dose)>0
    @test sum(a.w.fields.deposit)>0
    @test minimum(a.w.fields.nutrient)>=0
    @test Float64(material_mass(a.w))≈a.initial rtol=2e-6
    @test maximum(abs.(a.neural.latent1.-b.neural.latent1))>1e-3
    @test a.w.time≈.4f0
    # Irreversible exposure persists while the shutter closes and material can dissolve.
    remembered=copy(a.writer.dose); deposited=sum(a.w.fields.deposit)
    step!(a.w,Deposition(rate=0,dissolution=1),.01f0)
    @test a.writer.dose==remembered
    @test sum(a.w.fields.deposit)<deposited
    @test Float64(material_mass(a.w))≈a.initial rtol=2e-6
end

@testset "acoustic excitation reaches conservative particle exchange" begin
    w=World(Grid(16,16)); p=Interior.AcousticField(w;source_amplitude=0,
        activity_gain=.3,damping=.1)
    Interior.seed_acoustic!(w,p;dt=.01f0,radius=.2)
    carrier=ParticleCarrier(reshape(Float32[7.5,7.5,0],3,1))
    initial=Float64(material_mass(w))+Float64(total_particle_nutrient(carrier))
    for _ in 1:20
        step!(w,(p,Deposition(threshold=.0001),carrier),.01f0)
    end
    @test w.fields.pressure==p.pressure
    @test maximum(w.fields.activity)>0
    @test sum(w.fields.deposit)>0
    @test sum(carrier.storage)>0
    @test all(isfinite,carrier.positions)
    @test Float64(material_mass(w))+Float64(total_particle_nutrient(carrier))≈initial rtol=2e-6
    before=deepcopy(w.fields); hidden=copy(carrier.hidden)
    @test_throws ArgumentError step!(w,(Deposition(),carrier,p),.005f0)
    @test w.fields==before
    @test carrier.hidden==hidden
end

@testset "electrical event cadence and field feedback" begin
    function electrical_trial(dt)
        w=World(Grid(13,13,13);seed=31)
        p=Interior.DielectricGrowth(w;event_interval=.1,deposit_per_cell=.5)
        initial=Float64(material_mass(w)); potential=copy(p.potential)
        for _ in 1:round(Int,.4/dt); step!(w,p,dt); end
        (w=w,p=p,initial=initial,potential=potential)
    end
    a=electrical_trial(.1f0); b=electrical_trial(.025f0)
    @test a.p.events==b.p.events==4
    @test a.p.occupied==b.p.occupied
    @test a.p.potential==b.p.potential
    @test a.p.potential!=a.potential
    @test a.p.last_residual<=a.p.tolerance
    @test Float64(material_mass(a.w))≈a.initial rtol=2e-6
    @test all(isfinite,a.w.fields.activity)
    # Exchange after electrical events spends material without creating inventory.
    step!(a.w,Deposition(threshold=0),.01f0)
    @test Float64(material_mass(a.w))≈a.initial rtol=2e-6
end
