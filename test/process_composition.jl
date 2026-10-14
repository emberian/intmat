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

@testset "written dose changes a frozen optical path" begin
    w=World(Grid(12,12,12)); writer=Interior.ThermalWriter(w;path=:stationary,
        width=.25,heat=16,conversion=0,threshold=0,stop=.4)
    for _ in 1:40; step!(w,writer,.01f0); end
    before=deepcopy(w.fields); dose=copy(writer.dose); time=w.time
    # Dose-to-index is an explicit optical interpretation of accumulated exposure.
    lens=Interior.GradedIndex(1f0 .+ .4f0.*writer.dose;extent=w.grid.extent)
    ray=Interior.trace_grin(lens,(.2,0.,-1.),(0.,0.,1.);ds=.01)
    @test ray.status==:reached
    @test ray.position[1]<.2-1e-4
    @test ray.optical_length>ray.arc_length
    @test w.fields==before
    @test writer.dose==dose
    @test w.time==time
end

@testset "preflight compares the timestep the kernel receives" begin
    w=World(Grid(5,7,1;extent=(10,7,9))); fill!(w.fields.activity,1)
    phase=Interior.PhaseSeparation(w;mobility=0,activity_relaxation=0,nutrient_diffusion=2)
    limit=stability_limit(w,phase)
    @test Float32(limit)>limit
    before=deepcopy(w.fields)
    @test_throws ArgumentError step!(w,(Deposition(),phase),limit)
    @test w.fields==before
    @test w.time==0
end

@testset "mechanical carrier keeps separate history in a shared experiment" begin
    w=World(Grid(8,8)); fill!(w.fields.activity,.5f0)
    graph=Interior.BondNetwork(w;size=(3,3,1),jitter=0,stretch_rate=.05)
    e=Interior.Experiment(w,(solid=graph,conversion=Deposition());dt=.001f0)
    initial=Float64(material_mass(w)); reference=copy(graph.reference)
    report=Interior.run!(e;steps=8,capture_every=4,
        observe=(w,ps,k)->Interior.mechanical_optical_fields(w,ps[1]))
    @test graph.elapsed≈Float64(w.time) rtol=1e-6
    @test graph.position!=reference
    @test graph.reference==reference
    @test all(s->isapprox(s.diagnostics.material_mass,initial;rtol=2e-6),report.samples)
    before=deepcopy(w.fields); position=copy(graph.position); elapsed=graph.elapsed
    graph.load=zeros(2,length(graph.node_ids))
    @test_throws DimensionMismatch step!(w,(Deposition(),graph),.001f0)
    @test w.fields==before
    @test graph.position==position
    @test graph.elapsed==elapsed
end

@testset "shared diagnostics include every material carrier" begin
    w=World(Grid(8,8)); p=ParticleCarrier(reshape(Float32[3.5,3.5,0],3,1))
    ecology=Interior.Ecology(reshape(Float32[3.5,3.5,0],3,1);capacity=4)
    initial=Float64(material_mass(w))+sum(Float64,p.storage)+sum(Float64,ecology.storage)+ecology.waste
    e=Interior.Experiment(w,(mobile=p,colony=ecology,conversion=Deposition());dt=.01f0)
    report=Interior.run!(e;steps=8,capture_every=4)
    @test all(s->isapprox(s.diagnostics.material_mass,initial;rtol=2e-6),report.samples)
    @test report.samples[end].diagnostics.owned_mass>0
    @test all(s->s.diagnostics.finite,report.samples)
end

@testset "phase interfaces meet a finite-material ecology" begin
    w=World(Grid(8,8)); phase=Interior.PhaseSeparation(w;bulk=1,kappa=.001,
        mobility=.01,nutrient_diffusion=.01,interface_activity=.8)
    ecology=Interior.Ecology(Float32[3.5 3.5;3.5 3.5;0 0];seed=11,
        speed=.05,lifespan=10,capacity=12)
    initial=Interior.ecology_mass(w,ecology); composition=sum(Float64,phase.phi)
    for _ in 1:30
        step!(w,(phase,ecology,Deposition(threshold=0,rate=.1)),.01f0)
    end
    @test Interior.ecology_mass(w,ecology)≈initial rtol=2e-6
    @test sum(Float64,phase.phi)≈composition atol=1e-6
    @test minimum(w.fields.nutrient)>=0
    @test minimum(w.fields.deposit)>=0
    @test sum(w.fields.deposit)>0
    @test all(isfinite,phase.phi)
    @test all(isfinite,ecology.hidden)
    # Invalid mutable controls in a later process must fail before deposition.
    ecology.uptake_rate=-1
    before=deepcopy(w.fields)
    @test_throws ArgumentError step!(w,(Deposition(),ecology),.01f0)
    @test w.fields==before
end

@testset "experiment observations preserve integration" begin
    function acoustic_experiment()
        w=World(Grid(12,12)); acoustic=Interior.AcousticField(w;source_amplitude=0)
        Interior.seed_acoustic!(w,acoustic;dt=.01f0,radius=.2)
        Interior.Experiment(w,(sound=acoustic,conversion=Deposition(threshold=.0001));dt=.01f0)
    end
    a=acoustic_experiment(); b=acoustic_experiment()
    ra=Interior.run!(a;steps=12,capture_every=3,
        observe=(w,ps,k)->w.fields.activity)
    rb=Interior.run!(b;steps=12,capture_every=5)
    @test a.world.fields==b.world.fields
    @test a.world.time==b.world.time≈.12f0
    @test map(s->s.step,ra.samples)==[0,3,6,9,12]
    @test map(s->s.step,rb.samples)==[0,5,10,12]
    @test all(iszero,ra.samples[1].observation)
    @test ra.samples[end].observation==a.world.fields.activity
    @test ra.samples[end].observation!==a.world.fields.activity
    before=deepcopy(a.world.fields)
    @test_throws ArgumentError Interior.Experiment(a.world,a.processes;dt=.005f0)
    @test a.world.fields==before
    another=Interior.run!(a;steps=0)
    @test length(another.samples)==1
    @test another.metadata.start_time==another.metadata.end_time
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
