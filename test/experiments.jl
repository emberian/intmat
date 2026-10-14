@testset "experiment schedule and detached observations" begin
    function experiment_trial(cadence)
        w=World(Grid(12,12)); seed_activity!(w)
        p=NeuralField(w;diffusion=0)
        e=Experiment(w,(neural=p,reaction=Deposition(),transport=Diffusion(:nutrient,.001));dt=.01)
        r=run!(e;steps=17,capture_every=cadence,
            observe=(w,ps,k)->(activity=w.fields.activity,latent=ps[1].latent1))
        (w=w,p=p,r=r)
    end
    a=experiment_trial(1); b=experiment_trial(6)
    @test a.w.fields==b.w.fields
    @test a.p.latent1==b.p.latent1
    @test a.w.time==b.w.time
    @test getproperty.(b.r.samples,:step)==[0,6,12,17]
    @test b.r.metadata.process_order==(:neural,:reaction,:transport)
    saved=copy(b.r.samples[end].observation.activity)
    fill!(b.w.fields.activity,99)
    @test b.r.samples[end].observation.activity==saved
    @test b.r.samples[1].observation.activity!=saved
    e=Experiment(World(Grid(4,4)),();dt=.01)
    @test length(run!(e;steps=0).samples)==1
    @test_throws ArgumentError run!(e;steps=-1)
    @test_throws ArgumentError run!(e;steps=1,capture_every=0)
end

@testset "experiment preflight and explicit publishing order" begin
    w=World(Grid(8,8)); p=NeuralField(w)
    before=deepcopy(w.fields)
    @test_throws ArgumentError Experiment(w,(excitable=Excitable(),neural=p);dt=.01)
    @test w.fields==before && w.time==0
    e=Experiment(w,(excitable=Excitable(),neural=p);dt=.01,allow_overwrite=true)
    @test e.names==(:excitable,:neural)
    @test_throws ArgumentError Experiment(w,(reaction=Deposition(),bad=Diffusion(:absent,1));dt=.01)
    @test w.fields==before && w.time==0
    @test_throws ArgumentError Experiment(w,(:same=>Deposition(),:same=>Deposition());dt=.01)
    @test process_ports(p).publishes==(:activity,)
end

@testset "writer feedback ablation through shared experiment" begin
    function written_experiment(feedback)
        w=World(Grid(12,12))
        writer=ThermalWriter(w;path=:stationary,width=.18,heat=12,
            threshold=.02,conversion=3,start=0,stop=.2)
        neural=NeuralField(w;initial_noise=0,material_feedback=feedback,
            nutrient_feedback=0,diffusion=0,response_time=.2)
        e=Experiment(w,(writer=writer,neural=neural,
            reaction=Deposition(threshold=0,rate=.2),
            transport=Diffusion(:nutrient,.001));dt=.01)
        r=run!(e;steps=40,capture_every=13,
            observe=(w,ps,k)->(dose=sum(ps[1].dose),latent=ps[2].latent1))
        (e=e,r=r)
    end
    coupled=written_experiment(2); ablated=written_experiment(0)
    @test coupled.r.samples[end].observation.dose>0
    @test coupled.r.samples[end].diagnostics.deposit_mass>0
    @test maximum(abs.(coupled.e.processes[2].latent1.-ablated.e.processes[2].latent1))>1e-3
    @test coupled.r.samples[end].diagnostics.material_mass≈144 rtol=2e-6
    @test all(s->s.diagnostics.finite,coupled.r.samples)
    optical=OpticalMap((w,ps)->(density=copy(w.fields.deposit),activity=copy(w.fields.activity)))
    t=coupled.e.world.time
    @test size(optical(coupled.e).density)==(12,12,1)
    @test coupled.e.world.time==t
end
