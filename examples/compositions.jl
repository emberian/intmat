using Interior

# One schedule drives writing, response, conversion and redistribution.
function written_neural(;feedback=2f0,steps=80,capture_every=20)
    w=World(Grid(24,24);seed=31)
    writer=ThermalWriter(w;path=:stationary,width=.18,heat=12,
        threshold=.02,conversion=3,start=0,stop=.2)
    neural=NeuralField(w;initial_noise=0,material_feedback=feedback,
        nutrient_feedback=0,diffusion=0,response_time=.2)
    experiment=Experiment(w,(writer=writer,neural=neural,
        reaction=Deposition(threshold=0,rate=.2),
        transport=Diffusion(:nutrient,.001));dt=.01)
    report=run!(experiment;steps,capture_every,
        observe=(w,ps,k)->(dose=sum(ps[1].dose),latent=copy(ps[2].latent1)))
    (;experiment,report)
end

# Neural reads the pressure published earlier in each complete step.
function acoustic_neural(;feedback=.8f0,steps=80,capture_every=20)
    w=World(Grid(24,24);seed=31)
    acoustic=AcousticField(w;source_amplitude=0,activity_gain=0,damping=.1)
    seed_acoustic!(w,acoustic;dt=.01,radius=.2)
    neural=NeuralField(w;initial_noise=0,pressure_feedback=feedback,
        nutrient_feedback=0,material_feedback=0,diffusion=0)
    experiment=Experiment(w,(acoustic=acoustic,neural=neural,
        reaction=Deposition(threshold=0),transport=Diffusion(:nutrient,.001));dt=.01)
    report=run!(experiment;steps,capture_every)
    (;experiment,report)
end

if abspath(PROGRAM_FILE)==@__FILE__
    coupled=written_neural(); ablated=written_neural(feedback=0)
    println("Writer feedback changes neural state: ",
        maximum(abs.(coupled.experiment.processes[2].latent1.-ablated.experiment.processes[2].latent1)))
    wave=acoustic_neural(); silent=acoustic_neural(feedback=0)
    println("Pressure feedback changes neural state: ",
        maximum(abs.(wave.experiment.processes[2].latent1.-silent.experiment.processes[2].latent1)))
    println("Final writing diagnostics: ",coupled.report.samples[end].diagnostics)
end
