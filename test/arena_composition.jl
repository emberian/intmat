# Focused counterexamples for the coupled compartment environment.
@testset "active compartments conserve inventories in 3D and thin grids" begin
    for dims in ((7,5,4),(1,2,1),(1,1,1))
        w=World(Grid(dims...;extent=Float32.(dims)))
        phase=Interior.ImplicitPhaseSeparation(w;mineral_pinning=0,
            nutrient_diffusion=0,activity_relaxation=0,interface_activity=0)
        medium=Interior.ActiveMedium(w,phase;nutrient_diffusion=.3,signal_diffusion=.2)
        w.fields.nutrient .= reshape(Float32.(1:prod(dims)),dims)
        medium.signal[1]=1
        w.fields.deposit[1]=2
        mass=sum(Float64,w.fields.nutrient)+sum(Float64,w.fields.deposit)
        composition=sum(Float64,phase.phi)
        medium.signal_source[end]=.2
        medium.remodeling[1]=.5
        for _ in 1:15
            step!(w,(medium,phase),.05f0)
        end
        @test sum(Float64,w.fields.nutrient)+sum(Float64,w.fields.deposit)≈mass rtol=3e-6
        @test sum(Float64,phase.phi)≈composition atol=3e-6
        @test minimum(w.fields.nutrient)>=0
        @test all(isfinite,medium.permeability)
        @test all(x->0<=x<=1,medium.permeability)
        @test all(iszero,medium.signal_source)
        @test all(iszero,medium.remodeling)
    end
end

@testset "funded membrane changes exchange; signals reopen it" begin
    function transport_trial(deposit,signal)
        w=World(Grid(7,1,1;extent=(7,1,1)))
        phase=Interior.ImplicitPhaseSeparation(w;noise=0,mineral_pinning=0,
            nutrient_diffusion=0,activity_relaxation=0,interface_activity=0)
        medium=Interior.ActiveMedium(w,phase;nutrient_diffusion=.2,signal_diffusion=0,signal_decay=0)
        fill!(w.fields.nutrient,0); w.fields.nutrient[4]=1
        fill!(w.fields.deposit,deposit); fill!(medium.signal,signal)
        step!(w,medium,.1f0)
        w.fields.nutrient[4]
    end
    open=transport_trial(0,0); closed=transport_trial(3,0); signaled=transport_trial(3,10)
    @test open<signaled<closed
end

@testset "empty ecology remains a valid finite inventory carrier" begin
    w=World(Grid(2,2,2;extent=(2,2,2)))
    colony=Interior.NeuralEcology(zeros(Float32,3,0);capacity=0)
    initial=ecology_mass(w,colony)
    step!(w,colony,.1f0)
    @test isempty(colony.ids)
    @test ecology_mass(w,colony)==initial
    @test ecology_statistics(colony).energy_error==0
end

@testset "coupled checkpoint preserves aliases and deterministic continuation" begin
    w=World(Grid(6,5,3;extent=(6,5,3)))
    phase=Interior.ImplicitPhaseSeparation(w;mineral_pinning=0,
        nutrient_diffusion=0,activity_relaxation=0,interface_activity=0)
    medium=Interior.ActiveMedium(w,phase)
    colony=Interior.NeuralEcology(Float32[1 2; 1 2; 0 1];medium=medium,storage=2,energy=5,capacity=12,mutation=.05)
    # Force reproducible birth events, then checkpoint the actual FFT-backed process.
    fill!(colony.weights,0); colony.weights[size(colony.hidden,1)+8,end,:].=2
    a=Interior.Arena(w,(colony,medium,phase,Interior.Chemostat());dt=.02)
    for _ in 1:3; Interior.arena_step!(a); end
    mktempdir() do dir
        path=joinpath(dir,"restart.jls")
        Interior.checkpoint!(a,path)
        b=Interior.restore_arena(path)
        @test b.processes[2].membrane===b.world.fields.deposit
        @test b.processes[2].phase===b.processes[3]
        @test b.processes[1].medium===b.processes[2]
        for _ in 1:8; Interior.arena_step!(a); Interior.arena_step!(b); end
        @test a.world.fields==b.world.fields
        @test a.world.time==b.world.time
        @test a.processes[1].ids==b.processes[1].ids
        @test a.processes[1].weights==b.processes[1].weights
        @test a.processes[1].plastic_marks==b.processes[1].plastic_marks
        @test a.processes[1].hidden==b.processes[1].hidden
        @test a.processes[1].positions==b.processes[1].positions
        @test a.processes[3].phi==b.processes[3].phi
        snap=Interior.arena_snapshot(a;max_side=3,max_agents=1)
        @test length(snap["agents"])<=1
        @test all(f->all(<=(3),f["shape"]),values(snap["fields"]))
        @test abs(snap["metrics"]["material_error"])<1e-4
    end
end

@testset "inhabitant cues use physical time and require funded work" begin
    function cue_trial(dt;work_rate=.1)
        w=World(Grid(1,1,1))
        phase=Interior.ImplicitPhaseSeparation(w;noise=0,mineral_pinning=0,
            nutrient_diffusion=0,activity_relaxation=0,interface_activity=0)
        medium=Interior.ActiveMedium(w,phase;signal_diffusion=0,signal_decay=0,forcing_relaxation=0)
        p=Interior.NeuralEcology(zeros(Float32,3,1);medium=medium,uptake_rate=0,
            basal_cost=0,work_rate=work_rate,secretion_rate=0,speed=0,mutation=0)
        fill!(p.weights,0); p.weights[size(p.hidden,1)+9,end,1]=1
        for _ in 1:round(Int,.1/dt); step!(w,(p,medium),dt); end
        (signal=sum(medium.signal),forcing=sum(medium.forcing))
    end
    a=cue_trial(.1); b=cue_trial(.01)
    @test a.signal≈b.signal rtol=1e-5
    @test a.forcing≈b.forcing rtol=1e-5
    unfunded=cue_trial(.1;work_rate=0)
    @test unfunded.signal==unfunded.forcing==0
end

@testset "death returns cargo to final periodic cell after crossing" begin
    w=World(Grid(4,2,1;extent=(4,2,1))); fill!(w.fields.nutrient,0)
    p=Interior.NeuralEcology(reshape(Float32[3.9,.2,0],3,1);storage=2,energy=5,
        uptake_rate=0,basal_cost=0,work_rate=.1,secretion_rate=0,speed=2,lifespan=.1,
        plasticity_rate=0)
    fill!(p.weights,0); p.weights[size(p.hidden,1)+1,end,1]=4
    mass=ecology_mass(w,p)
    step!(w,p,.1f0)
    @test isempty(p.ids)
    @test w.fields.nutrient[1,1,1]==2
    @test w.fields.nutrient[4,1,1]==0
    @test ecology_mass(w,p)==mass
    @test abs(ecology_statistics(p).energy_error)<1e-6
end

@testset "closed fuel bounds cumulative work through reproduction and death" begin
    w=World(Grid(1,1,1)); fill!(w.fields.nutrient,4)
    p=Interior.NeuralEcology(zeros(Float32,3,1);storage=1,energy=3,capacity=8,
        uptake_rate=4,basal_cost=.2,work_rate=.1,secretion_rate=.1,speed=0,
        birth_material=.15,birth_energy=.1,mutation=0,lifespan=.8,
        catabolic_fraction=.5,energy_yield=2,plasticity_rate=0)
    fill!(p.weights,0); h=size(p.hidden,1)
    p.weights[h+4,end,:].=3; p.weights[h+8,end,:].=3
    initial_mass=ecology_mass(w,p); initial_work=p.energy_input
    errors=zeros(Float64,4)
    for _ in 1:500
        step!(w,p,.05f0)
        errors[1]=max(errors[1],abs(ecology_mass(w,p)-initial_mass))
        errors[2]=max(errors[2],abs(p.energy_input-initial_work-p.energy_yield*p.waste))
        errors[3]=max(errors[3],abs(sum(Float64,p.energy;init=0.)+p.heat-p.energy_input))
        errors[4]=max(errors[4],p.energy_input-initial_work-p.energy_yield*initial_mass)
    end
    @test errors[1]<3e-5
    @test errors[2]<3e-5
    @test errors[3]<3e-5
    @test errors[4]<3e-5
    @test p.births>0 && p.deaths>0
    @test isempty(p.ids)
    exhausted_work=p.energy_input; exhausted_heat=p.heat
    for _ in 1:20; step!(w,p,.05f0); end
    @test p.energy_input==exhausted_work && p.heat==exhausted_heat
    # Pure assimilation cannot power work, even while uptake stores matter.
    v=World(Grid(1,1,1))
    q=Interior.NeuralEcology(zeros(Float32,3,1);storage=0,energy=0,
        catabolic_fraction=0,energy_yield=100,basal_cost=0,work_rate=.1)
    before=ecology_mass(v,q); step!(v,q,.1f0)
    @test q.energy_input==q.heat==q.waste==0
    @test isempty(q.ids)
    @test ecology_mass(v,q)≈before atol=1e-6
end
