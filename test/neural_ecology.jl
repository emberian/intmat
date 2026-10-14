@testset "NeuralEcology inherited controller and material ledgers" begin
    w=World(Grid(8,8;extent=(8,8,1)))
    p=NeuralEcology(Float32[1 3;1 3;0 0];storage=2,energy=5,capacity=8,mutation=.1)
    h=size(p.hidden,1); p.weights[h+8,:,:].=0; p.weights[h+8,end,:].=2
    mass=ecology_mass(w,p); genotypes=copy(p.weights)
    step!(w,p,.1)
    @test p.births==2 && p.parents[3:4]==[1,2] && p.generations[3:4]==[1,1]
    @test p.weights[:,:,1:2]==genotypes
    @test p.weights[:,:,3]!=p.weights[:,:,1]
    @test all(iszero,p.plastic_marks[:,:,3:4]) && all(iszero,p.hidden[:,3:4])
    @test any(!iszero,p.plastic_marks[:,:,1:2])
    @test ecology_mass(w,p)≈mass rtol=1e-6
    @test abs(ecology_statistics(p).energy_error)<2e-6
    q=deepcopy(p); v=deepcopy(w); reorder_ecology!(q,reverse(eachindex(q.ids)))
    step!(w,p,.1); step!(v,q,.1)
    @test p.ids==q.ids && p.positions==q.positions && p.weights==q.weights && p.plastic_marks==q.plastic_marks
    @test w.fields.nutrient==v.fields.nutrient && w.fields.deposit==v.fields.deposit
    p.age.=p.lifespan
    before=ecology_mass(w,p); step!(w,p,.1)
    @test isempty(p.ids) && p.deaths>0
    @test ecology_mass(w,p)≈before rtol=1e-6
    @test abs(ecology_statistics(p).energy_error)<1e-5
end
@testset "NeuralEcology zero mutation and plasticity separation" begin
    w=World(Grid(3,3)); p=NeuralEcology(zeros(Float32,3,1);storage=2,energy=5,mutation=0,capacity=2)
    h=size(p.hidden,1); p.weights[h+8,:,:].=0; p.weights[h+8,end,:].=3
    p.plastic_marks.=.2
    genotype=copy(p.weights[:,:,1]); step!(w,p,.1)
    @test p.weights[:,:,2]==genotype && p.weights[:,:,1]==genotype
    @test all(iszero,p.plastic_marks[:,:,2])
    @test all(iszero,p.hidden[:,2])
    before=copy(w.fields.nutrient); p.plasticity_rate=-1
    @test_throws ArgumentError step!(w,p,.1)
    @test w.fields.nutrient==before
end
function neural_test_medium(dims;permeability=1f0)
    (phi=zeros(Float32,dims),interface=zeros(Float32,dims),permeability=fill(Float32(permeability),dims),signal=zeros(Float32,dims),signal_source=zeros(Float32,dims),remodeling=zeros(Float32,dims))
end
function neural_fixed_controller(positions;kw...)
    p=NeuralEcology(positions;uptake_rate=0,basal_cost=0,work_rate=.1,secretion_rate=0,plasticity_rate=0,energy=5,kw...)
    fill!(p.weights,0); h=size(p.hidden,1)
    p.weights[h+8,end,:].=-3 # No spontaneous division in movement checks.
    p
end
@testset "NeuralEcology directional 3D sensing and unique periodic neighbors" begin
    w=World(Grid(8,8,8;extent=(8,8,8))); fill!(w.fields.nutrient,0)
    w.fields.nutrient[3,3,4]=4
    p=neural_fixed_controller(reshape(Float32[2.1,2.1,2.1],3,1)); h=size(p.hidden,1)
    p.weights[h+3,h+12,1]=4 # nutrient central z difference
    before=p.positions[3,1]; step!(w,p,.2)
    @test p.positions[3,1]>before && p.last_actions[3,1]>.5
    # A neighbor above contributes a positive z cue, including across periodic seam.
    q=neural_fixed_controller(Float32[2.1 2.1;2.1 2.1;7.7 .2]); q.weights[h+3,h+31,:].=4
    step!(w,q,.1)
    @test q.last_actions[3,1]>0 && q.last_actions[3,2]<0
    tiny=World(Grid(1,1)); r=neural_fixed_controller(zeros(Float32,3,2)); r.weights[h+1,h+28,:].=1
    step!(tiny,r,.1)
    @test r.last_actions[1,1]≈tanh(.5f0) # one neighbor, never repeated periodic aliases
    @test all(iszero,r.positions[3,:])
end
@testset "NeuralEcology environmental interventions alter policy and movement" begin
    w=World(Grid(6,6;extent=(6,6,1))); med=neural_test_medium(w.grid.size)
    p=neural_fixed_controller(reshape(Float32[2.2,2.2,0],3,1);medium=med,speed=1)
    h=size(p.hidden,1); p.weights[h+1,h+3,1]=3 # local signal controls x movement
    q=deepcopy(p); v=deepcopy(w); q.medium.signal[3,3,1]=2
    step!(w,p,.2); step!(v,q,.2)
    @test q.last_actions[1,1]>p.last_actions[1,1] && q.positions[1,1]>p.positions[1,1]
    open=neural_fixed_controller(reshape(Float32[2.2,2.2,0],3,1);medium=neural_test_medium(w.grid.size),speed=1)
    open.weights[h+1,end,1]=3
    blocked=deepcopy(open); blocked.medium.permeability.=0
    a=deepcopy(w); b=deepcopy(w); step!(a,open,.2); step!(b,blocked,.2)
    @test open.positions[1,1]>2.2f0 && blocked.positions[1,1]==2.2f0
    # Very high requested speed still cannot hop several cells over a thin barrier.
    fast=neural_fixed_controller(reshape(Float32[2.2,2.2,0],3,1);speed=1e5)
    fast.weights[h+1,end,1]=3; step!(w,fast,1)
    @test fast.positions[1,1]-2.2f0<=.45001f0
    slow=neural_fixed_controller(reshape(Float32[2.2,2.2,0],3,1);speed=1e5,
        medium=neural_test_medium(w.grid.size;permeability=.1))
    slow.weights[h+1,end,1]=3; step!(w,slow,1)
    @test slow.positions[1,1]-2.2f0<=.04501f0
end
@testset "NeuralEcology membrane deposition erosion and funded signaling" begin
    w=World(Grid(4,4;extent=(4,4,1))); med=neural_test_medium(w.grid.size)
    p=NeuralEcology(reshape(Float32[1,1,0],3,1);medium=med,storage=2,energy=5,
        uptake_rate=0,basal_cost=0,secretion_rate=1,plasticity_rate=0,speed=0)
    h=size(p.hidden,1); fill!(p.weights,0)
    p.weights[h+5,end,1]=5; p.weights[h+6,end,1]=-10
    p.weights[h+7,end,1]=3; p.weights[h+8,end,1]=-5; p.weights[h+9,end,1]=3
    mass=ecology_mass(w,p); step!(w,p,.2)
    @test w.fields.deposit[2,2,1]>.19 && p.storage[1]<1.81
    @test med.signal_source[2,2,1]>0 && med.remodeling[2,2,1]>0
    @test ecology_mass(w,p)≈mass rtol=1e-6
    deposited=w.fields.deposit[2,2,1]; stored=p.storage[1]
    p.weights[h+5,end,1]=-10; p.weights[h+6,end,1]=5
    step!(w,p,.1)
    @test w.fields.deposit[2,2,1]<deposited && p.storage[1]>stored
    @test ecology_mass(w,p)≈mass rtol=1e-6
    @test abs(ecology_statistics(p).energy_error)<1e-6
    p.energy.=0; p.uptake_rate=0; before=copy(med.signal_source); step!(w,p,.1)
    @test med.signal_source==before && isempty(p.ids)
end
@testset "NeuralEcology source increments track physical time and paid work" begin
    function emit(dt,steps;work_rate=.1)
        w=World(Grid(4,4;extent=(4,4,1))); med=neural_test_medium(w.grid.size)
        p=NeuralEcology(zeros(Float32,3,1);medium=med,uptake_rate=0,basal_cost=0,
            secretion_rate=0,plasticity_rate=0,speed=0,energy=5,work_rate)
        h=size(p.hidden,1); fill!(p.weights,0)
        p.weights[h+7,end,1]=3; p.weights[h+8,end,1]=-3; p.weights[h+9,end,1]=2
        for _ in 1:steps; step!(w,p,dt); end
        med.signal_source[1,1,1],med.remodeling[1,1,1]
    end
    @test collect(emit(.2,1))≈collect(emit(.1,2))
    @test emit(.2,1;work_rate=0)==(0f0,0f0)
end
@testset "NeuralEcology actual membrane-to-permeability feedback" begin
    w=World(Grid(6,6;extent=(6,6,1))); phase=ImplicitPhaseSeparation(w)
    fill!(phase.phi,0); fill!(phase.interface,0)
    medium=ActiveMedium(w,phase;interface_barrier=0,signal_opening=0)
    p=NeuralEcology(reshape(Float32[2.2,2.2,0],3,1);medium,energy=5,storage=2,
        uptake_rate=0,basal_cost=0,plasticity_rate=0,secretion_rate=1,speed=1,
        signal_rate=0,remodeling_rate=0)
    h=size(p.hidden,1); fill!(p.weights,0)
    p.weights[h+1,end,1]=3; p.weights[h+5,end,1]=5
    p.weights[h+6,end,1]=-10; p.weights[h+8,end,1]=-3
    step!(w,p,.2); medium_observe!(medium)
    resistance=medium.permeability[3,3,1]
    @test resistance<.5
    open=deepcopy(p); open.medium.permeability.=1
    start=p.positions[1,1]; freeworld=deepcopy(w)
    step!(w,p,.1); step!(freeworld,open,.1)
    @test p.positions[1,1]-start<open.positions[1,1]-start
    before=w.fields.deposit[3,3,1]
    p.weights[h+5,end,1]=-10; p.weights[h+6,end,1]=5
    step!(w,p,.2); medium_observe!(medium)
    @test w.fields.deposit[3,3,1]<before
end
@testset "NeuralEcology catabolism spends finite fuel into inert matter" begin
    for fraction in (0f0,.5f0,1f0)
        w=World(Grid(2,2)); p=NeuralEcology(zeros(Float32,3,1);storage=.4,energy=1,
            catabolic_fraction=fraction,energy_yield=3,uptake_rate=1,basal_cost=0,
            work_rate=0,secretion_rate=0,plasticity_rate=0,speed=0)
        h=size(p.hidden,1); fill!(p.weights,0)
        p.weights[h+4,end,1]=4; p.weights[h+8,end,1]=-4
        initial_mass=ecology_mass(w,p); initial_energy=p.energy_input
        nutrient=w.fields.nutrient[1,1,1]; structural=p.storage[1]
        step!(w,p,.2); taken=nutrient-w.fields.nutrient[1,1,1]
        @test p.storage[1]-structural≈taken*(1-fraction) atol=1e-7
        @test p.waste≈taken*fraction atol=1e-7
        @test p.energy_input-initial_energy≈taken*fraction*3 atol=1e-7
        @test ecology_mass(w,p)≈initial_mass rtol=1e-6
        @test abs(ecology_statistics(p).energy_error)<1e-6
        fraction==0 && @test p.energy[1]==1
    end
    @test_throws ArgumentError NeuralEcology(zeros(Float32,3,1);catabolic_fraction=1.1)
    @test_throws ArgumentError NeuralEcology(zeros(Float32,3,1);catabolic_fraction=-.1)
    @test_throws ArgumentError NeuralEcology(zeros(Float32,3,1);energy_yield=-1)
    w=World(Grid(4,4;extent=(4,4,1))); fill!(w.fields.nutrient,.4)
    p=NeuralEcology(Float32[0 1 2 3;0 1 2 3;0 0 0 0];storage=1,energy=.8,
        capacity=16,lifespan=1,uptake_rate=4,basal_cost=.15,work_rate=.1,
        secretion_rate=0,speed=0,plasticity_rate=0,mutation=0)
    h=size(p.hidden,1); fill!(p.weights,0)
    p.weights[h+4,end,:].=4; p.weights[h+8,end,:].=4
    initial_mass=ecology_mass(w,p); initial_energy=p.energy_input
    for _ in 1:500; step!(w,p,.1); end
    @test p.births>0 && p.deaths>4 && isempty(p.ids)
    @test p.waste>0
    @test p.energy_input-initial_energy≈p.energy_yield*p.waste rtol=1e-6
    @test p.energy_input-initial_energy<=p.energy_yield*initial_mass+1e-5
    @test ecology_mass(w,p)≈initial_mass rtol=1e-6
    @test abs(ecology_statistics(p).energy_error)<1e-5
end
@testset "NeuralEcology proportional competition has no feeding-ID priority" begin
    function compete(swapped=false;equal=false,erosion=false)
        w=World(Grid(2,2)); fill!(w.fields.nutrient,0)
        if erosion; w.fields.deposit[1,1,1]=.6
        else; w.fields.nutrient[1,1,1]=.6; end
        p=NeuralEcology(zeros(Float32,3,2);storage=0,energy=1,uptake_rate=erosion ? 0 : 4,
            basal_cost=0,work_rate=erosion ? .1 : 0,secretion_rate=erosion ? 4 : 0,
            speed=0,plasticity_rate=0,catabolic_fraction=.5)
        h=size(p.hidden,1); fill!(p.weights,0); p.weights[h+8,end,:].=-5
        p.weights[h+5,end,:].=-10 # suppress builds while testing contested erosion
        channel=erosion ? h+6 : h+4
        p.weights[channel,end,1]=atanh(-.5f0)
        p.weights[channel,end,2]=equal ? atanh(-.5f0) : atanh(.5f0)
        swapped && (p.ids=reverse(p.ids))
        mass=ecology_mass(w,p); step!(w,p,1)
        @test ecology_mass(w,p)≈mass rtol=1e-6
        @test abs(ecology_statistics(p).energy_error)<1e-6
        @test all(x->x>=0,w.fields.nutrient) && all(x->x>=0,w.fields.deposit)
        if equal
            @test p.storage[1]==p.storage[2]
            return p.storage
        end
        low=findfirst(x->x<0,p.last_actions[erosion ? 6 : 4,:]); high=3-low
        expected=erosion ? Float32[.15,.45] : Float32[.075,.225]
        @test p.storage[[low,high]]≈expected atol=1e-7
        p.storage[[low,high]],p.energy[[low,high]]
    end
    @test compete(false)==compete(true)
    @test compete(false;erosion=true)==compete(true;erosion=true)
    compete(;equal=true); compete(;equal=true,erosion=true)
    # Fuel returned by an earlier ID's death is unavailable until the next step.
    w=World(Grid(2,2)); fill!(w.fields.nutrient,0)
    p=NeuralEcology(zeros(Float32,3,2);uptake_rate=4,storage=0,energy=1,
        basal_cost=0,work_rate=0,secretion_rate=0,plasticity_rate=0,speed=0)
    h=size(p.hidden,1); fill!(p.weights,0); p.weights[h+8,end,:].=-5
    p.storage[1]=2; p.age[1]=p.lifespan
    before=ecology_mass(w,p); step!(w,p,.1)
    @test p.ids==[2] && p.storage==[0f0] && p.energy==[1f0]
    @test w.fields.nutrient[1,1,1]==2
    @test ecology_mass(w,p)≈before rtol=1e-6
    step!(w,p,.1)
    @test p.storage[1]>0 && p.energy[1]>1
end
