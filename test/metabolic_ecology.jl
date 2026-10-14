@testset "Metabolic ecology finite fuel, specialization and inherited traits" begin
    function metabolic_fixture(; initial=(1,0), allocation=Float32[1 0;0 1],kwargs...)
        w=World(Grid(3,3)); r=ResourcePool(w;names=(:a,:b),potentials=(4,4),initial,diffusion=(0,0))
        p=MetabolicEcology(zeros(Float32,3,2),r;allocation,storage=0,energy=1,basal_cost=0,work_rate=0,secretion_rate=0,speed=0,plasticity_rate=0,kwargs...)
        c=p.controller; h=size(c.hidden,1); fill!(c.weights,0); c.weights[h+4,end,:].=4; c.weights[h+8,end,:].=-4
        w,p
    end
    w,p=metabolic_fixture(); mass=ecology_mass(w,p); potential=resource_potential(p.resources)
    q=deepcopy(p); v=deepcopy(w); reorder_ecology!(q,[2,1])
    ecology_step!(w,p,.2); ecology_step!(v,q,.2)
    @test p.controller.storage[1]>p.controller.storage[2]
    @test p.controller.storage==q.controller.storage && p.allocation==q.allocation
    @test p.resources.fields==q.resources.fields
    @test ecology_mass(w,p)≈mass atol=1e-6
    released=potential-resource_potential(p.resources)
    @test released≈p.chemical_heat+p.controller.energy_input-2 atol=1e-6
    @test abs(ecology_statistics(p).energy_error)<1e-6
    w2,p2=metabolic_fixture(initial=(0,1)); ecology_step!(w2,p2,.2)
    @test p2.controller.storage[2]>p2.controller.storage[1]
    # Structural corpses and deposits cannot recharge exhausted chemical fuel.
    for field in p.resources.fields; fill!(field,0); end
    fill!(w.fields.nutrient,100); fill!(w.fields.deposit,100)
    input=p.controller.energy_input; ecology_step!(w,p,.2)
    @test p.controller.energy_input==input
    # Zero mutation inherits genotype exactly; recurrent and plastic memory reset.
    wb,pb=metabolic_fixture(trait_mutation=0,mutation=0,capacity=4)
    cb=pb.controller; cb.storage.=2; cb.energy.=5; cb.plastic_marks.=.1
    h=size(cb.hidden,1); cb.weights[h+8,end,:].=4
    traits=copy(pb.allocation); ecology_step!(wb,pb,.1)
    @test cb.births==2 && pb.allocation[:,3:4]==traits
    @test pb.efficiency[3:4]==pb.efficiency[1:2]
    @test all(iszero,cb.hidden[:,3:4]) && all(iszero,cb.plastic_marks[:,:,3:4])
    wm,pm=metabolic_fixture(trait_mutation=.1,capacity=4)
    cm=pm.controller; cm.storage.=2; cm.energy.=5; hm=size(cm.hidden,1); cm.weights[hm+8,end,:].=4
    ecology_step!(wm,pm,.1)
    @test pm.allocation[:,3:4]!=pm.allocation[:,1:2]
    @test all(x->0<=x<=1,pm.allocation) && all(isapprox.(vec(sum(pm.allocation;dims=1)),1;atol=1e-6))
    # Efficiency costs throughput, so it is not a free gain in every trait.
    wa,pa=metabolic_fixture(); we,pe=metabolic_fixture()
    pa.efficiency.=0; pe.efficiency.=1
    ecology_step!(wa,pa,.2); ecology_step!(we,pe,.2)
    @test pa.consumed[1]>pe.consumed[1]
    @test pe.controller.energy[1]>pa.controller.energy[1]
end
@testset "Crossfeeding retains finite chemical potential" begin
    w=World(Grid(2,2)); r=ResourcePool(w;initial=(1,0,0),potentials=(4,3,1),diffusion=(0,0,0))
    p=MetabolicEcology(zeros(Float32,3,2),r;allocation=Float32[1 0;0 0;0 1],assimilation=0,efficiency=1,byproduct_fraction=1,storage=0,energy=1,basal_cost=0,work_rate=0,secretion_rate=0,speed=0,plasticity_rate=0)
    c=p.controller; h=size(c.hidden,1); fill!(c.weights,0); c.weights[h+4,end,:].=4; c.weights[h+8,end,:].=-4
    mass=ecology_mass(w,p); potential=resource_potential(r); initial=c.energy_input
    ecology_step!(w,p,.2)
    @test r.fields[3][1,1,1]>0 && p.consumed[3]==0 # same-step exports unavailable
    ecology_step!(w,p,.2)
    @test p.consumed[3]>0 && c.energy[2]>1
    for _ in 1:200; ecology_step!(w,p,.2); end
    @test ecology_mass(w,p)≈mass atol=2e-6
    @test resource_potential(r)+p.chemical_heat+c.energy_input-initial≈potential atol=2e-6
    @test c.energy_input-initial<=potential+1e-6
    @test c.waste>=0 && abs(ecology_statistics(p).energy_error)<2e-5
end
