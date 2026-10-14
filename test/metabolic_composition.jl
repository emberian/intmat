# Coupled ledgers: audit reservoirs once, regardless of shared object references.
function metabolic_audit_fixture(;trait_mutation=0,capacity=8,lifespan=.8,uptake_rate=2)
    w=World(Grid(3,2,1;extent=(3,2,1)))
    phase=ImplicitPhaseSeparation(w;noise=0,mineral_pinning=0,nutrient_diffusion=0,activity_relaxation=0)
    medium=ActiveMedium(w,phase;nutrient_diffusion=.1,signal_diffusion=0)
    boundary=BoundaryMaintenance(w,medium;work_per_mass=2,decay=.2)
    resources=ResourcePool(w;medium,initial=(.5,.3,.1),diffusion=(.1,.05,.02),
        source_rates=(.2,0,0),source_budgets=(.4,0,0),sink_rates=(0,0,.1),sink_budgets=(0,0,.1))
    p=MetabolicEcology(Float32[0 1;0 1;0 0],resources;boundary,medium,
        storage=2,energy=3,capacity,lifespan,trait_mutation,uptake_rate,
        basal_cost=.2,work_rate=.1,secretion_rate=0,plasticity_rate=0,speed=0,mutation=0,
        birth_material=.2,birth_energy=.1,boundary_rate=.2)
    c=p.controller; fill!(c.weights,0); h=size(c.hidden,1)
    c.weights[h+4,end,:].=3; c.weights[h+5,end,:].=2; c.weights[h+8,end,:].=3
    w,p,resources,boundary,phase
end
metabolic_audit_energy(p)=resource_potential(p.resources)+sum(Float64,p.controller.energy;init=0.)+p.controller.heat+p.chemical_heat
@testset "metabolic composed material and chemical/work reservoirs" begin
    w,p,r,b,phase=metabolic_audit_fixture()
    m0=ecology_mass(w,p); e0=metabolic_audit_energy(p)
    mass_error=0.; energy_error=0.
    for _ in 1:180
        step!(w,(p,b,r,phase),.05)
        mass_error=max(mass_error,abs(ecology_mass(w,p)-m0-sum(r.supplied)+sum(r.removed)))
        rs=resource_statistics(r)
        energy_error=max(energy_error,abs(metabolic_audit_energy(p)-e0-rs.supplied_potential+rs.removed_potential))
    end
    @test mass_error<5e-5
    @test energy_error<5e-5
    @test p.controller.births>0 && p.controller.deaths>0
    @test isempty(p.controller.ids)
    @test b.contributed_material>0 && b.spent_work>0 && b.recycled_material>0
    @test all(a->minimum(a)>=0,r.fields)
    @test abs(ecology_statistics(p).energy_error)<5e-5
    before=p.controller.energy_input
    for _ in 1:20; step!(w,(p,b,r,phase),.05); end
    @test p.controller.energy_input==before
end
@testset "metabolic inherited traits survive reorder birth and death" begin
    w,p,r,b,phase=metabolic_audit_fixture(;lifespan=10,uptake_rate=0)
    p.allocation[:,1].=Float32[1,0,0]; p.allocation[:,2].=Float32[0,1,0]
    p.assimilation.=Float32[.2,.8]; p.efficiency.=Float32[.3,.7]; p.cooperation.=Float32[.1,.9]; p.byproduct_fraction.=Float32[.1,.8]
    v=deepcopy(w); q=deepcopy(p); reorder_ecology!(q,[2,1])
    step!(w,p,.05); step!(v,q,.05)
    @test p.controller.ids==q.controller.ids
    @test p.allocation==q.allocation && p.efficiency==q.efficiency && p.cooperation==q.cooperation
    @test p.controller.births==2
    for i in 3:4
        parent=findfirst(==(p.controller.parents[i]),p.controller.ids)
        @test p.allocation[:,i]==p.allocation[:,parent]
        @test p.assimilation[i]==p.assimilation[parent]
        @test p.efficiency[i]==p.efficiency[parent]
        @test p.cooperation[i]==p.cooperation[parent]
        @test p.byproduct_fraction[i]==p.byproduct_fraction[parent]
        @test all(iszero,p.controller.hidden[:,i])
    end
    p.controller.age[1]=p.controller.lifespan
    step!(w,p,.05)
    @test !(1 in p.controller.ids)
    @test size(p.allocation,2)==length(p.controller.ids)==length(p.cooperation)
end
@testset "downhill crossfeeding cannot regenerate spent potential" begin
    w=World(Grid(1,1,1)); r=ResourcePool(w;initial=(1,0,0),potentials=(4,2,0),diffusion=(0,0,0))
    cell=CartesianIndex(1,1,1); e0=resource_potential(r); m0=resource_mass(r)
    first=resource_convert!(r,1,2,cell,1;efficiency=.5)
    second=resource_convert!(r,2,3,cell,1;efficiency=.5)
    @test resource_mass(r)==m0
    @test resource_potential(r)+first.energy+first.heat+second.energy+second.heat==e0
    @test_throws ArgumentError resource_convert!(r,3,1,cell,1)
    @test resource_convert!(r,2,3,cell,1).energy==0
end
@testset "metabolic checkpoint preserves coupled aliases and recordings" begin
    w,p,r,b,phase=metabolic_audit_fixture(;lifespan=10)
    a=Arena(w,(p,b,r,phase);dt=.05)
    arena_step!(a)
    mktempdir() do dir
        record=joinpath(dir,"capture"); record_state(w,a.processes,record;step=1)
        manifest=Interior.JSON3.read(read(joinpath(record,"state.json"),String))
        descriptor=manifest.processes[1].state.allocation
        captured=copy(reinterpret(Float32,read(joinpath(record,String(descriptor.file)))))
        @test captured==vec(p.allocation)
        for (index,key,expected) in ((2,:film,b.film),(3,:field_1,r.fields[1]))
            d=getproperty(manifest.processes[index].state,key)
            @test reinterpret(Float32,read(joinpath(record,String(d.file))))==vec(expected)
        end
        path=joinpath(dir,"restart.jls"); checkpoint!(a,path); q=restore_arena(path)
        qp,qb,qr,qphase=q.processes
        @test qp.resources===qr && qp.boundary===qb
        @test qp.controller.medium===qb.medium===qr.medium
        @test qb.medium.phase===qphase
        @test qb.medium.membrane===q.world.fields.deposit
        for _ in 1:8; arena_step!(a); arena_step!(q); end
        @test a.world.fields==q.world.fields
        @test p.allocation==qp.allocation && p.controller.weights==qp.controller.weights
        @test r.fields==qr.fields && b.film==qb.film
        @test reinterpret(Float32,read(joinpath(record,String(descriptor.file))))==captured
    end
end
@testset "shared film is funded once from caller reservoirs" begin
    w,p,r,b,phase=metabolic_audit_fixture(;uptake_rate=0,lifespan=10,capacity=2)
    c=p.controller; c.basal_cost=0; c.work_rate=0; p.boundary_rate=1
    fill!(c.weights,0); h=size(c.hidden,1); c.weights[h+5,end,:].=3; c.weights[h+8,end,:].=-3
    material=sum(Float64,c.storage); work=sum(Float64,c.energy); heat=c.heat
    step!(w,p,.1)
    @test material-sum(Float64,c.storage)≈b.contributed_material atol=2e-7
    @test work-sum(Float64,c.energy)≈b.spent_work atol=2e-7
    @test c.heat-heat≈b.spent_work atol=2e-7
    @test b.spent_work≈2b.contributed_material atol=2e-7
    @test boundary_material_mass(b)≈b.contributed_material atol=2e-7
    @test r.withdrawn==zeros(3)
    @test Interior.arena_inventory(w,(p,b,r,phase))≈ecology_mass(w,p) atol=2e-7
    p.cooperation.=0
    before=boundary_material_mass(b); step!(w,p,.1)
    @test boundary_material_mass(b)==before
end
@testset "owned stages cannot be scheduled twice" begin
    w,p,r,b,phase=metabolic_audit_fixture()
    baseline=deepcopy(w.fields); clock=w.time; film=copy(b.film)
    for stages in ((b,b.medium),(p,p.controller),(p,p),(r,r),(b,b))
        @test_throws ArgumentError step!(w,stages,.05)
        @test w.fields==baseline && w.time==clock && b.film==film
    end
    # Borrowed resources and boundary remain explicit independent stages.
    @test isnothing(Interior.validate_composition((p,b,r,phase)))
    # Immutable parameterized operators may legitimately repeat.
    @test isnothing(Interior.validate_composition((Deposition(),Deposition())))
end
