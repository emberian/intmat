# These fixtures silence motion and background turnover only when isolating a
# chemical or hereditary mechanism. Integration checks below use the real step.
function cellular_fixture(;positions=reshape([2.,2.,0.],3,1),initial=(.05,.05,0.),kw...)
    w=World(Grid(16,16;extent=(8,8,1));seed=37)
    resources=ResourcePool(w;potentials=(4,3,0),initial,diffusion=(0,0,0))
    defaults=(;capacity=8,speed=0.,repulsion=0.,stiffness=0.,basal_rate=0.,neural_cost=0.,
        plasticity_rate=0.,mutation=0.,developmental_partition=0.,division_mass=10.)
    p=Interior.CellularEcology(w,positions,resources;merge(defaults,(;kw...))...)
    fill!(p.weights,0)
    p.weights[size(p.hidden,1)+10,end,:].=-2 # Explicitly gate division off.
    w,p
end

cellular_test_mass(p)=resource_mass(p.resources)+(sum(p.biomass)+sum(p.metabolites)+p.inert_residual)
cellular_test_energy(p)=resource_potential(p.resources)+sum(p.energy)+p.work_spent+p.chemical_heat+
    sum(p.metabolites[s,i]*p.resources.potentials[s] for s in 1:2,i in eachindex(p.ids);init=0.)

function cellular_test_same_state(p,q)
    for name in (:ids,:parents,:lineages,:generations,:positions,:biomass,:metabolites,:energy,
            :hidden,:adaptation,:weights,:plastic_marks,:last_actions,:age,:next_id,:births,:deaths,
            :blocked_births,:work_spent,:chemical_heat,:inert_residual,:consumed,:junction_transferred,
            :synthesized,:oxidized,:shed)
        @test getproperty(p,name)==getproperty(q,name)
    end
    @test p.bonds.endpoints==q.bonds.endpoints
    @test p.bonds.rest_lengths==q.bonds.rest_lengths
    @test p.resources.fields==q.resources.fields
    @test p.resources.withdrawn==q.resources.withdrawn
    @test p.resources.returned==q.resources.returned
end

@testset "Cellular complementary physical niche geometry" begin
    cfg=Dict{String,Any}("model"=>"cellular","grid"=>Dict("shape"=>[16,12,1],"extent"=>[8.,6.,1.]),
        "cellular"=>Dict("population"=>0),"resources"=>Dict{String,Any}("niche_waves"=>[2,3,0],"niche_contrast"=>.9))
    a=build_arena(cfg); r=a.processes[1]
    @test r.source_masks[1].+r.source_masks[2]≈fill(2.,16,12,1)
    @test sum(Float64,r.source_masks[1])/192≈1 atol=1e-7
    @test r.source_masks[1][1,1,1]!=r.source_masks[1][1,3,1]
    @test minimum(r.source_masks[1])>=.1 && maximum(r.source_masks[1])<=1.9
    delete!(cfg["resources"],"niche_waves"); a=build_arena(cfg)
    @test a.processes[1].source_masks[1][:,1,:]==a.processes[1].source_masks[1][:,7,:]
    cfg["resources"]["niche_waves"]=[1,0,1]
    @test_throws ArgumentError build_arena(cfg)
end

@testset "Cellular terminal bodies and borrowed energy ownership" begin
    w,p=cellular_fixture(;metabolites=0.,energy=0.,uptake_rate=0.,growth_rate=0.,catabolic_rate=0.,
        basal_rate=1.,starvation_rate=10.)
    mass=cellular_test_mass(p); energy=cellular_test_energy(p)
    step!(w,p,.2)
    @test isempty(p.ids) && p.deaths==1
    @test cellular_test_mass(p)≈mass atol=1e-12
    @test cellular_test_energy(p)≈energy atol=1e-12
    w,p=cellular_fixture()
    @test Interior.arena_energy((p,))==Interior.arena_energy((p.resources,p))==Interior.arena_energy((p,p.resources))
    @test Interior.arena_energy((p,))≈cellular_test_energy(p)
end

@testset "Cellular immediate spatial permutation observations" begin
    w,p=cellular_fixture(;positions=[1. 6.;1. 6.;0. 0.])
    p.biomass.=[.5,1.5]
    p.resources.fields[1][1:8,:,:].=1
    before=sample_concentrations(p.resources,p.stencil)
    reorder_ecology!(p,[2,1])
    @test sample_concentrations(p.resources,p.stencil)≈before[:,[2,1]]
end

@testset "Cellular uptake, complementary growth, and chemical accounting" begin
    w,p=cellular_fixture(;initial=(1.,1.,0.),metabolites=0.,growth_rate=0.,catabolic_rate=0.,uptake_rate=100.)
    mass=cellular_test_mass(p); chemical=cellular_test_energy(p)
    dissolved=[sum(Float64,a) for a in p.resources.fields]
    for _ in 1:10; step!(w,p,.1); end
    @test all(p.metabolites.>0)
    # Per-site Float32 withdrawals can leave a sub-ulp amount of unused room.
    @test sum(p.metabolites)≈p.parameters.store_capacity*p.biomass[1] atol=1e-6
    @test all(vec(sum(p.metabolites;dims=1)).<=p.parameters.store_capacity.*p.biomass.+1e-14)
    @test p.consumed≈p.resources.withdrawn[1:2] atol=1e-14
    @test vec(p.metabolites)≈dissolved[1:2].-[sum(Float64,a) for a in p.resources.fields[1:2]] atol=1e-13
    @test cellular_test_mass(p)≈mass atol=1e-12
    @test cellular_test_energy(p)≈chemical atol=1e-12
    @test p.energy==[1.] && p.oxidized==0 && p.synthesized==0

    # One precursor alone cannot construct biomass, regardless of available work.
    w,p=cellular_fixture(;uptake_rate=0.,catabolic_rate=0.,metabolites=0.,growth_rate=2.)
    p.metabolites[:,1].=[.3,0]
    mass=cellular_test_mass(p); chemical=cellular_test_energy(p)
    step!(w,p,.2)
    @test p.biomass==[1.] && p.synthesized==0
    @test p.metabolites[:,1]==[.3,0]
    @test cellular_test_mass(p)≈mass atol=1e-12
    @test cellular_test_energy(p)≈chemical atol=1e-12
    # Supply the complement internally and observe stoichiometry and paid work.
    p.metabolites[2,1]=.2
    before=copy(p.metabolites); battery=p.energy[1]
    mass=cellular_test_mass(p); chemical=cellular_test_energy(p)
    step!(w,p,.2)
    built=p.biomass[1]-1
    @test built>0
    @test before.-p.metabolites≈fill(built/2,2,1) atol=1e-15
    @test battery-p.energy[1]≈built*p.parameters.growth_work atol=1e-15
    @test p.work_spent≈built*p.parameters.growth_work atol=1e-15
    @test p.chemical_heat≈built/2*7 atol=1e-15
    @test p.synthesized≈built atol=1e-15
    @test cellular_test_mass(p)≈mass atol=1e-12
    @test cellular_test_energy(p)≈chemical atol=1e-12
end

@testset "Cellular oxidation, finite battery, starvation, and death" begin
    w,p=cellular_fixture(;uptake_rate=0.,growth_rate=0.,metabolites=.3,energy=.01,
        energy_capacity=.1,catabolic_rate=100.)
    mass=cellular_test_mass(p); chemical=cellular_test_energy(p)
    step!(w,p,.2)
    @test p.oxidized>0
    @test sum(Float64,p.resources.fields[3])≈p.oxidized atol=1e-14
    @test p.resources.returned[3]≈p.oxidized atol=1e-14
    @test sum(p.metabolites)+p.oxidized≈.6 atol=1e-14
    @test 0<=p.energy[1]<=p.parameters.energy_capacity*p.biomass[1]
    @test p.energy[1]>.01 && p.chemical_heat>0
    @test cellular_test_mass(p)≈mass atol=1e-12
    @test cellular_test_energy(p)≈chemical atol=1e-12

    # Planned synthesis can overestimate work when one precursor is absent.
    # Excess respiratory yield must become heat rather than an unbounded battery.
    w,p=cellular_fixture(;uptake_rate=0.,growth_rate=100.,catabolic_rate=100.,
        metabolites=.3,energy=.01,energy_capacity=.1)
    p.metabolites[:,1].=[.6,0.]
    mass=cellular_test_mass(p); chemical=cellular_test_energy(p)
    step!(w,p,.2)
    @test p.energy[1]≈p.parameters.energy_capacity*p.biomass[1] atol=1e-14
    @test p.synthesized==0 && p.oxidized>.5
    @test p.chemical_heat>(1-p.parameters.efficiency)*4p.oxidized
    @test cellular_test_mass(p)≈mass atol=1e-12
    @test cellular_test_energy(p)≈chemical atol=1e-12

    # Shrinking structure must not leave intracellular stores over capacity.
    w,p=cellular_fixture(;uptake_rate=0.,growth_rate=0.,catabolic_rate=0.,
        metabolites=.39,energy=0.,basal_rate=1.,starvation_rate=1.)
    mass=cellular_test_mass(p); chemical=cellular_test_energy(p)
    step!(w,p,.2)
    @test p.shed>0 && p.biomass[1]<1
    @test sum(p.metabolites)<=p.parameters.store_capacity*p.biomass[1]+1e-14
    @test cellular_test_mass(p)≈mass atol=1e-12
    @test cellular_test_energy(p)≈chemical atol=1e-12

    # Force a death in a normal field, then with deposits below the field's ulp.
    # The latter must leave an explicit inert remainder, not an immortal tiny cell.
    for waste in (0.,1e10)
        w,p=cellular_fixture(;initial=(0.,0.,waste),biomass=.26,energy=0.,metabolites=.02,
            uptake_rate=0.,growth_rate=0.,catabolic_rate=0.,basal_rate=1.,starvation_rate=1.)
        # In the high-background case, start below the death threshold because
        # the field cannot represent even the preceding structural damage.
        waste>0 && (p.biomass[1]=.24)
        owned=sum(p.biomass)+sum(p.metabolites)
        mass=cellular_test_mass(p); chemical=cellular_test_energy(p)
        before_returned=p.resources.returned[3]
        step!(w,p,.2)
        @test isempty(p.ids) && p.deaths==1
        @test p.resources.returned[3]-before_returned+p.inert_residual≈owned atol=1e-14
        @test cellular_test_mass(p)≈mass atol=1e-12
        @test cellular_test_energy(p)≈chemical atol=1e-12
        waste>0 && @test p.inert_residual>0
        @test size(p.hidden,2)==size(p.metabolites,2)==size(p.weights,3)==0
        step!(w,p,.2)
        @test p.deaths==1 && isempty(p.ids)
        @test cellular_test_energy(p)≈chemical atol=1e-12
    end
end

@testset "Cellular developmental inheritance and sparse mutation" begin
    w,p=cellular_fixture(;biomass=2.,energy=1.6,metabolites=.2,division_mass=1.8)
    p.hidden[:,1].=range(-.4,.4;length=size(p.hidden,1))
    p.adaptation[:,1].=.13f0; p.plastic_marks[:,:,1].=.07f0
    genotype=copy(p.weights[:,:,1]); hidden=copy(p.hidden[:,1])
    adaptation=copy(p.adaptation[:,1]); marks=copy(p.plastic_marks[:,:,1])
    mass=cellular_test_mass(p); chemical=cellular_test_energy(p)
    Interior._cellular_divide!(w,p,[1])
    @test p.ids==[1,2] && p.parents==[0,1] && p.lineages==[1,1] && p.generations==[0,1]
    @test p.next_id==3 && p.births==1
    @test p.biomass==[1.,1.] && p.energy==[.8,.8]
    @test p.metabolites==fill(.1,2,2)
    @test p.weights[:,:,1]==p.weights[:,:,2]==genotype
    @test p.hidden[:,1]==p.hidden[:,2]==hidden
    @test p.adaptation[:,1]==p.adaptation[:,2]==adaptation
    @test p.plastic_marks[:,:,1]==p.plastic_marks[:,:,2]==marks
    @test p.bonds.endpoints==[(1,2)]
    @test p.positions[:,1]!=p.positions[:,2]
    @test all(0 .<= p.positions .< reshape(collect(w.grid.extent),3,1))
    @test cellular_test_mass(p)≈mass atol=1e-12
    @test cellular_test_energy(p)≈chemical atol=1e-12

    w,p=cellular_fixture(;biomass=2.,division_mass=1.8,mutation=.2,mutation_probability=.04)
    genotype=copy(p.weights[:,:,1]); q=deepcopy(p)
    Interior._cellular_divide!(w,p,[1]); Interior._cellular_divide!(w,q,[1])
    changed=p.weights[:,:,2].!=genotype
    @test 0<count(changed)<length(changed)÷4
    @test p.weights[:,:,1]==genotype
    @test p.weights[:,:,2][.!changed]==genotype[.!changed]
    @test p.weights==q.weights && p.positions==q.positions

    # The public lifecycle updates development once, then transmits that state.
    w,p=cellular_fixture(;biomass=2.,division_mass=1.8,uptake_rate=0.,growth_rate=0.,catabolic_rate=0.,
        plasticity_rate=.01)
    h=size(p.hidden,1); p.hidden[:,1].=.2; p.weights[1:h,end,1].=.4
    p.weights[h+10,end,1]=2
    mass=cellular_test_mass(p); chemical=cellular_test_energy(p)
    step!(w,p,.1)
    @test p.births==1 && length(p.ids)==2
    @test p.hidden[:,1]==p.hidden[:,2] && any(!iszero,p.hidden)
    @test p.plastic_marks[:,:,1]==p.plastic_marks[:,:,2] && any(!iszero,p.plastic_marks)
    @test p.adaptation[:,1]==p.adaptation[:,2] && any(!iszero,p.adaptation)
    @test cellular_test_mass(p)≈mass atol=1e-12
    @test cellular_test_energy(p)≈chemical atol=1e-12
end

@testset "Cellular stable IDs, permutation parity, and explicit capacity" begin
    positions=[2. 5.;2. 5.;0. 0.]
    w,p=cellular_fixture(;positions,biomass=2.,division_mass=1.8,capacity=4,
        mutation=.2,mutation_probability=.04)
    h=size(p.hidden,1); p.weights[h+10,end,:].=2
    p.hidden[:,1].=.1; p.hidden[:,2].=-.2
    v,q=deepcopy((w,p)); reorder_ecology!(q,[2,1])
    @test q.ids==[2,1] && q.hidden[:,1]==p.hidden[:,2]
    @test_throws ArgumentError reorder_ecology!(q,[1,1])
    for _ in 1:3; step!(w,p,.1); step!(v,q,.1); end
    @test p.births==2 && p.ids==[1,2,3,4] && p.parents==[0,0,1,2]
    cellular_test_same_state(p,q)
    @test w.time==v.time
    w,p=cellular_fixture(;positions,biomass=2.,division_mass=1.8,capacity=2)
    p.weights[size(p.hidden,1)+10,end,:].=2
    mass=cellular_test_mass(p)
    step!(w,p,.1)
    @test p.births==0 && p.blocked_births==2 && length(p.ids)==2
    @test cellular_test_mass(p)≈mass atol=1e-12
end

@testset "Cellular arena ledgers and exact checkpoint continuation" begin
    w,p=cellular_fixture(;biomass=2.,division_mass=1.8,capacity=5,mutation=.1,mutation_probability=.03)
    h=size(p.hidden,1); p.weights[h+10,end,:].=2
    p.medium=ActiveMedium(w,ImplicitPhaseSeparation(w))
    a=Arena(w,(p.resources,p);dt=.1)
    chemical=cellular_test_energy(p)
    @test Interior.arena_energy(a.processes)≈chemical atol=1e-12
    @test Interior.arena_inventory(w,a.processes)≈
        sum(Float64,w.fields.nutrient)+sum(Float64,w.fields.deposit)+cellular_test_mass(p) atol=1e-12
    @test Interior.owned_material_mass(p)≈sum(p.biomass)+sum(p.metabolites)+p.inert_residual atol=1e-15
    for _ in 1:3; arena_step!(a); end
    @test Interior.arena_inventory(w,a.processes)≈a.initial_mass atol=1e-11
    @test Interior.arena_energy(a.processes)≈chemical atol=1e-11
    mktempdir() do directory
        file=joinpath(directory,"cellular.jls")
        checkpoint!(a,file); b=restore_arena(file)
        resources,cells=b.processes
        @test cells.resources===resources
        @test resources!==p.resources && cells!==p
        @test cells.stencil.grid===b.world.grid
        @test cells.medium.membrane===b.world.fields.deposit
        @test cells.resources.medium===nothing
        for _ in 1:5; arena_step!(a); arena_step!(b); end
        cellular_test_same_state(p,cells)
        @test a.steps==b.steps && a.world.time==b.world.time
        @test a.world.fields==b.world.fields
        @test Interior.arena_inventory(b.world,b.processes)≈b.initial_mass atol=1e-11
        @test Interior.arena_energy(b.processes)≈chemical atol=1e-11
    end
end

@testset "Contiguous lifecycle growth preserves borrowed old arrays" begin
    w,p=cellular_fixture(positions=[2. 5.;2. 5.;0. 0.],biomass=2.,
        mutation_probability=0.,developmental_partition=.03)
    rng=MersenneTwister(842)
    p.weights.=randn(rng,Float32,size(p.weights)).*.1f0
    p.hidden.=rand(rng,Float32,size(p.hidden)).*.2f0
    p.adaptation.=rand(rng,Float32,size(p.adaptation)).*.1f0
    old=(positions=p.positions,metabolites=p.metabolites,hidden=p.hidden,
        adaptation=p.adaptation,weights=p.weights,plastic_marks=p.plastic_marks,last_actions=p.last_actions)
    saved=deepcopy(old)
    Interior._cellular_divide!(w,p,[1,2])
    @test length(p.ids)==4 && p.parents[3:4]==[1,2]
    for name in keys(old)
        @test isequal(getproperty(old,name),getproperty(saved,name))
        @test getproperty(p,name)!==getproperty(old,name)
    end
    @test p.weights[:,:,1]==p.weights[:,:,3]
    @test p.weights[:,:,2]==p.weights[:,:,4]
    parent=p.weights[1,1,1]
    p.weights[1,1,3]+=1f0
    @test p.weights[1,1,1]==parent
    @test all(isfinite,p.positions) && all(isfinite,p.hidden) && all(isfinite,p.adaptation)
    @test size(p.weights,3)==size(p.plastic_marks,3)==length(p.ids)
end
