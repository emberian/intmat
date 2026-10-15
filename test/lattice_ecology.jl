using Test,Interior
include(joinpath(@__DIR__,"..","studies","lattice_ecology.jl"))
using .LatticeEcologyStudy
const LE=LatticeEcologyStudy

function pair_fixture(;energy=1f0,matrix=0f0,mutation=1f0)
    w=lattice_world(Grid(2,1,1;extent=(2,1,1));seed=11)
    q=lattice_parameters(move_rate=1,division_rate=1,mutation_probability=mutation)
    p=LatticeEcology(w;parameters=q,loci=4,founder_genotypes=1,occupancy=0)
    p.genotype[1]=1;p.lineage[1]=1;p.ids[1]=1;p.cells[1,:].=(1,energy,7)
    p.hidden[1,:].=.4;p.environment[2,4]=matrix
    w,p
end
function contact!(w,p;tick=1)
    LE.matching_kernel!(CPU(),128)(p.cells,p.hidden,p.actions,p.environment,p.genotype,p.lineage,p.generation,p.ids,p.parents,p.events,p.heat,p.work,p.dims,p.parameters,1,1,p.seed,tick,p.loci,1f0;ndrange=2)
end
@testset "funded lattice inheritance and physical occupancy" begin
    w,p=pair_fixture();p.actions[1,3]=1
    before=lattice_statistics(w,p);weights=copy(p.weights);contact!(w,p);after=lattice_statistics(w,p)
    @test sum(p.cells[:,1])==1f0
    @test sum(p.cells[:,2])+sum(p.heat)≈1f0 atol=2e-7
    @test p.parents[2]==p.ids[1] && p.ids[2]!=p.ids[1]
    @test p.generation==[0,1] && p.lineage==[1,1]
    @test p.hidden[1,:]==fill(.4f0,4) && all(iszero,p.hidden[2,:])
    @test count_ones(UInt32(p.genotype[1]-1)⊻UInt32(p.genotype[2]-1))==1
    @test p.events[2,1]==1 && p.events[2,4]==1 && p.weights==weights
    @test after.material≈before.material atol=2e-6
    for (energy,matrix) in ((.01f0,0f0),(1f0,1f0))
        w,p=pair_fixture(;energy,matrix);p.actions[1,3]=1;p.actions[1,2]=0
        before=copy(p.cells);contact!(w,p)
        @test p.cells==before && p.genotype[2]==0 && sum(p.events[:,1])==0
        @test sum(p.heat)==0
    end
    w,p=pair_fixture(mutation=0);p.actions[1,3]=1;contact!(w,p)
    @test p.genotype==[1,1]
end
@testset "paid hopping carries identity and acquired state" begin
    w,p=pair_fixture();p.actions[1,2]=1
    contact!(w,p)
    @test p.genotype==[0,1] && p.ids==[0,1] && p.lineage==[0,1]
    @test p.cells[:,1]==[0,1] && p.cells[2,3]==7
    @test p.hidden[2,:]==fill(.4f0,4) && all(iszero,p.hidden[1,:])
    @test p.cells[2,2]+sum(p.heat)≈1f0 atol=2e-7
    @test sum(p.events[:,3])==1 && sum(p.events[:,1])==0
end
@testset "starvation cannot fund growth, work or descendants" begin
    w=lattice_world(Grid(1,1,1;extent=(1,1,1)))
    p=LatticeEcology(w;loci=4,founder_genotypes=1,occupancy=1,precursor=1,fuel=0,reserve=0)
    p.cells[1,2]=0;before=lattice_statistics(w,p);step!(w,p,.1);after=lattice_statistics(w,p)
    @test after.population==0 && after.births==0 && after.deaths==1
    @test after.matrix==0 && after.energy_accounted==0 && sum(after.work)==0
    @test after.material≈before.material atol=1e-7
end
@testset "dense composition, ledgers, identity and exact continuation" begin
    for dims in ((8,8,1),(4,4,4))
        w=lattice_world(Grid(dims...;extent=Float32.(dims));seed=29)
        p=LatticeEcology(w;loci=6,occupancy=.25,founder_genotypes=4)
        initial=lattice_statistics(w,p)
        for _ in 1:1000;step!(w,p,.1);end
        final=lattice_statistics(w,p)
        @test final.finite && final.minimum_physical>=0 && final.hidden_maxabs<=1.000001
        @test final.capacity_excess<=2e-6
        @test final.population==initial.population+final.births-final.deaths
        @test final.births>0 && final.deaths>0 && final.hops>0 && final.mutations>0
        @test abs(final.material-initial.material)<1e-4initial.material
        @test abs(final.energy_accounted-initial.energy_accounted)<1e-4initial.energy_accounted
        live=findall(>(0),p.genotype)
        @test length(unique(p.ids[live]))==length(live) && all(>(0),p.ids[live])
        @test all(p.cells[findall(==(0),p.genotype),1].==0)
        @test Interior.experiment_diagnostics(w,(p,)).material_mass≈final.material
        prior=copy(p.cells);tick=p.tick
        @test_throws ArgumentError step!(w,p,10.)
        @test p.cells==prior && p.tick==tick
        @test_throws ArgumentError step!(w,(p,p),.1)
        mktempdir() do dir
            path=joinpath(dir,"state.jls");lattice_checkpoint(w,p,path);v,q=restore_lattice(path)
            @test v.time==w.time && q.tick==p.tick && q.weights==p.weights
            for _ in 1:13;step!(w,p,.1);step!(v,q,.1);end
            for name in (:environment,:cells,:hidden,:genotype,:ids,:parents,:heat,:events)
                @test getfield(p,name)==getfield(q,name)
            end
            manifest=record_state(w,p,joinpath(dir,"native");step=p.tick)
            @test manifest["processes"][1]["state"]["genotype"]["shape"]==[prod(dims)]
        end
    end
    @test_throws ArgumentError LatticeEcology(lattice_world(Grid(7,8,1)))
    @test_throws ArgumentError lattice_parameters(efficiency=2)
    @test_throws ArgumentError lattice_parameters(division_mass=.2)
end
@testset "declared injury transfers body/matrix/work rather than deleting it" begin
    w=lattice_world(Grid(8,8,1;extent=(8,8,1)));p=LatticeEcology(w;loci=4,occupancy=.5,founder_genotypes=2)
    for _ in 1:100;step!(w,p,.1);end
    before=lattice_statistics(w,p);lattice_damage!(w,p;radius=2.);after=lattice_statistics(w,p)
    @test after.population<before.population && after.matrix<before.matrix
    @test abs(after.material-before.material)<1e-5before.material
    @test abs(after.energy_accounted-before.energy_accounted)<1e-5before.energy_accounted
    @test after.deaths-before.deaths==before.population-after.population
end
@testset "invalid representable scales and owner aliases fail before a step" begin
    @test_throws ArgumentError lattice_parameters(growth_work=1e-100)
    @test_throws ArgumentError lattice_parameters(lifespan=1e-100)
    @test_throws ArgumentError lattice_genotypes(loci=4,strength=Inf)
    w=lattice_world(Grid(4,4,1;extent=(4,4,1)));p=LatticeEcology(w;loci=4)
    p.actions=p.cells;before=copy(p.cells)
    @test_throws ArgumentError step!(w,p,.1)
    @test p.cells==before && p.tick==0 && w.time==0
end
@testset "funding and template regressions at representational boundaries" begin
    @test_throws ArgumentError LatticeEcology(lattice_world(Grid(2,2,1));reserve=1e100)
    @test_throws ArgumentError LatticeEcology(lattice_world(Grid(2,2,1));reserve=2e38,niches=true)
    w=lattice_world(Grid(1,1,1));p=LatticeEcology(w;loci=1,founder_genotypes=1,occupancy=0,fuel=0,reserve=8)
    step!(w,p,1f-8)
    @test all(iszero,p.environment[:,1:2]) && p.environment[:,6:7]==fill(8f0,1,2)
    q=lattice_parameters(energy_capacity=.1)
    p=LatticeEcology(w;parameters=q,loci=1,founder_genotypes=1,occupancy=1)
    @test p.cells[1,2]<=q.energy_capacity*p.cells[1,1]
    q=lattice_parameters(growth=100,uptake=0,feed=0,basal_work=0,neural_work=0,construction=0,erosion=0)
    p=LatticeEcology(w;parameters=q,loci=1,founder_genotypes=1,occupancy=1,fuel=0,reserve=0,precursor=1)
    fill!(p.weights,0);p.cells[1,2]=.0063005444f0
    before=lattice_statistics(w,p);step!(w,p,.1);after=lattice_statistics(w,p)
    @test all(>=(0),p.cells) && all(>=(0),p.work) && p.work[1,3]==0
    @test abs(after.energy_accounted-before.energy_accounted)<2e-8
    for (loci,q,mass) in ((1,lattice_parameters(),.32),(4,lattice_parameters(minimum_mass=.1,division_mass=.22,capacity=.25),.2))
        w=lattice_world(Grid(2,2,1));p=LatticeEcology(w;parameters=q,loci,founder_genotypes=1,initial_mass=mass)
        mktempdir() do dir
            path=joinpath(dir,"small.jls");lattice_checkpoint(w,p,path);v,r=restore_lattice(path)
            @test r.cells==p.cells && r.genotype==p.genotype && r.weights==p.weights
        end
    end
    amount,after,spent=LE.funded_flux(1f0,.0063005444f0,.35f0)
    @test amount>=0 && after>=0 && spent>=0 && spent<=.0063005444f0
    @test amount*.35f0<=spent
end
@testset "distinct wrappers cannot alias kernel ownership" begin
    w=lattice_world(Grid(4,4,1;extent=(4,4,1)));p=LatticeEcology(w;loci=4)
    p.scratch=reshape(vec(p.environment),size(p.environment))
    @test p.scratch!==p.environment && Base.mightalias(p.scratch,p.environment)
    @test_throws ArgumentError step!(w,p,.1)
    @test w.time==0 && p.tick==0
end
@testset "paid neural response remains bounded across work/time quantization" begin
    w=lattice_world(Grid(1,1,1));q=lattice_parameters(basal_work=0,neural_work=1e-40,growth=0,construction=0,erosion=0,uptake=0,feed=0)
    p=LatticeEcology(w;parameters=q,loci=1,founder_genotypes=1,occupancy=1,fuel=0,reserve=0)
    fill!(p.weights,0);p.weights[:,end,1:4].=1
    step!(w,p,.1)
    @test 0<maximum(p.hidden)<=1 && minimum(p.work)>=0
    w=lattice_world(Grid(1,1,1));q=lattice_parameters(response_time=1e7,growth=0,construction=0,erosion=0,uptake=0,feed=0)
    p=LatticeEcology(w;parameters=q,loci=1,founder_genotypes=1,occupancy=1,fuel=0,reserve=0)
    fill!(p.weights,0);p.weights[:,end,1:4].=1;step!(w,p,.1)
    @test 0<maximum(p.hidden)<1e-7
    @test_throws ArgumentError lattice_parameters(efficiency=1e-100)
end
@testset "stable transport scales before multiplying a large concentration jump" begin
    w=lattice_world(Grid(2,1,1;extent=(2,1,1)))
    q=lattice_parameters(fuel_diffusion=1e30,precursor_diffusion=0,feed=0)
    p=LatticeEcology(w;parameters=q,loci=1,founder_genotypes=1,occupancy=0,fuel=0,reserve=0)
    p.environment[1,1]=1f20;before=lattice_statistics(w,p);step!(w,p,1f-31);after=lattice_statistics(w,p)
    @test after.finite && after.minimum_physical>=0
    @test abs(after.material-before.material)<1e-7before.material
end
