@testset "Shared anabolism has one material and work budget" begin
    for (stores,energy) in ((.05,2.),(.3,.005),(0.,2.),(.3,0.),(1e-18,1e-18))
        w,r,c,f,h=catalytic_fixture()
        c.metabolites[:,1].=stores;c.energy[1]=energy
        p=ConstructiveEcology(w,c,h)
        allocator=Interior.CatalyticAnabolism(p)
        mass=Interior.arena_inventory(w,(p,));potential=Interior.arena_energy((p,));body=only(c.biomass)
        predicted=Interior._cellular_planned_work(.16,w,c,1,.2,c.parameters,allocator)
        @test predicted≈.17 atol=1e-14
        Interior._cellular_anabolic_build!(w,c,1,.2,.4,c.parameters,allocator)
        @test Interior.arena_inventory(w,(p,))≈mass atol=1e-12
        @test Interior.arena_energy((p,))≈potential atol=1e-12
        @test all(>=(0),c.metabolites) && all(>=(0),c.energy)
        @test p.growth_planned+sum(p.construction_planned)<=2stores
        @test c.parameters.growth_work*p.growth_planned+h.parameters.secretion_work*sum(p.construction_planned)<=energy
        if stores>1e-10 && energy>0
            @test sum(h.secreted)>0 && only(c.biomass)>body
            @test sum(p.construction_planned)/p.growth_planned≈.125 atol=1e-14
        else
            @test sum(h.secreted)<=2stores
        end
    end
end

@testset "Construction ownership, terminal cells and observations" begin
    w,r,c,f,h=catalytic_fixture()
    p=ConstructiveEcology(w,c,h)
    @test_throws ArgumentError Arena(w,(r,f,p,h);dt=.1)
    @test_throws ArgumentError Arena(w,(r,f,p,c);dt=.1)
    c.energy.=0;c.metabolites.=0
    c.parameters=merge(c.parameters,(basal_rate=1.,starvation_rate=10.))
    a=Arena(w,(f,r,p);dt=.2);arena_step!(a)
    @test isempty(c.ids) && isempty(h.ids) && c.deaths==1
    @test isfinite(Interior.arena_energy(a.processes))
    @test abs(arena_snapshot(a)["metrics"]["material_error"])<1e-8
end

@testset "Shared construction full composition and exact restart" begin
    cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","arena","sorptive-refuges.toml"))
    cfg["grid"]["shape"]=[12,12,1];cfg["cellular"]["population"]=6
    cfg["catalysis"]["shared_anabolism"]=true
    a=build_arena(cfg)
    for _ in 1:100;arena_step!(a);end
    p=only(filter(x->x isa ConstructiveEcology,a.processes));s=arena_snapshot(a)
    @test p.cells===p.habitat.cells && last(a.processes).habitat===p.habitat
    @test sum(p.habitat.secreted)>0 && sum(p.construction_planned)>0
    @test haskey(s["metrics"],"anabolism") && all(haskey(x,"catalysis") for x in s["agents"])
    @test abs(s["metrics"]["material_error"])<1e-4 && abs(s["metrics"]["chemical_energy_error"])<1e-4
    mktempdir() do directory
        path=joinpath(directory,"shared.jls");checkpoint!(a,path);b=restore_arena(path)
        for _ in 1:50;arena_step!(a);arena_step!(b);end
        sa=arena_snapshot(a);sb=arena_snapshot(b)
        delete!(sa,"generated_at_unix");delete!(sb,"generated_at_unix")
        @test sa==sb
    end
end

@testset "mixed death division and reorder identity" begin
 w=World(Grid(8,8;extent=(4.,4.,1.));seed=73)
 r=ResourcePool(w;potentials=(4,4,0),initial=(0,0,0),diffusion=(0,0,0))
 c=CellularEcology(w,[1. 3.;1. 3.;0. 0.],r;energy=1.,metabolites=.1,capacity=8,
  speed=0.,repulsion=0.,stiffness=0.,basal_rate=0.,neural_cost=0.,uptake_rate=0.,growth_rate=0.,catabolic_rate=0.,
  division_mass=.5,minimum_mass=.1,mutation=0.,mutation_probability=0.,developmental_partition=0.)
 fill!(c.weights,0f0);c.weights[size(c.hidden,1)+10,end,:].=1f0
 f=ResourcePool(w;names=(:bound_a,:bound_b),potentials=(6,5),initial=(.1,.1),diffusion=(0,0))
 h=CatalyticHabitat(w,c,f;secretion_rate=0.,catalytic_rate=0.,decay_rate=0.,mutation=0.,mutation_probability=0.)
 firstgen=copy(h.weights[:,:,1]);survivorgen=copy(h.weights[:,:,2]);p=ConstructiveEcology(w,c,h)
 c.biomass[1]=.01;c.metabolites[:,1].=0;c.energy[1]=0
 reorder_ecology!(c,[2,1])
 step!(w,p,.1)
 @test c.ids==h.ids==[2,3]
 @test c.parents==[0,2]
 @test c.deaths==1 && c.births==1
 @test h.weights[:,:,1]==h.weights[:,:,2]==survivorgen
 @test h.weights[:,:,1]!=firstgen
 @test h.last_secretion==zeros(2,2)
 reorder_ecology!(c,[2,1]);step!(w,p,.1)
 @test issorted(c.ids) && c.ids==h.ids
 @test all(i->h.weights[:,:,i]==survivorgen,eachindex(c.ids))
end
