function catalytic_fixture(;habitat_options...)
    w=World(Grid(8,8;extent=(4.,4.,1.));seed=19);fill!(w.fields.nutrient,0)
    r=ResourcePool(w;potentials=(4,4,0),initial=(0,0,0),diffusion=(0,0,0))
    c=CellularEcology(w,reshape([2.,2.,0.],3,1),r;energy=2.,metabolites=.3,
        speed=0.,repulsion=0.,stiffness=0.,basal_rate=0.,neural_cost=0.,uptake_rate=0.,growth_rate=0.,catabolic_rate=0.,
        division_mass=10.,mutation=0.,developmental_partition=0.)
    f=ResourcePool(w;names=(:bound_a,:bound_b),potentials=(6,5),initial=(.5,.5),diffusion=(0,0))
    h=CatalyticHabitat(w,c,f;merge((;regulation=false,investment_bias=0.,secretion_rate=1.,catalytic_rate=0.,decay_rate=0.),(;habitat_options...))...)
    w,r,c,f,h
end
catalytic_mass(w,r,c,f,h)=Interior.arena_inventory(w,(f,r,c,h))
catalytic_energy(r,c,f,h)=Interior.arena_energy((f,r,c,h))

@testset "Paid catalytic material and complete ledgers" begin
    w,r,c,f,h=catalytic_fixture()
    mass=catalytic_mass(w,r,c,f,h);energy=catalytic_energy(r,c,f,h)
    stores=copy(c.metabolites);battery=only(c.energy)
    step!(w,h,.2)
    deposited=sum(h.secreted)
    @test deposited>0
    @test sum(sum,h.enzymes)≈deposited atol=1e-14
    @test stores.-c.metabolites≈fill(deposited/2,2,1) atol=1e-14
    @test battery-only(c.energy)≈deposited*h.parameters.secretion_work atol=1e-14
    @test c.chemical_heat≈4deposited atol=1e-14
    @test catalytic_mass(w,r,c,f,h)≈mass atol=1e-12
    @test catalytic_energy(r,c,f,h)≈energy atol=1e-12
    @test Interior.arena_inventory(w,(h,))≈mass atol=1e-12
    @test Interior.arena_energy((h,))≈energy atol=1e-12
    @test Interior.arena_inventory(w,(h,c,r,f))≈mass atol=1e-12
    @test Interior.arena_energy((h,c,r,f))≈energy atol=1e-12

    # No precursor or no work means no free catalyst construction.
    for lacking in (:precursor,:work)
        w,r,c,f,h=catalytic_fixture()
        lacking==:precursor ? (c.metabolites[1,1]=0.) : (c.energy[1]=0.)
        step!(w,h,.2)
        @test sum(h.secreted)==0
        @test all(>=(0),c.metabolites) && all(>=(0),c.energy)
    end
end

@testset "Catalytic reservoirs are independent finite inventories" begin
    w,r,c,f,h=catalytic_fixture()
    f.fields[1]=r.fields[1]
    @test_throws ArgumentError step!(w,h,.2)
    @test sum(h.secreted)==0 && w.time==0
    w,r,c,f,h=catalytic_fixture()
    h.release[1]=h.enzymes[1]
    @test_throws ArgumentError step!(w,h,.2)
    w,r,c,f,h=catalytic_fixture()
    f.fields[1]=w.fields.nutrient
    @test_throws ArgumentError step!(w,h,.2)
    w,r,c,f,h=catalytic_fixture()
    r.potentials[1:2].=3f38;f.potentials.=3f38
    before=catalytic_energy(r,c,f,h)
    step!(w,h,.2)
    @test isfinite(c.chemical_heat)
    @test catalytic_energy(r,c,f,h)≈before rtol=1e-14
end

@testset "Catalytic secretion never overdraws tiny donors" begin
    w,r,c,f,h=catalytic_fixture(;secretion_rate=100.)
    for stock in (1e-18,5.421010862427522e-20,.000123456789,.3)
        c.metabolites[:,1].=stock;c.energy[1]=stock
        step!(w,h,.2)
        @test all(>=(0),c.metabolites)
        @test all(>=(0),c.energy)
    end
end

@testset "Catalytic unlocking and inert decay" begin
    w,r,c,f,h=catalytic_fixture(;secretion_rate=0.,catalytic_rate=100.)
    mass=catalytic_mass(w,r,c,f,h);step!(w,h,.2)
    @test sum(h.solubilized)==0 && sum(sum,r.fields)==0
    # A known material catalyst unlocks matching chemistry without cell uptake.
    h.enzymes[1][4,4,1]=1.
    mass=catalytic_mass(w,r,c,f,h);energy=catalytic_energy(r,c,f,h)
    step!(w,h,.2)
    @test h.solubilized[1]>.49 && h.solubilized[2]==0
    @test f.fields[1][4,4,1]==0
    @test c.chemical_heat≈2h.solubilized[1] atol=1e-12
    @test catalytic_mass(w,r,c,f,h)≈mass atol=1e-12
    @test catalytic_energy(r,c,f,h)≈energy atol=1e-12
    @test only(c.energy)==2.
    h.parameters=merge(h.parameters,(catalytic_rate=0.,decay_rate=2.))
    before=sum(sum,h.enzymes);step!(w,h,.2)
    @test sum(sum,h.enzymes)≈before*exp(-2Float64(Float32(.2))) atol=1e-14
    @test sum(h.degraded)>0 && sum(r.fields[3])>0
    @test catalytic_mass(w,r,c,f,h)≈mass atol=1e-12
    @test catalytic_energy(r,c,f,h)≈energy atol=1e-12
end

@testset "Inherited catalytic regulation survives birth and reorder" begin
    w,r,c,f,h=catalytic_fixture(;regulation=true,mutation=0.)
    original=copy(h.weights)
    Interior._cellular_divide!(w,c,[1]);rebuild_cell_stencil!(c.stencil,c.positions)
    step!(w,h,.1)
    @test h.ids==c.ids==[1,2]
    @test h.weights[:,:,1]==original[:,:,1]==h.weights[:,:,2]
    reorder_ecology!(c,[2,1]);step!(w,h,.1)
    @test h.ids==[2,1] && h.weights[:,:,1]==original[:,:,1]
    h.parameters=merge(h.parameters,(mutation=.1,mutation_probability=1.))
    Interior._cellular_divide!(w,c,[1]);rebuild_cell_stencil!(c.stencil,c.positions)
    step!(w,h,.1)
    @test h.ids==c.ids
    @test h.weights[:,:,3]!=h.weights[:,:,1]
    @test h.weights[:,:,1]==original[:,:,1]
end

@testset "Catalytic composition, observations and exact restart" begin
    cfg=Dict{String,Any}("model"=>"catalytic","seed"=>51,"grid"=>Dict("shape"=>[16,16,1],"extent"=>[8.,8.,1.]),
        "cellular"=>Dict("population"=>8,"capacity"=>64),"feedstock"=>Dict("potentials"=>[6.,5.]))
    a=build_arena(cfg)
    for _ in 1:20;arena_step!(a);end
    s=arena_snapshot(a)
    @test length(s["agents"])==s["metrics"]["population"]
    @test all(haskey(s["fields"],key) for key in ("enzyme_a","enzyme_b","release_a","release_b","bound_a","bound_b"))
    @test all(haskey(ag,"catalysis") for ag in s["agents"])
    @test abs(s["metrics"]["material_error"])<1e-4
    @test abs(s["metrics"]["chemical_energy_error"])<1e-4
    @test sum(s["metrics"]["catalytic"]["secreted"])>0
    mktempdir() do directory
        path=joinpath(directory,"habitat.jls");checkpoint!(a,path);b=restore_arena(path)
        @test b.processes[end].cells===b.processes[3]
        @test b.processes[end].substrates===b.processes[1]
        for _ in 1:15;arena_step!(a);arena_step!(b);end
        x=arena_snapshot(a);y=arena_snapshot(b)
        @test x["metrics"]==y["metrics"]
        @test x["agents"]==y["agents"]
        @test x["fields"]==y["fields"]
    end
end
