using TOML

@testset "Finite-capacity sorption kinetics" begin
    f=Interior.sorption_target
    for (M,C,kon,koff) in ((2.,1.,2.,.1),(1.,1.,2.,0.),(1e-8,8.,1.,3.),(3.,2.,0.,.2))
        B=min(M,C)*.7
        half=f(B,M,C,kon,koff,.25,.3)
        @test 0<=half<=min(M,C)
        @test f(half,M,C,kon,koff,.25,.3)≈f(B,M,C,kon,koff,.25,.6) atol=1e-13
        eq=f(B,M,C,kon,koff,.25,1e5)
        @test abs(kon/.25*(M-eq)*(C-eq)-koff*eq)<1e-8
    end
end

@testset "Sorption owns real chemical inventory" begin
    w,r,c,f,h=catalytic_fixture(;secretion_rate=0.,catalytic_rate=0.)
    h.enzymes[1][4,4,1]=.2;r.fields[1][4,4,1]=.5
    p=SorptiveMatrix(w,h;capacity_per_enzyme=2.,adsorption_rate=5.,desorption_rate=.1)
    a=Arena(w,(f,r,c,h,p);dt=.2)
    mass=a.initial_mass;energy=Interior.arena_energy(a.processes)
    step!(w,p,.2)
    @test 0<p.bound[1][4,4,1]<=.4
    @test r.fields[1][4,4,1]<.5
    @test Interior.arena_inventory(w,a.processes)≈mass atol=1e-12
    @test Interior.arena_energy(a.processes)≈energy atol=1e-12
    @test Interior.arena_inventory(w,(p,))≈mass atol=1e-12
    @test Interior.arena_energy((p,))≈energy atol=1e-12
    @test Interior.arena_energy((p,h,c,r,f))≈energy atol=1e-12
    # The catalyst is removed here to isolate capacity collapse, not destroyed
    # in the simulation law (normal degradation deposits the same mass as waste).
    h.enzymes[1][4,4,1]=0
    before=p.bound[1][4,4,1];step!(w,p,.2)
    @test p.displaced[1]>before-1e-7
    @test p.bound[1][4,4,1]<1e-7
    @test p.desorbed[1]>0 && all(>=(0),p.bound[1])
    @test Interior.arena_inventory(w,a.processes)≈mass-.2 atol=1e-12
    @test Interior.arena_energy(a.processes)≈energy atol=1e-12
end

@testset "Constructed resistance composes without accumulating" begin
    w,r,c,f,h=catalytic_fixture(;secretion_rate=0.,catalytic_rate=0.)
    h.enzymes[1].=.1
    p=SorptiveMatrix(w,h;capacity_per_enzyme=0.,transport_barrier=5.)
    @test maximum(p.permeability)<1 # existing material observed at construction
    Interior.cellular_porosity!(w,c);body=copy(r.permeability)
    step!(w,p,.1);once=copy(r.permeability)
    @test once≈body./3 atol=1e-7 # site volume=.25, resistance=1+5*.1/.25
    step!(w,p,.1)
    @test r.permeability==once
    @test sum(sum,p.bound)==0
    p.parameters=merge(p.parameters,(transport_barrier=0.,))
    step!(w,p,.1)
    @test r.permeability==body
    p.bound[1]=h.enzymes[1]
    @test_throws ArgumentError step!(w,p,.1)
    p=SorptiveMatrix(w,h;transport_barrier=1e300)
    step!(w,p,.1);step!(w,p,.1)
    @test all(iszero,p.permeability)
    p.permeability=f.scratch
    @test_throws ArgumentError step!(w,p,.1)
end

@testset "Full catalytic-matrix composition and checkpoint" begin
    cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","arena","catalytic-mosaic.toml"))
    cfg["grid"]["shape"]=[12,12,1];cfg["cellular"]["population"]=6
    cfg["matrix"]=Dict("capacity_per_enzyme"=>4.,"transport_barrier"=>10.)
    a=build_arena(cfg)
    for _ in 1:30;arena_step!(a);end
    snap=arena_snapshot(a);p=last(a.processes)
    @test p isa SorptiveMatrix
    @test sum(sum,p.bound)>0 && minimum(p.permeability)<1
    @test all(haskey(snap["fields"],name) for name in ("retained_a","retained_b","enzyme_a","precursor_a","matrix_permeability"))
    @test abs(snap["metrics"]["material_error"])<1e-5
    @test abs(snap["metrics"]["chemical_energy_error"])<1e-4
    mktempdir() do directory
        checkpoint!(a,joinpath(directory,"state.jls"));b=restore_arena(joinpath(directory,"state.jls"))
        @test last(b.processes).habitat===b.processes[end-1]
        for _ in 1:20;arena_step!(a);arena_step!(b);end
        sa=arena_snapshot(a);sb=arena_snapshot(b)
        delete!(sa,"generated_at_unix");delete!(sb,"generated_at_unix")
        @test sa==sb
    end
end

@testset "Capacity loss cannot leave food on absent scaffold" begin
    w,r,c,f,h=catalytic_fixture(;secretion_rate=0.,catalytic_rate=0.)
    p=SorptiveMatrix(w,h);p.bound[1][1]=1.;r.fields[1][1]=1f20
    before=Interior.arena_energy((p,));step!(w,p,.1)
    @test p.bound[1][1]==0
    @test p.inert_residual==1
    @test c.chemical_heat==4
    @test Interior.arena_energy((p,))==before
end
