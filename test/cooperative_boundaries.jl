@testset "Funded shared boundary film" begin
    function fixture()
        w=World(Grid(5,5;extent=(5,5,1)))
        phase=ImplicitPhaseSeparation(w;mean=.2,noise=0,mineral_pinning=0,
            nutrient_diffusion=0,activity_relaxation=0)
        medium=ActiveMedium(w,phase;interface_barrier=0,signal_diffusion=0,
            nutrient_diffusion=.2,membrane_barrier=6,membrane_pinning=.5)
        b=BoundaryMaintenance(w,medium;work_per_mass=2,decay=.1)
        w,phase,medium,b
    end
    w,p,m,b=fixture()
    @test maintain_boundary!(b,(3,3,1),1,0,1)==(material=0.0,work=0.0)
    debit=maintain_boundary!(b,(3,3,1),1,1,1)
    @test debit.material==.5 && debit.work==1
    @test boundary_material_mass(b)==.5
    @test local_boundary_benefit(b,(4,3,1))>0
    @test local_boundary_benefit(b,(1,1,1))==0
    initial=material_mass(w)+boundary_material_mass(b)
    deposit=copy(w.fields.deposit)
    step!(w,b,.1)
    @test w.fields.deposit==deposit
    @test material_mass(w)+boundary_material_mass(b)≈initial atol=1e-5
    @test b.recycled_material≈.5-boundary_material_mass(b)
    @test m.permeability[3,3,1]<m.permeability[1,1,1]
    @test p.chemical_forcing[3,3,1]>p.chemical_forcing[1,1,1]
    before=sum(p.phi)
    step!(w,p,.1)
    @test sum(p.phi)≈before atol=1e-5
    @test maximum(p.phi)>minimum(p.phi)
    # Retention benefits the neighboring location as well as the payer's cell.
    wf,pf,mf,bf=fixture(); wu,pu,mu,bu=fixture()
    maintain_boundary!(bf,(3,3,1),1,2,1)
    fill!(wf.fields.nutrient,0); fill!(wu.fields.nutrient,0)
    wf.fields.nutrient[4,3,1]=1; wu.fields.nutrient[4,3,1]=1
    step!(wf,bf,.1); step!(wu,bu,.1)
    @test wf.fields.nutrient[4,3,1]>wu.fields.nutrient[4,3,1]
    @test bf.spent_work==2
    # Allocation is limited by representable addition, never charged for zero.
    bf.film[3,3,1]=1f8
    @test maintain_boundary!(bf,(3,3,1),1,2,1).material==0
    @test_throws ArgumentError maintain_boundary!(bf,(3,3,1),1,-1,1)
end
