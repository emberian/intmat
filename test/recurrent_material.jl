@testset "recurrent material composition" begin
    a=World(Grid(8,8,8;extent=(8,8,8));seed=31)
    b=World(a.grid;seed=31)
    p=RecurrentMaterial(a);q=RecurrentMaterial(b;feedback=0)
    original=a.fields.nutrient
    mass=sum(Float64,a.fields.nutrient)
    for _ in 1:100
        advance!(a,p,.05);advance!(b,q,.05)
    end
    @test all(isfinite,p.state)
    @test maximum(abs.(p.state.-q.state))>1e-3
    @test minimum(a.fields.nutrient)>=0 && minimum(a.fields.deposit)>=0
    @test isapprox(sum(Float64,a.fields.nutrient)+sum(Float64,a.fields.deposit),mass;atol=.001)
    @test a.fields.nutrient===original # even number of buffer swaps
    @test maximum(a.fields.deposit)>0
    @test record_process_fields(p).state===p.state
    @test process_ports(p).owned==(:state,)
    fine=World(Grid(8,8,8;extent=(.1,.1,.1)))
    r=RecurrentMaterial(fine)
    @test stability_limit(fine,r)<.001
    @test_throws ArgumentError advance!(fine,r,.05)
    invalid=World(Grid(2,2,2)); bad=RecurrentMaterial(invalid)
    bad.scratch=bad.state
    @test_throws ArgumentError advance!(invalid,bad,.05)
    invalid=World(Grid(2,2,2)); bad=RecurrentMaterial(invalid)
    invalid.scratch=merge(invalid.scratch,(nutrient=invalid.fields.nutrient,))
    @test_throws ArgumentError advance!(invalid,bad,.05)
    a.fields=merge(a.fields,(nutrient=copy(a.fields.nutrient),))
    @test_throws ArgumentError advance!(a,p,.05)
end

struct RecurrentAuditWrapper{P}
    child::P
end
Interior.process_children(p::RecurrentAuditWrapper)=(p.child,)

@testset "recurrent swaps reject fixed borrowers before mutation" begin
    for kind in (:medium,:boundary,:neural,:pool,:pool_deposit,:metabolic,:cellular,:wrapped_pool), reverse in (false,true)
        w=World(Grid(3,3,2;extent=(3,3,2)))
        p=RecurrentMaterial(w)
        phase=ImplicitPhaseSeparation(w;mineral_pinning=0,nutrient_diffusion=0,
            interface_activity=0,activity_relaxation=0)
        medium=ActiveMedium(w,phase)
        pool=ResourcePool(w;potentials=(4,3,0))
        # A reservoir may explicitly borrow a structural world array.
        pool.fields[1]=kind==:pool_deposit ? w.fields.deposit : w.fields.nutrient
        positions=zeros(Float32,3,0)
        borrower=kind==:medium ? medium : kind==:boundary ? BoundaryMaintenance(w,medium) :
            kind==:neural ? NeuralEcology(positions;medium=medium) : kind in (:pool,:pool_deposit) ? pool :
            kind==:metabolic ? MetabolicEcology(positions,pool) :
            kind==:cellular ? CellularEcology(w,positions,pool) : RecurrentAuditWrapper(pool)
        processes=reverse ? (borrower,p) : (p,borrower)
        before=deepcopy(w.fields); state=copy(p.state)
        nutrient=w.fields.nutrient; deposit=w.fields.deposit
        @test_throws ArgumentError step!(w,processes,.01f0)
        @test w.fields==before && p.state==state && w.time==0
        @test w.fields.nutrient===nutrient && w.fields.deposit===deposit
    end
    # Dynamic world-field readers and independent reservoirs remain composable.
    w=World(Grid(3,3,2;extent=(3,3,2)));p=RecurrentMaterial(w)
    independent=ResourcePool(w)
    before=sum(Float64,w.fields.nutrient)+sum(Float64,w.fields.deposit)
    for stages in ((p,Deposition()),(Deposition(),p),(p,Diffusion(:nutrient,.01)),(p,independent))
        step!(w,stages,.01f0)
    end
    @test w.time≈4Float64(.01f0)
    @test sum(Float64,w.fields.nutrient)+sum(Float64,w.fields.deposit)≈before atol=1e-5
end
