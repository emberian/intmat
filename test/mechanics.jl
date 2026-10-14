@testset "Deformable persistent bond graph" begin
    w=World(Grid(12,12,12));p=BondNetwork(w;size=(4,4,4),jitter=0,damping=2)
    ids=copy(p.bond_ids);nodes=copy(p.node_ids);rest=copy(p.position)
    @test maximum(abs,bond_forces(p))<1e-12
    step!(w,p,.001)
    @test p.position==rest
    p.position[1,22]+=.006
    @test maximum(abs,sum(bond_forces(p);dims=2))<1e-12
    energy=mechanical_energy(p).total
    for i in 1:500;step!(w,p,.001);end
    state=mechanical_energy(p)
    @test state.elastic+state.kinetic<energy
    @test state.damping>0
    @test abs(state.total-energy)<energy*.015
    @test p.bond_ids==ids && p.node_ids==nodes
    q=BondNetwork(w;size=(3,3,1),jitter=0,failure_strain=.01)
    q.position[1,:].*=1.15
    before=mechanical_energy(q).total
    break_bonds!(q)
    @test count(==(1),q.damage)>0
    @test mechanical_energy(q).total≈before
    damaged=copy(q.damage);q.position.=q.reference;break_bonds!(q)
    @test q.damage==damaged
    @test q.damage_energy>0
    @test record_process_fields(q).bond_ids==q.bond_ids
    @test_throws ArgumentError step!(w,p,2stability_limit(w,p))
    p.mass[1]=-1
    @test_throws ArgumentError step!(w,p,.001)
    p.mass[1]=.02;p.damage[1]=1.1
    @test_throws ArgumentError step!(w,p,.001)
end
