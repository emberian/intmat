using Test, Interior
include(joinpath(@__DIR__,"..","studies","multirate_resolved_mechanics.jl"))
using .MultirateResolvedMechanics
@testset "Resolved contacts with slow biological clock" begin
    x=resolved_comparison(size=8,founders=2,capacity=32,fast_steps=8,ratios=(1,2,4))
    @test first(x.rows).all_arrays_exact
    @test all(r->r.audit.finite,x.rows)
    @test all(r->r.actual_mechanical_substep==first(x.rows).actual_mechanical_substep,x.rows)
    @test all(r->r.bond_capture_events==8,x.rows)
    @test all(r->abs(r.material_error)<1e-4,x.rows)
    @test all(r->abs(r.potential_error)<1e-4,x.rows)
end
@testset "Fast recurrent clock and biological event boundary" begin
    x=resolved_comparison(size=8,founders=2,capacity=32,fast_steps=8,ratios=(1,2,4),fast_controller=true)
    @test first(x.rows).all_arrays_exact
    @test all(r->r.neural_updates==8,x.rows)
    @test all(r->r.audit.finite,x.rows)
    r=resolved_fixture(size=8,founders=2,capacity=16,ratio=4,fast_controller=true)
    s=r.process.cells.state
    s.mass[1]=.1f0;s.mass[2]=2f0;s.energy[2]=2f0
    before=Interior.device_habitat_audit(Interior.device_habitat_snapshot(r.process))
    resolved_step!(r)
    snapshot=Interior.device_habitat_snapshot(r.process)
    after=Interior.device_habitat_audit(snapshot)
    alive=findall(==(1),snapshot.cells.state.active)
    @test sort(snapshot.cells.state.ids[alive])==[2,3]
    child=only(filter(i->snapshot.cells.state.ids[i]==3,alive))
    @test snapshot.cells.state.parents[child]==2
    @test all(iszero,snapshot.arrays.last_secretion[child,:])
    @test abs(after.material-before.material)<1e-4
    @test abs(after.potential+after.work+after.heat-before.potential-before.work-before.heat)<1e-4
end
@test_throws ErrorException resolved_comparison(fast_steps=0)
@testset "Fast bonded neural and funded mechanical coupling" begin
    r=resolved_fixture(size=8,founders=2,capacity=16,ratio=4,fast_controller=true)
    s=r.process.cells.state;h=Int(r.process.cells.geometry.hidden)
    s.position[1,:].=Float32[1,1,1];s.position[2,:].=Float32[1.32,1,1]
    for (i,j) in ((1,2),(2,1))
        s.bonds[i,1]=j;s.bond_ids[i,1]=s.ids[j];s.rest[i,1]=.3f0
        s.weights[i,:,h+8].=0;s.weights[i,end,h+8]=2
    end
    s.hidden[1,1]=.5f0;s.hidden[2,1]=-.5f0
    before=Interior.device_habitat_audit(Interior.device_habitat_snapshot(r.process))
    resolved_step!(r)
    after=Interior.device_habitat_audit(Interior.device_habitat_snapshot(r.process))
    @test s.bonds[1,1]==2&&s.bonds[2,1]==1
    @test s.hidden[1,1]!=.5f0&&s.hidden[2,1]!=-.5f0
    @test all(>=(0),s.energy)
    @test after.work>before.work
    @test abs(after.material-before.material)<1e-4
    @test abs(after.potential+after.work+after.heat-before.potential-before.work-before.heat)<1e-4
end
