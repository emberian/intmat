# Standalone until experimental DeviceCellular is integrated into Interior.
using Interior, Test
isdefined(Interior,:DeviceCellular) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_cellular.jl"))
function device_test_fixture(n=0)
    w=World(Grid(2,2;extent=(2,2,1)));fill!(w.fields.nutrient,0)
    r=ResourcePool(w;potentials=(4,4,0),initial=(.2,.3,.4),diffusion=(0,0,0))
    c=CellularEcology(w,zeros(3,n),r;capacity=2)
    w,Interior.DeviceCellular(w,c;mechanical_substeps=1)
end
@testset "Device cellular contracts and failure boundary" begin
    w,p=device_test_fixture()
    @test_throws ArgumentError Interior.device_cellular_step!(w,p,1e-100)
    original=p.state
    p.state=merge(original,(position=original.position_next,))
    @test_throws ArgumentError Interior.device_cellular_step!(w,p,.01)
    p.state=original
    p.environment.status[1]=8
    before=Interior.device_cellular_snapshot(p;check_status=false)
    for _ in 1:3;Interior.device_cellular_step!(w,p,.01;check_status=false);end
    after=Interior.device_cellular_snapshot(p;check_status=false)
    @test all(getproperty(before.state,k)==getproperty(after.state,k) for k in keys(before.state))
    @test all(getproperty(before.environment,k)==getproperty(after.environment,k) for k in keys(before.environment))
    @test after.environment.metadata[6]==0
    @test p.attempted_steps==3
end
@testset "Device precision and antisymmetric periodic forces" begin
    w,p=device_test_fixture();e=p.environment
    e.resource_ledger[1,1]=1f6
    for _ in 1:100;Interior.dc_resource_add!(e,1,1,-.01f0);end
    @test abs(Float64(e.resource_ledger[1,1])+Float64(e.resource_ledger_low[1,1])-(1e6-100Float64(.01f0)))<1e-8
    @test e.status[1]==0
    e.fields[1,3]=3f38;e.returns[1]=3f38
    Interior.dc_commit_waste!(CPU(),128)(e;ndrange=4)
    @test e.status[1]==8
    @test isfinite(e.fields[1,3])
    positions=(position=Float32[0 0 0;1 0 0],)
    g=(extent=(2f0,2f0,1f0),nd=2,)
    @test Interior.dc_delta(positions,g,1,2,1)==-Interior.dc_delta(positions,g,2,1,1)
end
@testset "Device complete checkpoint and deterministic continuation" begin
    w,p=device_test_fixture(1)
    Interior.device_cellular_step!(w,p,.01)
    checkpoint=Interior.device_cellular_checkpoint(w,p)
    w2=World(w.grid);fill!(w2.fields.nutrient,0)
    p2=Interior.device_cellular_restore(w2,checkpoint)
    before=Interior.device_cellular_snapshot(p)
    restored=Interior.device_cellular_snapshot(p2)
    @test before==restored
    badstate=merge(checkpoint.snapshot.state,(active=fill(Int32(2),length(p.state.active)),))
    bad=merge(checkpoint,(snapshot=merge(checkpoint.snapshot,(state=badstate,)),))
    @test_throws ArgumentError Interior.device_cellular_restore(w2,bad)
    badgeom=merge(checkpoint.snapshot.geometry,(spacing=(.1f0,.1f0,.1f0),))
    bad=merge(checkpoint,(snapshot=merge(checkpoint.snapshot,(geometry=badgeom,)),))
    @test_throws ArgumentError Interior.device_cellular_restore(w2,bad)
    @test_throws ArgumentError Interior.device_cellular_restore(w2,merge(checkpoint,(julia_version=v"0.0.0",)))
    for _ in 1:3
        Interior.device_cellular_step!(w,p,.01)
        Interior.device_cellular_step!(w2,p2,.01)
    end
    @test Interior.device_cellular_snapshot(p)==Interior.device_cellular_snapshot(p2)
    p.environment.status[1]=8
    @test_throws ArgumentError Interior.device_cellular_checkpoint(w,p)
end
@testset "Device mechanics timestep convergence" begin
    w=World(Grid(8,8;extent=(8,8,1)));fill!(w.fields.nutrient,0)
    r=ResourcePool(w;potentials=(4,4,0),initial=(0,0,0),diffusion=(0,0,0))
    c=CellularEcology(w,Float64[2 3.5;2 2;0 0],r;capacity=2)
    add_cell_bond!(c.bonds,1,2,1.)
    p=Interior.DeviceCellular(w,c);s=p.state;e=p.environment;g=p.geometry
    q=merge(p.parameters,(speed=0f0,repulsion=0f0,stiffness=1f0))
    initial=copy(s.position);errors=Float64[]
    target=.5exp(-.2)
    for substeps in (2,4,8)
        copyto!(s.position,initial)
        for _ in 1:substeps
            Interior.device_cellular_launch!(p,Interior.dc_clear_heads!,e;ndrange=length(e.heads))
            Interior.device_cellular_launch!(p,Interior.dc_geometry!,s,e,q,g;ndrange=2)
            Interior.device_cellular_launch!(p,Interior.dc_forces!,s,e,q,g,.1f0/substeps;ndrange=2)
            Interior.device_cellular_launch!(p,Interior.dc_move!,s,e,q,g,.1f0/substeps;ndrange=2)
            Interior.device_cellular_launch!(p,Interior.dc_commit_positions!,s,e;ndrange=2)
        end
        Interior.device_cellular_status!(p)
        push!(errors,abs(Float64(s.position[2,1]-s.position[1,1]-1)-target))
    end
    @test errors[3]<errors[2]<errors[1]
    @test errors[3]<.0011
end
@testset "Projected complete-step execution" begin
    w=World(Grid(8,8;extent=(8,8,1)));fill!(w.fields.nutrient,0)
    r=ResourcePool(w;potentials=(4,4,0),initial=(.3,.3,0),diffusion=(.02,.02,.01))
    c=CellularEcology(w,Float64[2 3.2 6;2 2 6;0 0 0],r;capacity=8,biomass=2,mutation=.01)
    c.biomass[3]=.1
    h=size(c.hidden,1);c.weights[h+8,end,:].=2;c.weights[h+9,end,:].=2;c.weights[h+10,end,:].=2
    a=Interior.DeviceCellular(w,c)
    saved=Interior.device_cellular_checkpoint(w,a)
    b=Interior.device_cellular_restore(w,saved)
    for _ in 1:20
        Interior.device_cellular_step!(w,a,.01;execution=:staged)
        Interior.device_cellular_step!(w,b,.01;execution=:projected)
    end
    original=Interior.device_cellular_snapshot(a)
    @test Interior.device_cellular_snapshot(b)==original
    for execution in (:staged,:projected)
        copyto!(b.environment.status,Int32[8])
        before=Interior.device_cellular_snapshot(b;check_status=false)
        Interior.device_cellular_step!(w,b,.01;execution=execution,check_status=false)
        after=Interior.device_cellular_snapshot(b;check_status=false)
        @test before.state==after.state && before.environment==after.environment
    end
end
@testset "Rounded device support layout and empty World" begin
    grid=Grid(8,8;extent=(8,8,1));w=World(grid,CPU(),NamedTuple(),NamedTuple(),0.,17)
    r=ResourcePool(w;potentials=(4,4,0),initial=(.2,.3,.4),diffusion=(0,0,0))
    c=CellularEcology(w,zeros(3,1),r;capacity=2,exchange_radius=1+1e-12)
    p=Interior.DeviceCellular(w,c)
    @test p.geometry.candidates==(4,4,1)
    saved=Interior.device_cellular_checkpoint(w,p)
    restored=Interior.device_cellular_restore(w,saved)
    @test Interior.device_cellular_snapshot(restored)==saved.snapshot
    c.resources.potentials[3]=1
    @test_throws ArgumentError Interior.DeviceCellular(w,c)
end
