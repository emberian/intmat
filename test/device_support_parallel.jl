using Interior,KernelAbstractions,Test

function parallel_support_fixture(shape=(8,8,8);execution=:serial,inventory_type=Float64)
    w=World(Grid(shape...;extent=(2.,2.,shape[3]==1 ? 1. : 2.)),CPU(),NamedTuple(),NamedTuple(),0.,721)
    xyz=Float64[0 .9999 1.9999;0 .51 1.9999;0 (shape[3]==1 ? 0 : 1.1) (shape[3]==1 ? 0 : 1.9999)]
    p=device_cellular_initialize(w,xyz;capacity=8,hidden=2,seed=721,biomass=2.,energy=2.,
        energy_type=Float64,inventory_type,support_strategy=:materialized,support_execution=execution,
        exchange_radius=.18,mechanical_substeps=2,radius_bound=1.,speed=0.,
        repulsion=0.,stiffness=0.,mutation=0.,neural_cost=0.,basal_rate=.001,
        minimum_mass=.1,division_mass=1.7,initial=(.1,.1,0.),diffusion=(.001,.001,.001))
    w,p
end
function support_arrays_equal(a,b)
    aa=device_cellular_snapshot(a);bb=device_cellular_snapshot(b)
    isequal(aa,bb)
end
@testset "Native ordered parallel materialized support" begin
    for shape in ((4,4,1),(4,4,4),(16,16,16)),T in (Float32,Float64)
        wa,a=parallel_support_fixture(shape;inventory_type=T)
        wb,b=parallel_support_fixture(shape;execution=:parallel,inventory_type=T)
        @test Interior.device_support_execution(a)===:serial
        @test Interior.device_support_workspace(a)===nothing
        @test Interior.device_support_execution(b)===:parallel
        @test support_arrays_equal(a,b)
        # Preserve stale tails and inactive slots, including fallback on coarse grids.
        a.state.sites[:,end].=19;b.state.sites[:,end].=19
        a.state.support[:,end].=.123f0;b.state.support[:,end].=.123f0
        for p in (a,b)
            Interior.device_cellular_launch!(p,Interior.dc_support!,p.state,p.environment,p.parameters,p.geometry;ndrange=8)
        end
        @test support_arrays_equal(a,b)
        for _ in 1:8
            device_cellular_step!(wa,a,.01);device_cellular_step!(wb,b,.01)
        end
        @test support_arrays_equal(a,b)
        @test a.environment.metadata[3]>0 # real division/composition, not a controller-only check
        cp=device_cellular_checkpoint(wb,b)
        @test cp.support_execution===:parallel
        wr=World(wb.grid,CPU(),NamedTuple(),NamedTuple(),wb.time,wb.seed)
        restored=device_cellular_restore(wr,cp)
        @test Interior.device_support_execution(restored)===:parallel
        @test Interior.device_support_workspace(restored)!==Interior.device_support_workspace(b)
        @test support_arrays_equal(b,restored)
        device_cellular_step!(wb,b,.01);device_cellular_step!(wr,restored,.01)
        @test support_arrays_equal(b,restored)
        @test_throws ArgumentError device_cellular_restore(wr,merge(cp,(support_execution=:invalid,)))
        @test_throws ArgumentError Interior.device_cellular_launch!(b,Interior.dc_support!,b.state,b.environment,b.parameters,b.geometry;ndrange=7)
    end
end
@testset "Support workspace ownership and invalid choices" begin
    w,p=parallel_support_fixture(;execution=:parallel)
    work=Interior.device_support_workspace(p);resource=Interior.device_resource_workspace(p);lifecycle=Interior.device_lifecycle_workspace(p)
    @test Interior.device_support_workspace_validate(work,p.backend,p.geometry)
    @test Interior.device_support_workspace_nonalias(work,p.state,p.environment,resource,lifecycle)
    aliased=merge(work.arrays,(sites=p.state.sites,))
    bad=Interior.DeviceSupportWorkspace(aliased,p.backend,(arrays=aliased,shapes=map(size,aliased),dims=work.contract.dims))
    @test_throws ArgumentError Interior.device_cellular_contract(p.state,p.environment,p.parameters,p.geometry;
        support_execution=:parallel,support_workspace=bad)
    @test_throws ArgumentError Interior.device_cellular_contract(p.state,p.environment,p.parameters,p.geometry;
        support_workspace=work)
    resize!(work.arrays.totals,7)
    @test_throws ArgumentError Interior.device_cellular_validate(w,p)
    @test_throws ArgumentError parallel_support_fixture(;execution=:invalid)
    @test_throws ArgumentError device_cellular_initialize(w,zeros(3,0);support_strategy=:axis,support_execution=:parallel)
end
