using Interior,Test,KernelAbstractions
isdefined(Interior,:DeviceProfileCarrier) || Base.include(Interior,joinpath(@__DIR__,"..","studies","device_launch_profile.jl"))
isdefined(Interior,:DeviceCatalyticHabitat) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_habitat.jl"))
isdefined(Interior,:DeviceSupportFastCarrier) || Base.include(Interior,joinpath(@__DIR__,"..","studies","device_support_fastpaths.jl"))
isdefined(Main,:constructive_device_fixture) || include(joinpath(@__DIR__,"..","scripts","device_habitat_fixture.jl"))

@testset "Exact-zero support shortcuts do not traverse support" begin
    # Any descriptor access would fail: zero export must return first.
    @test Interior.dh_support_fast_export!(nothing,1,1,0.)===0.
    # A zero-waste kernel can omit a deliberately inaccessible descriptor.
    s=(waste=zeros(Float32,3),);e=(status=zeros(Int32,1),)
    Interior.dc_support_fast_return_waste!(CPU(),32)(s,e,nothing;ndrange=3)
    KernelAbstractions.synchronize(CPU())
    @test e.status==[0] && all(iszero,s.waste)
end

@testset "Support fastpaths preserve complete coupled trajectories" begin
    for strategy in (:materialized,:axis)
        w,p=constructive_device_fixture()
        if strategy===:axis
            positions=permutedims(Array(p.cells.state.position)[1:2,:])
            axis=Interior.device_cellular_initialize(w,positions;capacity=16,hidden=Int(p.cells.geometry.hidden),
                seed=911,support_strategy=:axis,mechanical_substeps=1,biomass=2.,energy=2.,metabolites=.2,
                initial=(.01,.01,0.),diffusion=(.03,.03,.01),growth_rate=.3,
                speed=0.,repulsion=0.,stiffness=0.,basal_rate=.01,neural_cost=0.,mutation=0.,division_mass=1.7)
            r=ResourcePool(w;potentials=(4,4,0),initial=(0,0,0),diffusion=(0,0,0))
            reference=CellularEcology(w,positions,r;hidden=Int(axis.geometry.hidden),capacity=16)
            feed=ResourcePool(w;names=(:raw_a,:raw_b),potentials=(6,6),initial=(.3,.3),diffusion=(.01,.01))
            habitat=CatalyticHabitat(w,reference,feed;secretion_rate=.2,regulation=false,
                investment_bias=1.,mutation=0.)
            matrix=SorptiveMatrix(w,habitat)
            p=Interior.DeviceConstructiveEcology(w,axis,habitat;matrix)
        end
        begin
            checkpoint=Interior.device_habitat_checkpoint(w,p)
            w2=World(w.grid;seed=911);fill!(w2.fields.nutrient,0)
            candidate=Interior.device_habitat_restore(w2,checkpoint)
            for _ in 1:20
                Interior.device_cellular_step!(w,p.cells,.01;anabolism=p.habitat)
                Interior.device_support_fast_step!(w2,candidate,.01)
                @test Interior.device_habitat_snapshot(p)==Interior.device_habitat_snapshot(candidate)
            end
            p.cells.environment.status[1]=8;candidate.cells.environment.status[1]=8
            before=Interior.device_habitat_snapshot(candidate;check_status=false)
            @test_throws ArgumentError Interior.device_support_fast_step!(w2,candidate,.01)
            @test Interior.device_habitat_snapshot(candidate;check_status=false)==before
        end
    end
end

@testset "Support fastpath carrier composes with launch profiling" begin
    w,p=constructive_device_fixture();before=Interior.device_habitat_snapshot(p)
    ordinary=Interior.device_profile_plan(w,p.cells,.01;anabolism=p.habitat)
    wrapped=Interior.device_profile_plan(w,Interior.DeviceSupportFastCarrier(p.cells),.01;anabolism=p.habitat)
    @test wrapped==ordinary
    @test count(==(:dc_return_waste!),wrapped)==1 && count(==(:dc_metabolism!),wrapped)==1
    trace=Interior.DeviceProfileTrace(Symbol[])
    carrier=Interior.DeviceSupportFastCarrier(Interior.DeviceProfileCarrier(p.cells,trace))
    Interior.device_cellular_enqueue!(w,carrier,.01f0,:projected,p.habitat,UInt32(p.cells.attempted_steps+1))
    @test count(==(:dc_support_fast_return_waste!),trace.labels)==1
    @test count(==(:dc_support_fast_metabolism!),trace.labels)==1
    @test Interior.device_habitat_snapshot(p)==before
end
