using Interior, Test, KernelAbstractions
isdefined(Interior,:DeviceProfileCarrier) || Base.include(Interior,joinpath(@__DIR__,"..","studies","device_launch_profile.jl"))
isdefined(Main,:device_habitat_3d_fixture) || include(joinpath(@__DIR__,"..","scripts","device_habitat_gate.jl"))
@testset "Launch profiler observes the shared schedule without changing its law" begin
    for coupled in (false,true)
        w,p=device_habitat_3d_fixture(CPU();size=8)
        wr,reference=device_habitat_3d_fixture(CPU();size=8)
        hook=coupled ? p.habitat : nothing;reference_hook=coupled ? reference.habitat : nothing
        before=Interior.device_habitat_snapshot(p)
        plan=Interior.device_profile_plan(w,p.cells,.001;anabolism=hook)
        @test length(plan)==(coupled ? 122 : 112)
        @test Interior.device_habitat_snapshot(p)==before
        @test p.cells.attempted_steps==0 && w.time==0
        recorded=Int[];waits=Ref(0);count=2length(plan)
        runtime=(HIP=(record=(event)->push!(recorded,event),
            synchronize=(_)->(waits[]+=1),elapsed=(a,b)->Float64(b-a)*1e-6),)
        events=[(2i-1,2i) for i in 1:count]
        observer=Interior.DeviceProfileEvents(plan,events,runtime;steps=2)
        for _ in 1:2
            Interior.device_cellular_step!(wr,reference.cells,.001;anabolism=reference_hook)
            Interior.device_profile_step!(w,p.cells,observer,.001;anabolism=hook)
        end
        Interior.device_cellular_status!(p.cells)
        @test Interior.device_habitat_snapshot(p)==Interior.device_habitat_snapshot(reference)
        Interior.device_profile_wait!(observer)
        summary=Interior.device_profile_finish(observer)
        @test waits[]==1
        @test recorded==collect(1:2count)
        @test length(summary.samples)==count
        @test summary.kernels[:dc_support!][1]==4
        @test summary.categories[:support_sensing_gradient_gather][1]==2
        @test summary.categories[:support_uptake_scatter][1]==2
        @test summary.categories[:dense_diffusion][1]==(coupled ? 4 : 2)
        @test all(sample->sample.seconds==1e-6,summary.samples)
        @test_throws ArgumentError Interior.DeviceProfileEvents(plan,[(1,2)],runtime;steps=2)
        mismatched=Interior.DeviceProfileEvents([:wrong],[(1,2)],runtime)
        @test_throws ArgumentError Interior.device_profile_launch!(mismatched,p.cells,Interior.dc_support!,
            p.cells.state,p.cells.environment,p.cells.parameters,p.cells.geometry;ndrange=Int(p.cells.geometry.capacity))
        @test mismatched.cursor==0
    end
end

@testset "Profile carrier preserves a failed coupled step and its clock" begin
    w,p=device_habitat_3d_fixture(CPU();size=8)
    wr,reference=device_habitat_3d_fixture(CPU();size=8)
    plan=Interior.device_profile_plan(w,p.cells,.001;anabolism=p.habitat)
    runtime=(HIP=(record=(_)->nothing,synchronize=(_)->nothing,elapsed=(a,b)->1e-6),)
    observer=Interior.DeviceProfileEvents(plan,[(2i-1,2i) for i in eachindex(plan)],runtime)
    p.cells.environment.status[1]=8;reference.cells.environment.status[1]=8
    Interior.device_profile_step!(w,p.cells,observer,.001;anabolism=p.habitat)
    Interior.device_cellular_step!(wr,reference.cells,.001;check_status=false,anabolism=reference.habitat)
    @test Interior.device_habitat_snapshot(p;check_status=false)==Interior.device_habitat_snapshot(reference;check_status=false)
    @test p.cells.environment.metadata[6]==0 && p.cells.attempted_steps==1 && w.time==0
    @test length(Interior.device_profile_finish(observer).samples)==122
    @test_throws ArgumentError Interior.device_habitat_checkpoint(w,p)
end
