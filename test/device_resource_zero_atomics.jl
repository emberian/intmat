using Interior,Test,KernelAbstractions,Random
using Interior: JSON3
isdefined(Interior,:DeviceResourceZeroCarrier) || Base.include(Interior,joinpath(@__DIR__,"..","studies","device_resource_zero_atomics.jl"))
isdefined(Main,:device_habitat_pilot_fixture) || include(joinpath(@__DIR__,"..","scripts","device_habitat_pilot_fixture.jl"))
@testset "Resource zero atomic loads preserve mixed grants and failure checks" begin
    function exchange_case(rates,sinks,masks;invalid=false,latched=false)
        fields=reshape(Float32[.1,.2,.3,.4,.5,.6,.7,.8,.9,.1,.2,.3],4,3)
        ledger=zeros(Float32,3,6);ledger[:,1].=.15;ledger[:,2].=.12
        invalid && (ledger[2,5]=NaN32)
        e=(fields=copy(fields),source_masks=copy(masks),sink_masks=reverse(copy(masks);dims=1),
            quota=ones(Float32,3,2),resource_ledger=ledger,resource_ledger_low=zeros(Float32,3,6),status=Int32[latched ? 8 : 0])
        q=(source_rates=Float32.(rates),sink_rates=Float32.(sinks),source_totals=vec(sum(masks;dims=1)),sink_totals=vec(sum(masks;dims=1)))
        e,q
    end
    for (rates,sinks) in (((0,0,0),(0,0,0)),((.3,.2,0),(0,0,0)),((0,0,0),(.4,.2,0)),((.3,0,.2),(.4,.2,0)))
        masks=Float32[1 0 1;0 1 0;.5 0 1;0 .5 0]
        for invalid in (false,true),latched in (false,true)
            e,q=exchange_case(rates,sinks,masks;invalid,latched);candidate=deepcopy(e)
            # Skip quota reset only for deliberate corrupt accumulator coverage.
            if !invalid
                for state in (e,candidate)
                    Interior.dc_resource_quota!(CPU(),32)(state,q,.1f0;ndrange=3)
                end
            end
            Interior.dc_resource_exchange!(CPU(),32)(e,q,nothing,.1f0;ndrange=(4,3))
            Interior.dc_resource_zero_exchange!(CPU(),32)(candidate,q,nothing,.1f0;ndrange=(4,3))
            @test isequal(e,candidate)
            if !invalid && !latched
                for state in (e,candidate)
                    Interior.dc_resource_finish!(CPU(),32)(state;ndrange=3)
                end
                @test isequal(e,candidate)
                @test all(e.resource_ledger[:,1:2].>=0)
                @test all(e.fields.>=0)
                @test sum(e.resource_ledger[:,3])≈sum(e.resource_ledger[:,5])
                @test sum(e.resource_ledger[:,4])≈sum(e.resource_ledger[:,6])
            end
        end
    end
end
@testset "Resource zero guards preserve closed and fed coupled trajectories" begin
    for variant in ("closed","finite_sources")
        w,p=device_habitat_pilot_fixture(CPU();size=8,founders=2,capacity=32,support_strategy=:materialized,variant)
        checkpoint=Interior.device_habitat_checkpoint(w,p)
        w2=World(w.grid;seed=911);fill!(w2.fields.nutrient,0)
        candidate=Interior.device_habitat_restore(w2,checkpoint)
        for _ in 1:12
            Interior.device_cellular_step!(w,p.cells,.001;anabolism=p.habitat)
            Interior.device_resource_zero_step!(w2,candidate,.001)
            @test Interior.device_habitat_snapshot(p)==Interior.device_habitat_snapshot(candidate)
        end
        @test Interior.device_habitat_audit(Interior.device_habitat_snapshot(p)).finite
        if variant=="finite_sources"
            @test Interior.device_habitat_audit(Interior.device_habitat_snapshot(p)).supplied>0
        end
    end
end
