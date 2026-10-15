using Interior,Test,KernelAbstractions,Random
isdefined(Interior,:DeviceResourceReductionCarrier) || Base.include(Interior,joinpath(@__DIR__,"..","studies","device_resource_reduction.jl"))
isdefined(Main,:device_habitat_pilot_fixture) || include(joinpath(@__DIR__,"..","scripts","device_habitat_pilot_fixture.jl"))
function reduce_values(values;status=Int32[0])
    high=copy(values);low=nothing
    while true
        count=cld(size(high,1),128);next=zeros(Float32,count,size(high,2));nextlow=similar(next)
        Interior.dc_resource_reduce!(CPU(),32)(next,nextlow,high,low,status,size(high,1);ndrange=size(next))
        high,low=next,nextlow
        count==1 && return high,low,status
    end
end
@testset "Compensated resource tree retains small grants" begin
    for n in (1,127,128,129,16385)
        values=fill(2f0^-60,n,6);values[1,:].=2f0^60
        high,low,status=reduce_values(values)
        @test status==[0]
        for k in 1:6
            expected=sum(BigFloat,values[:,k])
            @test BigFloat(high[1,k])+BigFloat(low[1,k])==expected
        end
    end
    _,_,status=reduce_values(reshape(Float32[2f0^60,1f0,2f0^-60],3,1))
    @test status==[16]
    for x in (NaN32,Inf32)
        _,_,status=reduce_values(fill(x,129,6))
        @test status==[8]
    end
    _,_,status=reduce_values(fill(floatmax(Float32),129,6))
    @test status==[8]
end
function resource_environment(n;budget=.15f0,field=.1f0,latched=false)
    masks=reshape(Float32[i%3==0 ? 0 : i%5==0 ? .5 : 1 for i in 1:3n],n,3)
    ledger=zeros(Float32,3,6);ledger[:,1].=budget;ledger[:,2].=budget
    e=(fields=fill(field,n,3),source_masks=masks,sink_masks=reverse(masks;dims=1),quota=zeros(Float32,3,2),
        resource_ledger=ledger,resource_ledger_low=zeros(Float32,3,6),status=Int32[latched ? 8 : 0])
    q=(source_rates=(.3f0,0f0,.2f0),sink_rates=(.4f0,.2f0,0f0),source_totals=Tuple(vec(sum(masks;dims=1))),sink_totals=Tuple(vec(sum(masks;dims=1))))
    e,q
end
@testset "Reduced resource exchange conserves actual voxel grants" begin
    for n in (1,129,4096),budget in (.15f0,1f-20),field in (0f0,.1f0,1f20)
        e,q=resource_environment(n;budget,field);baseline=deepcopy(e);before=copy(e.fields)
        for state in (e,baseline)
            Interior.dc_resource_quota!(CPU(),32)(state,q,.1f0;ndrange=3)
        end
        grants=zeros(Float32,n,6)
        Interior.dc_resource_grants!(CPU(),32)(e,q,nothing,.1f0,grants;ndrange=(n,3))
        Interior.dc_resource_exchange!(CPU(),32)(baseline,q,nothing,.1f0;ndrange=(n,3))
        @test e.fields==baseline.fields
        high,low,status=reduce_values(grants;status=e.status)
        @test status==[0]
        Interior.dc_resource_reduced_finish!(CPU(),32)(e,high,low;ndrange=3)
        @test e.status==[0]
        for k in 1:3
            supplied=sum(BigFloat,grants[:,k]);removed=sum(BigFloat,grants[:,k+3])
            @test sum(BigFloat,e.fields[:,k])-sum(BigFloat,before[:,k])==supplied-removed
            @test BigFloat(high[1,k])+BigFloat(low[1,k])==supplied
            @test BigFloat(high[1,k+3])+BigFloat(low[1,k+3])==removed
            remaining=BigFloat(e.resource_ledger[k,1])+BigFloat(e.resource_ledger_low[k,1])
            @test remaining>=0
            @test abs(remaining+supplied-BigFloat(budget))<=BigFloat(eps(Float32))^2*max(BigFloat(budget),BigFloat(1f-40))
            @test all(e.fields[:,k].>=0)
        end
    end
    e,q=resource_environment(129);e.fields[1,1]=NaN32
    grants=zeros(Float32,129,6)
    Interior.dc_resource_grants!(CPU(),32)(e,q,nothing,.1f0,grants;ndrange=(129,3))
    @test e.status==[8]
    e,q=resource_environment(129;latched=true);before=deepcopy(e);fill!(grants,7f0)
    Interior.dc_resource_grants!(CPU(),32)(e,q,nothing,.1f0,grants;ndrange=(129,3))
    @test e==before && all(==(7f0),grants)
end
@testset "Resource tree preserves closed/fed complete coupled trajectories" begin
    for variant in ("closed","finite_sources")
        w,p=device_habitat_pilot_fixture(CPU();size=8,founders=2,capacity=32,support_strategy=:materialized,variant)
        cp=Interior.device_habitat_checkpoint(w,p)
        w2=World(w.grid;seed=911);fill!(w2.fields.nutrient,0)
        candidate=Interior.device_habitat_restore(w2,cp)
        carrier=Interior.DeviceResourceReductionCarrier(candidate.cells)
        for _ in 1:12
            Interior.device_cellular_step!(w,p.cells,.001;anabolism=p.habitat)
            Interior.device_resource_reduction_step!(w2,carrier,.001;anabolism=candidate.habitat)
            a=Interior.device_habitat_snapshot(p);b=Interior.device_habitat_snapshot(candidate)
            @test a.cells.state==b.cells.state
            @test a.feedstock.fields==b.feedstock.fields
            @test a.arrays==b.arrays
            @test a.cells.environment.metadata==b.cells.environment.metadata
            for env in (b.cells.environment,b.feedstock)
                @test all((Float64.(env.resource_ledger[:,1:2])+Float64.(env.resource_ledger_low[:,1:2])).>=0)
            end
        end
        # Scratch is reconstructible: strict underlying state restore must accept
        # compensated totals, then an explicitly recreated study carrier resumes.
        cp2=Interior.device_habitat_checkpoint(w2,candidate)
        w3=World(w.grid;seed=911);fill!(w3.fields.nutrient,0)
        restored=Interior.device_habitat_restore(w3,cp2)
        restored_carrier=Interior.DeviceResourceReductionCarrier(restored.cells)
        @test Interior.device_habitat_snapshot(restored)==Interior.device_habitat_snapshot(candidate)
        for _ in 1:3
            Interior.device_resource_reduction_step!(w2,carrier,.001;anabolism=candidate.habitat)
            Interior.device_resource_reduction_step!(w3,restored_carrier,.001;anabolism=restored.habitat)
            @test Interior.device_habitat_snapshot(restored)==Interior.device_habitat_snapshot(candidate)
        end
    end
end
