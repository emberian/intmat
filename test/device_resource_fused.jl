using Interior, Test, Random, KernelAbstractions
isdefined(Interior,:DeviceResourceFusedCarrier) || Base.include(Interior,joinpath(@__DIR__,"..","studies","device_resource_fused.jl"))
@testset "Error-free Float32 endpoint pairs" begin
    rng=MersenneTwister(7821)
    setprecision(BigFloat,512) do
        for _ in 1:1000
            a=reinterpret(Float32,rand(rng,UInt32)&0x7f7fffff)
            b=reinterpret(Float32,rand(rng,UInt32)&0x7f7fffff)
            request=reinterpret(Float32,rand(rng,UInt32)&0x7f7fffff)
            hi,lo=Interior.resource_endpoint_difference(a,b)
            @test BigFloat(hi)+BigFloat(lo)==BigFloat(a)-BigFloat(b)
            @test Interior.resource_endpoint_gt(a,b,request)==(BigFloat(a)-BigFloat(b)>BigFloat(request))
            rounded=max(0f0,a-b)
            @test Interior.resource_endpoint_gt(a,b,rounded)==(BigFloat(a)-BigFloat(b)>BigFloat(rounded))
        end
    end
end
function resource_fused_environment(values;source=1f0,sink=0f0,budget=2f0^27)
    fields=Float32.(values);n=size(fields,1)
    e=(fields=fields,source_masks=ones(Float32,n,3),sink_masks=ones(Float32,n,3),
       quota=ones(Float32,3,2),status=zeros(Int32,1),resource_ledger=zeros(Float32,3,6),resource_ledger_low=zeros(Float32,3,6))
    e.resource_ledger[:,1:2].=budget
    q=(source_rates=ntuple(_->Float32(source),3),sink_rates=ntuple(_->Float32(sink),3),
       source_totals=ntuple(_->Float32(n),3),sink_totals=ntuple(_->Float32(n),3))
    e,q
end
function resource_fused_direct!(e,q,dt;quota=true)
    n=size(e.fields,1);work=Interior.DeviceResourceWorkspace(CPU(),n)
    fill!(work.grants,NaN32)
    process=(backend=CPU(),ownership=(stages=(resource=work,),))
    carrier=Interior.DeviceResourceFusedCarrier(process)
    if quota
        Interior.device_resource_phase!(carrier,e,q,NamedTuple(),Float32(dt))
    else
        level=first(work.levels)
        Interior.dc_resource_fused_first!(CPU(),128)(e,q,Float32(dt),level.high,level.low;ndrange=(size(level.high,1),3))
        high=level.high;low=level.low
        for next in Base.tail(work.levels)
            Interior.dc_resource_tree_reduce!(CPU(),128)(next.high,next.low,high,low,e.status,size(high,1);ndrange=size(next.high))
            high=next.high;low=next.low
        end
        Interior.dc_resource_tree_finish!(CPU(),128)(e,high,low;ndrange=3)
    end
    KernelAbstractions.synchronize(CPU())
    e,work
end
@testset "Actual endpoint credits and conservative request clamps" begin
    e,q=resource_fused_environment(ones(Float32,1,3);source=2f0^28)
    before=copy(e.fields);resource_fused_direct!(e,q,.25;quota=false)
    @test e.status==Int32[0]
    @test e.fields==fill(2f0^26,1,3)
    for k in 1:3
        @test BigFloat(e.resource_ledger[k,3])+BigFloat(e.resource_ledger_low[k,3])==BigFloat(2)^26-1
        @test BigFloat(e.resource_ledger[k,1])+BigFloat(e.resource_ledger_low[k,1])==BigFloat(2)^26+1
    end
    # Rounded subtraction formerly accepted a grant above its exact request.
    request=prevfloat(prevfloat(1f0));initial=3f0*2f0^-25
    e,q=resource_fused_environment(fill(initial,1,3);source=request,budget=2f0)
    resource_fused_direct!(e,q,1;quota=false)
    @test all(x->BigFloat(x)-BigFloat(initial)<=BigFloat(request),e.fields)
    @test e.status==Int32[0]
    rng=MersenneTwister(98)
    for (n,source,sink) in ((1,0f0,0f0),(257,.12f0,0f0),(257,0f0,.07f0),(17000,.12f0,.07f0))
        e,q=resource_fused_environment(rand(rng,Float32,n,3);source,sink,budget=4096f0)
        before=copy(e.fields);e,work=resource_fused_direct!(e,q,.03)
        @test e.status==Int32[0]
        @test all(isnan,work.grants)
        for k in 1:3
            actual=sum(BigFloat.(e.fields[:,k]))-sum(BigFloat.(before[:,k]))
            credit=BigFloat(e.resource_ledger[k,3])+BigFloat(e.resource_ledger_low[k,3])-BigFloat(e.resource_ledger[k,4])-BigFloat(e.resource_ledger_low[k,4])
            @test actual==credit
            @test Float64(e.resource_ledger[k,1])+Float64(e.resource_ledger_low[k,1])>=0
            @test Float64(e.resource_ledger[k,2])+Float64(e.resource_ledger_low[k,2])>=0
        end
    end
end
@testset "Workspace traffic and failure boundary" begin
    e,q=resource_fused_environment(zeros(Float32,257,3);source=.2f0,budget=1f-7)
    e,work=resource_fused_direct!(e,q,.25)
    @test all(isnan,work.grants) # V×6 storage was never read or written.
    @test e.status==Int32[0]
    @test all(k->Float64(e.resource_ledger[k,1])+Float64(e.resource_ledger_low[k,1])>=0,1:3)
    e,q=resource_fused_environment(fill(NaN32,1,3))
    resource_fused_direct!(e,q,.25)
    @test e.status[1]&8==8
    e,q=resource_fused_environment(fill(1f30,1,3);source=1f-20)
    before=copy(e.fields);resource_fused_direct!(e,q,.25)
    @test e.fields==before && e.status[1]==0
end
include(joinpath(@__DIR__,"..","scripts","device_habitat_fixture.jl"))
@testset "Fused actual credits in complete paid ecology" begin
    for fed in (false,true),energy_type in (Float32,Float64)
        w,p=constructive_device_fixture(;energy_type)
        cp=Interior.device_habitat_checkpoint(w,p)
        if fed
            feed=deepcopy(cp.feedstock);feed.resource_ledger[1:2,1].=1f0
            cp=merge(cp,(feedstock=feed,feed_parameters=merge(cp.feed_parameters,(source_rates=(.02f0,.01f0,0f0),))))
        end
        wa=World(w.grid);fill!(wa.fields.nutrient,0)
        wb=World(w.grid);fill!(wb.fields.nutrient,0)
        a=Interior.device_habitat_restore(wa,cp);b=Interior.device_habitat_restore(wb,cp)
        carrier=Interior.DeviceResourceFusedCarrier(b.cells)
        fill!(b.cells.ownership.stages.resource.grants,NaN32)
        for _ in 1:20
            Interior.device_cellular_step!(wa,a.cells,.01;anabolism=a.habitat)
            Interior.device_resource_fused_step!(wb,carrier,.01;anabolism=b.habitat)
        end
        xa=Interior.device_habitat_snapshot(a);xb=Interior.device_habitat_snapshot(b)
        @test xa.cells==xb.cells
        @test xa.arrays==xb.arrays
        @test xa.feedstock==xb.feedstock
        @test all(isnan,b.cells.ownership.stages.resource.grants)
        saved=Interior.device_resource_fused_checkpoint(wb,b,carrier)
        wr=World(w.grid);fill!(wr.fields.nutrient,0)
        restored,restored_carrier=Interior.device_resource_fused_restore(wr,saved)
        @test Interior.device_habitat_snapshot(restored)==xb
        for _ in 1:3
            Interior.device_resource_fused_step!(wb,carrier,.01;anabolism=b.habitat)
            Interior.device_resource_fused_step!(wr,restored_carrier,.01;anabolism=restored.habitat)
        end
        @test Interior.device_habitat_snapshot(restored)==Interior.device_habitat_snapshot(b)
        @test_throws ArgumentError Interior.device_resource_fused_restore(wr,merge(saved,(source=(bad=true,),)))
        fill!(b.cells.environment.status,Int32(8))
        before=Interior.device_habitat_snapshot(b;check_status=false)
        Interior.device_resource_fused_step!(wb,carrier,.01;anabolism=b.habitat,check_status=false)
        after=Interior.device_habitat_snapshot(b;check_status=false)
        @test before.cells.state==after.cells.state && before.cells.environment==after.cells.environment
        @test before.arrays==after.arrays && before.feedstock==after.feedstock
        @test_throws ArgumentError Interior.device_resource_fused_checkpoint(wb,b,carrier)
    end
end
