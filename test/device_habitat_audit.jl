using Interior, Test, KernelAbstractions
isdefined(Interior,:DeviceCatalyticHabitat) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_habitat.jl"))
isdefined(Main,:constructive_device_fixture) || include(joinpath(@__DIR__,"..","scripts","device_habitat_fixture.jl"))

@testset "Habitat strict checkpoint rejects corrupted owned roles" begin
    w,p=constructive_device_fixture()
    original=Interior.device_habitat_checkpoint(w,p)
    for (group,role,index,value) in ((:feedstock,:resource_ledger,(1,1),NaN32),
            (:feedstock,:resource_ledger,(1,1),-1f0),
            (:feedstock,:source_masks,(1,1),-1f0),
            (:arrays,:expression,(1,1),NaN32),
            (:feedstock,:permeability,(1,),2f0))
        payload=deepcopy(original)
        getproperty(getproperty(payload,group),role)[index...]=value
        target=World(w.grid;seed=911);fill!(target.fields.nutrient,0)
        @test_throws ArgumentError Interior.device_habitat_restore(target,payload)
    end
    badtotals=merge(original,(feed_parameters=merge(original.feed_parameters,
        (source_totals=(original.feed_parameters.source_totals[1]/2,
            original.feed_parameters.source_totals[2],0f0),)),))
    target=World(w.grid;seed=911);fill!(target.fields.nutrient,0)
    @test_throws ArgumentError Interior.device_habitat_restore(target,badtotals)
    third=deepcopy(original);third.feedstock.fields[1,3]=1f0
    @test_throws ArgumentError Interior.device_habitat_restore(target,third)
    unsupported=deepcopy(original);unsupported.arrays.bound[1,1]=1f0
    unsupported.arrays.enzyme_bits[1,:].=reinterpret(UInt32,0f0)
    @test_throws ArgumentError Interior.device_habitat_restore(target,unsupported)
    negative=deepcopy(original)
    negative.feedstock.resource_ledger[1,1]=1f0
    negative.feedstock.resource_ledger_low[1,1]=-2f0
    @test_throws ArgumentError Interior.device_habitat_restore(target,negative)
    base_unnormalized=deepcopy(original)
    base_unnormalized.cells.snapshot.environment.resource_ledger[1,1]=1f0
    base_unnormalized.cells.snapshot.environment.resource_ledger_low[1,1]=-.5f0
    @test_throws ArgumentError Interior.device_habitat_restore(target,base_unnormalized)
    unnormalized=deepcopy(original)
    unnormalized.feedstock.resource_ledger[1,1]=1f0
    unnormalized.feedstock.resource_ledger_low[1,1]=-.5f0
    @test_throws ArgumentError Interior.device_habitat_restore(target,unnormalized)
    compensated=deepcopy(original)
    compensated.feedstock.resource_ledger[1,1]=1f0
    compensated.feedstock.resource_ledger_low[1,1]=-eps(1f0)/4
    restored=Interior.device_habitat_restore(target,compensated)
    @test Float64(restored.habitat.feedstock.arrays.resource_ledger[1,1])+Float64(restored.habitat.feedstock.arrays.resource_ledger_low[1,1])==1-Float64(eps(1f0))/4
end

@testset "Habitat generic failed preflight precedes sibling mutation" begin
    w,p=constructive_device_fixture()
    p.cells.environment.status[1]=8
    fill!(w.fields.activity,1);fill!(w.fields.nutrient,1)
    before=deepcopy(w.fields)
    @test_throws ArgumentError step!(w,(Deposition(threshold=0),p),.01)
    @test w.fields==before && w.time==0
end

@testset "Imported scaffold collapse preserves supported inventory and potential" begin
    a=Interior.device_habitat_arrays(CPU(),1,2,1)
    e=(fields=Float32[1 0 0;1f20 0 0],status=zeros(Int32,1))
    a.bound[:,1].=1f0
    q=(capacity_per_enzyme=0f0,cross_binding=0f0,food_potentials=(4f0,4f0,0f0))
    before=[BigFloat(e.fields[i,1])+BigFloat(a.bound[i,1]) for i in 1:2]
    Interior.dh_initial_collapse!(CPU(),32)(a,e,q;ndrange=2)
    KernelAbstractions.synchronize(CPU())
    @test e.status==[0] && all(iszero,a.bound)
    @test e.fields[:,1]==Float32[2,1f20]
    @test all(i->BigFloat(e.fields[i,1])+BigFloat(a.grid_inert[i])+BigFloat(a.grid_inert_low[i])==before[i],1:2)
    @test Float64(a.grid_ledger[1,9])+Float64(a.grid_ledger_low[1,9])==1
    @test Float64(a.grid_inert[2])+Float64(a.grid_inert_low[2])==1
    @test Float64(a.grid_ledger[2,11])+Float64(a.grid_ledger_low[2,11])==4
end
