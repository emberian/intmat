using Interior, Test, KernelAbstractions
isdefined(Interior,:DeviceCellular) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_cellular.jl"))
isdefined(Interior,:DeviceCatalyticHabitat) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_habitat.jl"))
@kernel function habitat_budget_probe!(out,arg,a,b,mass,energy,growth,growth_work)
    i=@index(Global,Linear)
    values=Interior.device_anabolic_build!(arg,i,a,b,mass,energy,growth,growth_work,4f0,4f0)
    for k in 1:7;out[i,k]=values[k];end
end
function habitat_budget_fixture(C;secretion_work=.2f0)
    arrays=Interior.device_habitat_arrays(CPU(),C,1,1)
    carrier=(geometry=nothing,state=(sites=ones(Int32,C,1),support=ones(Float32,C,1),support_count=ones(Int32,C)),
        environment=(status=zeros(Int32,1),))
    arrays.requests.=.03f0
    argument=Interior.device_habitat_argument(carrier,arrays,(;secretion_work))
    arrays,argument,zeros(Float32,C,7)
end
@testset "Paid device catalyst conservative scalar and atomic transfers" begin
    bits=zeros(UInt32,1,1);status=zeros(Int32,1)
    got=Interior.dh_catalyst_deposit!(bits,status,1,1,.2)
    @test 0<=got<=.2
    @test Float64(reinterpret(Float32,bits[1]))==got
    @test status==[0]
    arrays,arg,out=habitat_budget_fixture(1)
    habitat_budget_probe!(CPU(),64)(out,arg,.1f0,.2f0,1f0,.001f0,.02f0,1f0;ndrange=1)
    KernelAbstractions.synchronize(CPU())
    @test arg.status==[0]
    @test all(x->isfinite(x)&&x>=0,out)
    @test out[1,3]>1 && all(>(0),arrays.last_secretion)
    before=Float64(.1f0)+Float64(.2f0)+1
    after=sum(Float64,out[1,1:3])+sum(Float64,reinterpret(Float32,arrays.enzyme_bits))+
        sum(Float64,arrays.cell_inert)+sum(Float64,arrays.cell_inert_low)
    @test abs(after-before)<1e-10
    before_potential=4*(Float64(.1f0)+Float64(.2f0))+Float64(.001f0)
    after_potential=4sum(Float64,out[1,1:2])+Float64(out[1,4])+sum(Float64,out[1,5:6])+
        sum(Float64,arrays.cell_ledger[:,9:10])+sum(Float64,arrays.cell_ledger_low[:,9:10])
    @test abs(after_potential-before_potential)<1e-10
    # A field too large to accept the request must not debit either donor.
    arrays,arg,out=habitat_budget_fixture(1)
    arrays.enzyme_bits.=reinterpret(UInt32,1f20)
    habitat_budget_probe!(CPU(),64)(out,arg,.1f0,.2f0,1f0,1f0,0f0,1f0;ndrange=1)
    KernelAbstractions.synchronize(CPU())
    @test out[1,1:4]==[.1f0,.2f0,1f0,1f0]
    @test all(iszero,arrays.last_secretion) && all(iszero,arrays.cell_inert)
    # Contending donors share one recipient, but every successful CAS is paid.
    C=256;arrays,arg,out=habitat_budget_fixture(C;secretion_work=0f0)
    habitat_budget_probe!(CPU(),64)(out,arg,.1f0,.2f0,1f0,1f0,0f0,1f0;ndrange=C)
    KernelAbstractions.synchronize(CPU())
    after=sum(Float64,out[:,1:3])+sum(Float64,reinterpret(Float32,arrays.enzyme_bits))+
        sum(Float64,arrays.cell_inert)+sum(Float64,arrays.cell_inert_low)
    @test arg.status==[0]
    @test abs(after-C*(Float64(.1f0)+Float64(.2f0)+1))<1e-8
    @test all(x->isfinite(x)&&x>=0,out)
end
include(joinpath(@__DIR__,"..","scripts","device_habitat_fixture.jl"))
@testset "Resident shared anabolism and consequential sorptive habitat" begin
    w,p=constructive_device_fixture()
    first=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
    for _ in 1:50;advance!(w,p,.01);end
    snapshot=Interior.device_habitat_snapshot(p);last=Interior.device_habitat_audit(snapshot)
    @test last.finite && last.population>first.population
    @test all(>(0),last.secreted) && all(>(0),last.enzyme)
    @test all(>(0),last.solubilized) && all(>(0),last.bound)
    @test minimum(snapshot.arrays.permeability)<1
    @test abs(last.material-first.material-(last.supplied-first.supplied)+(last.removed-first.removed))<1e-4
    @test abs(last.potential-first.potential+last.work-first.work+last.heat-first.heat-
        (last.supplied_potential-first.supplied_potential)+(last.removed_potential-first.removed_potential))<2e-4
    @test snapshot.cells.environment.metadata[6]==50
    checkpoint=Interior.device_habitat_checkpoint(w,p)
    w2=World(w.grid;seed=911);fill!(w2.fields.nutrient,0)
    restored=Interior.device_habitat_restore(w2,checkpoint)
    @test restored.habitat.cells===restored.cells
    @test restored.habitat.argument.status===restored.cells.environment.status
    @test restored.habitat.feedstock.argument.status===restored.cells.environment.status
    @test Interior.device_habitat_snapshot(restored)==snapshot
    for _ in 1:5;advance!(w,p,.01);advance!(w2,restored,.01);end
    @test Interior.device_habitat_snapshot(p)==Interior.device_habitat_snapshot(restored)
    # Duplicate scheduling is rejected before mutation, including owned children.
    @test_throws ArgumentError Interior.validate_composition((p,p.cells))
end
@testset "Catalytic inheritance follows surviving parents through slot reuse" begin
    w,p=constructive_device_fixture();c=p.cells;a=p.habitat.arrays
    # Root1 dies while root2 divides into its vacated slot. Distinct readouts
    # make copying stale dead-slot genotype observable even without regulation.
    a.weights[1,:,:].=11f0;a.weights[2,:,:].=22f0
    c.state.mass[1]=.05f0;c.state.metabolites[1,:].=.2f0;c.state.energy[1]=1f0
    initial=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
    advance!(w,p,.01)
    final=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
    @test c.environment.metadata[4]==1
    @test c.state.ids[1]!=1 && c.state.parents[1]==2
    @test all(==(22f0),a.weights[1,:,:])
    @test all(iszero,a.last_secretion[1,:]) && all(iszero,a.expression[1,:])
    @test sum(a.cell_ledger[1,7:8])+sum(a.cell_ledger_low[1,7:8])>0 # old root paid before death
    @test abs(final.material-initial.material)<1e-4
    @test abs(final.potential-initial.potential+final.work-initial.work+final.heat-initial.heat)<2e-4
end

@testset "Device constructive generic inventory boundary" begin
    w,p=constructive_device_fixture()
    audit=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
    @test Interior.owned_material_mass(p)≈audit.material atol=1e-12
    @test Interior.process_state_finite(p)
    @test !Interior.arena_process_supported(p)
    @test !Interior.arena_process_supported(p.habitat)
end

@testset "Zero barrier remains transparent for finite extreme catalysts" begin
    w,p=constructive_device_fixture();h=p.habitat
    h.arrays.enzyme_bits.=reinterpret(UInt32,floatmax(Float32))
    q=merge(h.parameters,(transport_barrier=0f0,))
    Interior.dh_permeability!(CPU(),64)(h.arrays,p.cells.environment,p.cells.geometry,q;ndrange=length(h.arrays.permeability))
    KernelAbstractions.synchronize(CPU())
    @test all(==(1f0),h.arrays.permeability)
    @test p.cells.environment.status==[0]
end
