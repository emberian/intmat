using Interior, Test, KernelAbstractions
isdefined(@__MODULE__,:constructive_device_fixture) || include(joinpath(@__DIR__,"..","scripts","device_habitat_fixture.jl"))
@testset "Wide resident inventory and precision-aware strict restart" begin
    w,p=constructive_device_fixture(;energy_type=Float64,inventory_type=Float64)
    @test eltype(p.cells.state.mass)===Float64
    @test eltype(p.cells.state.metabolites)===Float64
    @test eltype(p.cells.environment.fields)===Float64
    @test eltype(p.cells.environment.demand)===Float64
    @test eltype(p.cells.environment.returns)===Float64
    @test eltype(p.habitat.feedstock.arrays.fields)===Float64
    @test eltype(p.cells.ownership.stages.resource.grants)===Float64
    @test eltype(p.cells.state.hidden)===Float32
    @test_throws ArgumentError constructive_device_fixture(;inventory_type=Float64)
    first=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
    for _ in 1:1000;advance!(w,p,.01);end
    last=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
    @test last.finite && last.population>first.population
    @test all(>(0),last.secreted) && all(>(0),last.solubilized)
    @test abs(last.material-first.material)<1e-8
    @test abs(last.potential+last.work+last.heat-first.potential)<1e-8
    checkpoint=Interior.device_habitat_checkpoint(w,p)
    rw=World(w.grid;seed=911);fill!(rw.fields.nutrient,0)
    restored=Interior.device_habitat_restore(rw,checkpoint)
    @test Interior.device_habitat_snapshot(restored)==Interior.device_habitat_snapshot(p)
    for _ in 1:7;advance!(w,p,.01);advance!(rw,restored,.01);end
    @test Interior.device_habitat_snapshot(restored)==Interior.device_habitat_snapshot(p)
    damaged=merge(checkpoint,(cells=merge(checkpoint.cells,(snapshot=merge(checkpoint.cells.snapshot,
        (geometry=merge(checkpoint.cells.snapshot.geometry,(inventory_bits=Int32(32),)),)),)),))
    @test_throws ArgumentError Interior.device_habitat_restore(rw,damaged)
end

@testset "Wide exchange budget retains transfers below Float32 ULP" begin
    w,p=constructive_device_fixture(;energy_type=Float64,inventory_type=Float64)
    e=p.cells.environment
    fill!(e.fields,1.);fill!(e.source_masks,1f0);fill!(e.sink_masks,0f0)
    fill!(e.resource_ledger,0f0);fill!(e.resource_ledger_low,0f0)
    e.resource_ledger[:,1].=1f0
    q=merge(p.cells.parameters,(source_rates=(1f-8,1f-8,1f-8),source_totals=(64f0,64f0,64f0),sink_rates=(0f0,0f0,0f0)))
    before=sum(e.fields)
    for _ in 1:50;Interior.device_resource_phase!(p.cells,e,q,p.cells.geometry,.001f0);end
    Interior.device_cellular_status!(p.cells)
    after=sum(e.fields);supplied=sum(Float64,e.resource_ledger[:,3])+sum(Float64,e.resource_ledger_low[:,3])
    @test after>before
    @test abs((after-before)-supplied)<1e-12
end

@testset "Species diffusion conserves wide physical inventory" begin
    residuals=Dict{DataType,Vector{Float64}}()
    for T in (Float32,Float64)
        n=8;V=n^3
        fields=[T(.01k*(1+.4sin(i*.73)+.2cos(i*.29))) for i in 1:V,k in 1:3]
        e=(fields=fields,scratch=similar(fields),permeability=Float32[.2+.8abs(sin(i*.63)) for i in 1:V],status=zeros(Int32,1))
        q=(diffusion=(.06f0,.03f0,.01f0),)
        g=(dims=(Int32(n),Int32(n),Int32(n)),spacing=(.25f0,.25f0,.25f0))
        before=vec(sum(Float64.(fields);dims=1))
        for _ in 1:4000
            Interior.dc_diffuse!(CPU(),128)(e,q,g,.001f0;ndrange=(V,3))
            Interior.dc_commit_resources!(CPU(),128)(e;ndrange=(V,3))
        end
        KernelAbstractions.synchronize(CPU())
        residuals[T]=vec(sum(Float64.(e.fields);dims=1)).-before
        @test all(isfinite,e.fields) && minimum(e.fields)>=0
    end
    @test maximum(abs,residuals[Float64])<1e-11
    @test maximum(abs,residuals[Float64])<maximum(abs,residuals[Float32])/1000
end

@testset "Wide waste scatter closes nonnormalized support weights" begin
    s=(waste=[1.],sites=reshape(Int32[1,2,3],1,3),support=reshape(Float32[.1,.2,.7],1,3),support_count=Int32[3])
    e=(returns=zeros(3),status=zeros(Int32,1))
    Interior.dc_return_waste!(CPU(),32)(s,e,nothing;ndrange=1)
    KernelAbstractions.synchronize(CPU())
    @test sum(Float64,s.support)!=1
    @test sum(e.returns)==1
    @test all(>=(0),e.returns)
end

@testset "Wide multi-level source/sink reductions keep endpoint credits" begin
    V=129^2
    work=Interior.DeviceResourceWorkspace(CPU(),V;inventory_type=Float64)
    @test length(work.levels)==3
    @test all(level->eltype(level.high)===Float64 && eltype(level.low)===Float64,work.levels)
    status=zeros(Int32,1)
    for k in 1:6,i in 1:V;work.grants[i,k]=1e-7*(1+.3sin(i*.71+k));end
    initial=vec(sum(work.grants;dims=1))
    high=work.grants;low=nothing
    for level in work.levels
        Interior.dc_resource_tree_reduce!(CPU(),128)(level.high,level.low,high,low,status,size(high,1);ndrange=size(level.high))
        high=level.high;low=level.low
    end
    KernelAbstractions.synchronize(CPU())
    @test only(status)==0
    @test vec(high.+low)≈initial rtol=1e-14 atol=1e-16
end
