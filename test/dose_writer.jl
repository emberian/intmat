using Test
isdefined(Interior,:DoseThermalWriter) || Base.include(Interior,joinpath(@__DIR__,"../src/dose_writer.jl"))
@testset "dose-owned fixed material writer" begin
    @test Interior.DoseThermalWriter(World(Grid(8,8,8))) isa Interior.DoseThermalWriter
    for pulse in (0f0,7f0)
        full=World(Grid(24,24,24));compact=World(full.grid,CPU(),NamedTuple(),NamedTuple(),0.,full.seed)
        original=ThermalWriter(full;path=:helix,width=.09,period=.2,stop=.2,pulse_rate=pulse,heat=40,conversion=.37)
        dose=Interior.DoseThermalWriter(compact;path=:helix,width=.09,period=.2,stop=.2,pulse_rate=pulse,heat=40,conversion=.37)
        for _ in 1:300;step!(full,original,.001);step!(compact,dose,.001);end
        @test original.temperature==dose.writer.temperature
        @test original.dose==dose.writer.dose
        n=similar(full.fields.nutrient);d=similar(n);Interior.writer_material!(n,d,dose)
        @test n ≈ full.fields.nutrient atol=5f-6 rtol=5f-6
        @test d ≈ full.fields.deposit atol=5f-6 rtol=5f-6
        @test all(x->x>=0,n) && all(x->x>=0,d)
        @test maximum(abs.(n.+d.-1))<=2eps(Float32)
        @test Interior.owned_material_mass(dose)==length(n)
        @test compact.time==full.time
    end
end
@testset "private reservoir and bounded readout" begin
    w=World(Grid(12,12,12));p=Interior.DoseThermalWriter(w;initial_nutrient=.7,initial_deposit=.2,path=:stationary,width=.2,heat=30,conversion=.37)
    before=copy(w.fields.nutrient);solid=copy(w.fields.deposit)
    for _ in 1:50;step!(w,p,.01);end
    @test w.fields.nutrient==before && w.fields.deposit==solid
    @test Interior.owned_material_mass(p) ≈ (.7f0+.2f0)*length(before) rtol=1f-6
    n=zeros(Float32,12,12,2);d=similar(n)
    Interior.writer_material!(n,d,p;dose=view(p.writer.dose,:,:,1:2))
    @test maximum(abs.(n.+d.-.9f0))<=2eps(Float32)
    @test all(isfinite,n) && all(isfinite,d)
    @test_throws ArgumentError step!(w,(p,p.writer),.01)
    @test_throws ArgumentError Interior.writer_material!(p.writer.dose,similar(p.writer.dose),p)
end
@testset "saturation and long-run rounding" begin
    full=World(Grid(8,8,8));private=World(full.grid,CPU(),NamedTuple(),NamedTuple(),0.,full.seed)
    original=ThermalWriter(full;path=:stationary,width=.15,heat=100,conversion=17,stop=2)
    compact=Interior.DoseThermalWriter(private;path=:stationary,width=.15,heat=100,conversion=17,stop=2)
    for _ in 1:2000;step!(full,original,.002);step!(private,compact,.002);end
    n=similar(full.fields.nutrient);d=similar(n);Interior.writer_material!(n,d,compact)
    nutrient_error=maximum(abs.(n.-full.fields.nutrient));deposit_error=maximum(abs.(d.-full.fields.deposit))
    @info "2000-step material rounding" nutrient_error deposit_error nutrient_minimum=minimum(n)
    @test original.temperature==compact.writer.temperature && original.dose==compact.writer.dose
    @test nutrient_error<5f-5 && deposit_error<5f-5
    @test minimum(n)<eps(Float32)
    @test maximum(abs.(n.+d.-1))<=2eps(Float32)
end
@testset "ping-pong ownership and borrowed handles" begin
    w=World(Grid(8,8,8),CPU(),NamedTuple(),NamedTuple(),0.,1)
    p=DoseThermalWriter(w);first=p.writer;second=p.alternate;dose=p.writer.dose
    step!(w,p,.01)
    @test p.writer===second && p.alternate===first
    @test p.writer.dose===dose
    step!(w,p,.01)
    @test p.writer===first && p.alternate===second
    @test_throws ArgumentError step!(w,(p,second),.01)
    for field in (p.writer.temperature,p.writer.scratch,view(p.writer.temperature,:,:,:))
        @test_throws ArgumentError writer_material!(field,similar(p.writer.dose),p)
    end
end
