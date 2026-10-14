@testset "volume transport" begin
    T,L=transport_segment(1.,0.,2.,0.,0.5)
    @test T ≈ exp(-1)
    @test L == 0
    T,L=transport_segment(1.,0.,2.,3.,0.5)
    @test L ≈ 1.5*(1-exp(-1))
    @test transport_segment(1.,0.,0.,3.,0.5) == (1.,1.5)
    # Subdivision invariance for a homogeneous medium.
    t=1.; l=0.
    for i in 1:100
        t,l=transport_segment(t,l,2.,3.,0.005)
    end
    @test t ≈ T
    @test l ≈ L
    field=reshape(Float32.(1:27),3,3,3)
    @test volume_sample(field,(0.,0.,0.)) == field[2,2,2]
    @test volume_sample(field,(2.,0.,0.)) == 0
    mktempdir() do dir
        constant=ones(Float32,4,4,4)
        settings=RenderSettings(width=8,height=8,steps=8,scattering=0.)
        pixels=render_volume((density=constant,);settings)
        @test size(pixels)==(3,8,8)
        @test all(isfinite,pixels)
        path=render_volume((density=constant,),joinpath(dir,"test.ppm");settings)
        @test isfile(path)
        @test_throws DimensionMismatch export_volume((density=constant,activity=zeros(3,3,3)),joinpath(dir,"bad"))
        manifest=export_volume((density=constant,),dir)
        @test isfile(manifest)
        @test filesize(joinpath(dir,"density.vol")) == 48+4length(constant)
    end
end
