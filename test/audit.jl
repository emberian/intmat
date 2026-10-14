@testset "composition preflight" begin
    w=World(Grid(4,4)); fill!(w.fields.activity,1)
    before=deepcopy(w.fields)
    @test_throws ArgumentError step!(w,(Deposition(),Diffusion(:missing,0.01)),0.01)
    @test w.fields==before
    @test w.time==0
    p=ParticleCarrier(reshape(Float32[1,1,0],3,1))
    p.hidden=zeros(Float32,3,1)
    @test_throws DimensionMismatch step!(w,(Deposition(),p),0.01)
    @test w.fields==before
    p=ParticleCarrier(reshape(Float64[1,1,0],3,1))
    @test_throws ArgumentError step!(w,(Deposition(),p),0.01)
    @test w.fields==before
    @test_throws ArgumentError step!(w,(),1e-100)
    @test_throws ArgumentError step!(w,(),1e100)
    @test w.time==0
end
@testset "reaction timestep guard" begin
    w=World(Grid(4,4)); seed_activity!(w)
    before=deepcopy(w.fields)
    @test_throws ArgumentError step!(w,(Deposition(),Excitable(gain=1f8)),0.01)
    @test w.fields==before
    fill!(w.fields.activity,1f4)
    @test_throws ArgumentError step!(w,Excitable(),0.01)
    fill!(w.fields.activity,NaN32)
    @test_throws ArgumentError step!(w,Excitable(),0.001)
    choir=build_scene(Dict("simulation"=>Dict("resolution"=>40),"particles"=>Dict("count"=>0),
        "excitable"=>Dict("diffusion"=>0.003,"gain"=>4,"threshold"=>0.09,"material_feedback"=>0.25)))
    @test stability_limit(choir.world,choir.processes[2])>=0.02
end
