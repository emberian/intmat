@testset "Semi-implicit conserved phase" begin
    for grid in (Grid(16,12;extent=(16,24,1)),Grid(8,7,6;extent=(8,14,6)))
        w=World(grid); p=ImplicitPhaseSeparation(w;nutrient_diffusion=0,activity_relaxation=0)
        mass=sum(Float64,p.phi); energy=phase_energy(w,p)
        for _ in 1:80
            step!(w,p,0.2)
            next=phase_energy(w,p)
            @test next<=energy+1e-5
            energy=next
        end
        @test sum(Float64,p.phi)≈mass atol=3e-6
        @test hasproperty(record_process_fields(p),:chemical_forcing)
        @test size(phase_optical_fields(w,p).density)==grid.size
        fill!(p.phi,0.3); fill!(p.chemical_forcing,4)
        step!(w,p,0.2)
        @test maximum(abs.(p.phi.-0.3f0))<1e-7
    end
    # Same spatial discretization, temporal refinement against explicit Euler.
    errors=Float64[]
    for dt in (0.002,0.001,0.0005)
        w=World(Grid(12,12;extent=(12,12,1))); v=World(w.grid)
        p=ImplicitPhaseSeparation(w;nutrient_diffusion=0,activity_relaxation=0)
        q=PhaseSeparation(v;nutrient_diffusion=0,activity_relaxation=0)
        for _ in 1:round(Int,0.02/dt); step!(w,p,dt); step!(v,q,dt); end
        push!(errors,maximum(abs.(p.phi.-q.phi)))
    end
    @test errors[2]<0.6errors[1]
    @test errors[3]<0.6errors[2]
    w=World(Grid(12,12;extent=(12,12,1)))
    p=ImplicitPhaseSeparation(w;noise=0,nutrient_diffusion=0,activity_relaxation=0)
    p.chemical_forcing[1,1,1]=1
    step!(w,p,0.1)
    @test sum(Float64,p.phi)≈0 atol=1e-7
    @test maximum(abs,p.phi)>0.001
    w.fields.deposit[1,1,1]=1
    @test_throws ArgumentError step!(w,p,0.1)
    # Uniform pinning retains the exact explicit constitutive mobility.
    v=World(w.grid); fill!(v.fields.deposit,4)
    q=ImplicitPhaseSeparation(v;nutrient_diffusion=0,activity_relaxation=0)
    @test Interior.implicit_mobility(v,q)≈1/9
end

@testset "Implicit phase checkpoint continuity" begin
    w=World(Grid(16,16;extent=(16,16,1)))
    p=ImplicitPhaseSeparation(w;nutrient_diffusion=0,activity_relaxation=0)
    step!(w,p,0.1)
    io=IOBuffer()
    Interior.Serialization.serialize(io,(w,p,p))
    seekstart(io)
    v,q,alias=Interior.Serialization.deserialize(io)
    @test q===alias
    step!(w,p,0.1); step!(v,q,0.1)
    @test p.phi==q.phi
    @test w.time==v.time
    singleton=World(Grid(1,1,1))
    uniform=ImplicitPhaseSeparation(singleton;noise=0,mean=0.3)
    step!(singleton,uniform,0.1)
    @test uniform.phi[1]==0.3f0
end
