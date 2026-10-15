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

@testset "Real FFT phase matches the full complex representation" begin
    # Independent old representation: all Fourier modes, including redundant
    # conjugates. Odd/singleton first axes exercise real-transform layout.
    for grid in (Grid(9,7,5;extent=(9,14,7)),
                 Grid(1,9,7;extent=(2,9,14)),
                 Grid(10,7;extent=(8,14,1)))
        w=World(grid;seed=63); fill!(w.fields.deposit,2)
        p=ImplicitPhaseSeparation(w;mineral_pinning=.7,stabilization=3,
            nutrient_diffusion=0,activity_relaxation=0)
        for I in CartesianIndices(p.chemical_forcing)
            p.chemical_forcing[I]=.15sin(sum(Tuple(I)))
        end
        phi=copy(p.phi); q=zeros(Float32,grid.size)
        for I in CartesianIndices(q)
            q[I]=sum(4sinpi((I[d]-1)/grid.size[d])^2/grid.spacing[d]^2 for d in 1:3 if grid.size[d]>1;init=0.0)
        end
        dt=.15f0; mobility=Interior.implicit_mobility(w,p)
        for _ in 1:30
            u=Float64.(phi)
            old=Interior.FFTW.fft(u)
            nonlinear=Interior.FFTW.fft(p.bulk .* (u.^3 .- u) .+ p.chemical_forcing)
            for I in eachindex(old)
                a=Float64(dt)*mobility*Float64(q[I])
                old[I]=((1+a*p.stabilization)*old[I]-a*nonlinear[I])/(1+a*(p.stabilization+p.kappa*Float64(q[I])))
            end
            phi.=real.(Interior.FFTW.ifft(old))
            step!(w,p,dt)
        end
        @test p.phi≈phi atol=2e-7 rtol=2e-6
        @test size(p.spectrum)==(grid.size[1]÷2+1,grid.size[2],grid.size[3])
        @test size(p.real_buffer)==grid.size
    end
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
