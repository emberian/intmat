@testset "Conserved phase material" begin
    w=World(Grid(16,16;extent=(16,16,1));seed=83)
    p=PhaseSeparation(w;activity_relaxation=0,nutrient_diffusion=0)
    clone=PhaseSeparation(w;activity_relaxation=0,nutrient_diffusion=0)
    @test p.phi==clone.phi
    mass=sum(Float64,p.phi); energy=phase_energy(w,p)
    dt=min(0.005,stability_limit(w,p))
    for _ in 1:200
        step!(w,p,dt)
        next=phase_energy(w,p)
        @test next<=energy+1e-6
        energy=next
    end
    @test sum(Float64,p.phi)≈mass atol=2e-6
    @test_throws ArgumentError advance!(w,p,2stability_limit(w,p))
    @test_throws DimensionMismatch Interior.validate(World(Grid(8,8)),p)
    # A wavelength-eight spinodal mode must grow, while reducing free energy.
    ws=World(Grid(16,16;extent=(16,16,1)))
    ps=PhaseSeparation(ws;noise=0,activity_relaxation=0,nutrient_diffusion=0)
    for j in 1:16,i in 1:16; ps.phi[i,j,1]=0.03f0*cos(2f0*Float32(pi)*i/8f0); end
    amplitude=maximum(ps.phi); initial_energy=phase_energy(ws,ps)
    for _ in 1:300; step!(ws,ps,0.005); end
    @test maximum(ps.phi)>1.3amplitude
    @test phase_energy(ws,ps)<initial_energy
    # Uniform compositions are stationary, including 3D and unequal grid spacings.
    w3=World(Grid(8,7,6;extent=(8,14,6)))
    uniform=PhaseSeparation(w3;mean=0.3,noise=0)
    before=copy(uniform.phi)
    step!(w3,uniform,stability_limit(w3,uniform)/2)
    @test uniform.phi==before
    # Mineral mobility changes the evolution; nutrient feedback is conservative.
    wa=World(Grid(12,12;extent=(12,12,1))); wb=World(wa.grid)
    pa=PhaseSeparation(wa); pb=PhaseSeparation(wb)
    fill!(wb.fields.deposit,4)
    wa.fields.nutrient[1,1,1]=2; wb.fields.nutrient[1,1,1]=2
    food=sum(Float64,wa.fields.nutrient)
    step!(wa,pa,0.005); step!(wb,pb,0.005)
    @test maximum(abs.(pa.phi.-pb.phi))>1e-5
    @test sum(Float64,wa.fields.nutrient)≈food atol=2e-6
    @test maximum(wa.fields.activity)>0
    @test hasproperty(Interior.record_process_fields(pa),:phi)
    # Translation equivariance of the full periodic constitutive update.
    wc=World(wa.grid); pc=PhaseSeparation(wc)
    wd=World(wa.grid); pd=PhaseSeparation(wd)
    copyto!(pd.phi,circshift(pc.phi,(3,2,0)))
    step!(wc,pc,0.005); step!(wd,pd,0.005)
    @test pd.phi≈circshift(pc.phi,(3,2,0)) atol=1e-7
end
