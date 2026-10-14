@testset "moving particle carrier" begin
    nutrient = fill(0.02f0,8,7,1)
    deposit = zeros(Float32,size(nutrient))
    positions = Float32[1.2 1.2 5.5; 2.3 2.3 4.5; 0 0 0]
    hidden = Float32[0.3 -0.2 0.1; -0.4 0.5 0.2; 0.1 0.2 -0.3; 0.2 -0.1 0.4]
    p = ParticleCarrier(positions;ids=[41,7,100],hidden=hidden)
    original_hidden = copy(p.hidden)
    original_positions = copy(p.positions)
    budget = sum(nutrient)+sum(deposit)+total_particle_nutrient(p)
    for _ in 1:100
        particle_step!(p,nutrient,deposit,0.1f0)
    end
    @test sum(nutrient)+sum(deposit)+total_particle_nutrient(p) ≈ budget rtol=2f-6
    @test minimum(nutrient) >= 0
    @test minimum(p.storage) >= 0
    @test sum(deposit) > 0
    @test p.hidden != original_hidden
    @test p.positions != original_positions
    @test all(iszero,p.positions[3,:])

    a = ParticleCarrier(positions;ids=[41,7,100],hidden=hidden)
    b = deepcopy(a)
    reorder_particles!(b,[3,1,2])
    na,nb = fill(0.1f0,8,7,1),fill(0.1f0,8,7,1)
    da,db = zeros(Float32,8,7,1),zeros(Float32,8,7,1)
    for _ in 1:20
        particle_step!(a,na,da,0.2f0)
        particle_step!(b,nb,db,0.2f0)
    end
    ia,ib = sortperm(a.ids),sortperm(b.ids)
    @test a.ids[ia] == b.ids[ib]
    @test a.hidden[:,ia] == b.hidden[:,ib]
    @test a.positions[:,ia] == b.positions[:,ib]
    @test a.storage[ia] == b.storage[ib]
    @test na == nb
    @test da == db
    @test_throws ArgumentError reorder_particles!(a,[1,1,3])

    p3 = ParticleCarrier(reshape(Float32[2,2,2],3,1);hidden=fill(0.3f0,4,1))
    n3,d3 = ones(Float32,4,4,4),zeros(Float32,4,4,4)
    particle_step!(p3,n3,d3,0.2f0)
    @test p3.positions[3,1] != 2f0
    @test all(0 .<= p3.positions .< 4)
    pzero = deepcopy(p3)
    particle_step!(p3,n3,d3,0)
    @test p3.positions == pzero.positions
    @test p3.hidden == pzero.hidden
    @test_throws ArgumentError particle_step!(p3,n3,d3,-1)

    # Two constitutions see the same nutrient but carry distinct latent memory.
    still = ParticleCarrier(positions;hidden=zeros(Float32,4,3))
    active = ParticleCarrier(positions;hidden=hidden)
    ns,nt = ones(Float32,8,7,1),ones(Float32,8,7,1)
    ds,dt = zeros(Float32,8,7,1),zeros(Float32,8,7,1)
    particle_step!(still,ns,ds,0.2f0)
    particle_step!(active,nt,dt,0.2f0)
    @test still.positions != active.positions
    @test ds != dt

    # Equal physical speed across spatial refinement with identical observations.
    coarse = ParticleCarrier(reshape(Float32[1,1,0],3,1);hidden=fill(0.3f0,4,1))
    fine = ParticleCarrier(reshape(Float32[2,2,0],3,1);hidden=fill(0.3f0,4,1))
    nc,dc = ones(Float32,4,4,1),zeros(Float32,4,4,1)
    nf,df = ones(Float32,8,8,1),zeros(Float32,8,8,1)
    particle_step!(coarse,nc,dc,0.1f0;spacing=(0.5,0.5,1))
    particle_step!(fine,nf,df,0.1f0;spacing=(0.25,0.25,1))
    @test coarse.positions .* 0.5f0 ≈ fine.positions .* 0.25f0

    passive = ParticleCarrier(reshape(Float32[1,1,0],3,1);rule=ParticleRule(speed=0))
    velocity = (fill(0.5f0,4,4,1),zeros(Float32,4,4,1),zeros(Float32,4,4,1))
    particle_step!(passive,nc,dc,0.2f0;spacing=(0.5,0.5,1),velocity=velocity)
    @test passive.positions[1,1] ≈ 1.2f0

    world = World(Grid(8,8))
    attached = ParticleCarrier(reshape(Float32[1,2,0],3,1))
    initial_mass = material_mass(world)+total_particle_nutrient(attached)
    step!(world,(attached,),0.1f0)
    @test material_mass(world)+total_particle_nutrient(attached) ≈ initial_mass
    @test world.time ≈ 0.1f0
end
