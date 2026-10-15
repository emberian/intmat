# Explicit paid 3D recipe. ResourcePool.initial below is MASS/voxel.
function device_habitat_pilot_fixture(backend;size=128,founders=512,capacity=8192,
        seed=911,variant="closed",source_budget=100.,source_rate=.5)
    variant in ("closed","finite_sources") || error("unknown habitat variant")
    side=ceil(Int,cbrt(founders));positions=zeros(3,founders)
    rng=MersenneTwister(seed+31)
    for (i,xyz) in enumerate(Iterators.product(0:side-1,0:side-1,0:side-1))
        i>founders && break
        positions[:,i].=4 .* (collect(xyz).+.5)./side .+ (.025/side).*(rand(rng,3).-.5)
    end
    grid=Grid(size,size,size;extent=(4,4,4));volume=prod(grid.spacing)
    cpu=World(grid,CPU(),NamedTuple(),NamedTuple(),0.,seed)
    food=ResourcePool(cpu;potentials=(4,4,0),initial=(.02volume,.02volume,0),
        diffusion=(.001,.001,.001))
    cells=CellularEcology(cpu,positions,food;capacity,seed,biomass=1.6,energy=1.,metabolites=.12,
        density=128.,exchange_radius=.18)
    raw_rates=variant=="closed" ? (0.,0.) : (source_rate/prod(grid.size),source_rate/prod(grid.size))
    raw_budget=variant=="closed" ? (0.,0.) : (source_budget,source_budget)
    raw=ResourcePool(cpu;names=(:raw_a,:raw_b),potentials=(6,6),initial=(4volume,4volume),
        diffusion=(.001,.001),source_rates=raw_rates,source_budgets=raw_budget)
    habitat=CatalyticHabitat(cpu,cells,raw;secretion_rate=.12,secretion_work=.2,
        allocation_exponent=2.,catalytic_rate=6.,half_saturation=.2,decay_rate=.008,
        regulation=true,investment_bias=-1.,seed)
    matrix=SorptiveMatrix(cpu,habitat;capacity_per_enzyme=6.,adsorption_rate=4.,
        desorption_rate=.12,cross_binding=0.,transport_barrier=20.)
    world=World(grid,backend,NamedTuple(),NamedTuple(),0.,seed)
    device=DeviceCellular(world,cells;capacity,support_strategy=:axis,mechanical_substeps=16)
    world,DeviceConstructiveEcology(world,device,habitat;matrix)
end
