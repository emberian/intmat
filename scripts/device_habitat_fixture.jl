# Shared bounded closed-resource fixture; initialization only.
function constructive_device_fixture(backend=CPU();support_strategy=:materialized,energy_type=Float32,inventory_type=Float32)
    w=World(Grid(8,8;extent=(8,8,1));seed=911);fill!(w.fields.nutrient,0)
    r=ResourcePool(w;potentials=(4,4,0),initial=(.01,.01,0),diffusion=(.03,.03,.01),
        source_rates=(0,0,0),source_budgets=(0,0,0),sink_rates=(0,0,0),sink_budgets=(0,0,0))
    c=CellularEcology(w,Float64[2 6;4 4;0 0],r;capacity=16,biomass=2.,energy=2.,division_mass=1.7,
        minimum_mass=.1,growth_rate=.3,mutation=0.,speed=0.,repulsion=0.,stiffness=0.,basal_rate=.01,neural_cost=0.)
    c.metabolites.=.2
    feed=ResourcePool(w;names=(:raw_a,:raw_b),potentials=(6,6),initial=(.3,.3),diffusion=(.01,.01),
        source_rates=(0,0),source_budgets=(0,0),sink_rates=(0,0),sink_budgets=(0,0))
    habitat=Interior.CatalyticHabitat(w,c,feed;secretion_rate=.2,regulation=false,investment_bias=1.,mutation=0.)
    matrix=Interior.SorptiveMatrix(w,habitat)
    target=World(w.grid;backend,seed=911);fill!(target.fields.nutrient,0)
    device=Interior.DeviceCellular(target,c;capacity=16,radius_bound=2.,mechanical_substeps=1,support_strategy,energy_type,inventory_type)
    target,Interior.DeviceConstructiveEcology(target,device,habitat;matrix)
end
