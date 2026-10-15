# Initialization and snapshot comparison only; numerical stages remain in Interior.
function constructive_batch_fixture(seed;flux=false,strategy=:axis,budget_scale=1.,backend=CPU(),
        shape=(6,6,1),extent=(6.,6.,1.),population=2,capacity=max(8,2population),hidden=2,mechanics=false)
    w=World(Grid(shape...;extent);seed);fill!(w.fields.nutrient,0)
    r=ResourcePool(w;potentials=(4,4,0),initial=(.01,.01,0),diffusion=(.01,.01,.01),
        source_rates=flux ? (.001,.002,0.) : (0.,0.,0.),source_budgets=flux ? (.5*budget_scale,.7*budget_scale,0.) : (0.,0.,0.))
    nd=shape[3]==1 ? 2 : 3;k=ceil(Int,population^(1/nd));positions=zeros(3,population)
    for i in 1:population
        index=i-1
        positions[1,i]=(mod(index,k)+.5)*extent[1]/k
        positions[2,i]=(mod(div(index,k),k)+.5)*extent[2]/k
        positions[3,i]=nd==3 ? (div(index,k*k)+.5)*extent[3]/k : 0.
    end
    c=CellularEcology(w,positions,r;seed,capacity,hidden,biomass=2.,energy=2.,
        division_mass=1.7,minimum_mass=.1,growth_rate=.3,mutation=.04,mutation_probability=.3,
        speed=mechanics ? .15 : 0.,repulsion=mechanics ? 8. : 0.,stiffness=mechanics ? 3. : 0.,basal_rate=.01)
    c.metabolites.=.2;c.weights[hidden+10,end,:].+=2
    f=ResourcePool(w;names=(:raw_a,:raw_b),potentials=(6,6),initial=(.3,.3),diffusion=(.01,.01),
        source_rates=flux ? (.002,.001) : (0.,0.),source_budgets=flux ? (.7*budget_scale,.5*budget_scale) : (0.,0.))
    h=Interior.CatalyticHabitat(w,c,f;seed=seed+1000,secretion_rate=.2,regulation=true,
        investment_bias=1.,mutation=.03,mutation_probability=.4)
    for enzyme in h.enzymes;fill!(enzyme,.01*prod(w.grid.spacing[d] for d in 1:nd));end
    m=Interior.SorptiveMatrix(w,h)
    target=backend==CPU() ? w : World(w.grid,backend,NamedTuple(),NamedTuple(),0.,Int(seed))
    d=Interior.DeviceCellular(target,c;mechanical_substeps=mechanics ? 16 : 1,radius_bound=2.,support_strategy=strategy)
    target,Interior.DeviceConstructiveEcology(target,d,h;matrix=m)
end
function constructive_snapshot_equal(a,b)
    isequal(a.cells,b.cells)&&isequal(a.arrays,b.arrays)&&isequal(a.feedstock,b.feedstock)&&
        a.parameters==b.parameters&&a.feed_parameters==b.feed_parameters
end
