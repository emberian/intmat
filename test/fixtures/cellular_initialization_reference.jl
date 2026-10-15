# Frozen pre-extraction founder constructor; exact allocation/RNG oracle.
function cellular_initialization_reference(w::World,positions::AbstractMatrix,resources::ResourcePool;
        medium=nothing,seed=w.seed,capacity=4096,hidden=8,biomass=1.,energy=1.,genome_storage=:dense,
        metabolites=.12,exchange_radius=1.,density=1.,store_capacity=.8,energy_capacity=3.,
        uptake_rate=1.2,half_saturation=.2,enzyme_exponent=1.4,
        growth_rate=.8,catabolic_rate=1.4,growth_work=.4,efficiency=.65,
        basal_rate=.035,neural_cost=.0002,starvation_rate=.16,
        minimum_mass=.25,division_mass=1.8,speed=.7,motor_cost=.5,
        junction_rate=.8,repulsion=8.,stiffness=3.,break_strain=.8,max_degree=6,
        response_time=.8,adaptation_time=6.,adaptation_strength=.3,
        plasticity_rate=.002,plasticity_decay=.02,mutation=.08,
        mutation_probability=.03,developmental_partition=.03,packing_barrier=.6)
    n=size(positions,2); h=Int(hidden)
    size(positions,1)==3 && all(isfinite,positions) || throw(ArgumentError("physical positions must be finite 3×N"))
    0<=n<=capacity && h>=2 || throw(ArgumentError("invalid cellular capacity/hidden size"))
    length(resources.fields)==3 && resources.potentials[3]==0 || throw(ArgumentError("cells require two energetic precursors and one zero-potential waste species"))
    pars=(;exchange_radius=Float64(exchange_radius),density=Float64(density),store_capacity=Float64(store_capacity),energy_capacity=Float64(energy_capacity),
        uptake_rate=Float64(uptake_rate),half_saturation=Float64(half_saturation),enzyme_exponent=Float64(enzyme_exponent),
        growth_rate=Float64(growth_rate),catabolic_rate=Float64(catabolic_rate),growth_work=Float64(growth_work),efficiency=Float64(efficiency),
        basal_rate=Float64(basal_rate),neural_cost=Float64(neural_cost),starvation_rate=Float64(starvation_rate),
        minimum_mass=Float64(minimum_mass),division_mass=Float64(division_mass),speed=Float64(speed),motor_cost=Float64(motor_cost),
        junction_rate=Float64(junction_rate),repulsion=Float64(repulsion),stiffness=Float64(stiffness),break_strain=Float64(break_strain),max_degree=Int(max_degree),
        response_time=Float64(response_time),adaptation_time=Float64(adaptation_time),adaptation_strength=Float64(adaptation_strength),
        plasticity_rate=Float64(plasticity_rate),plasticity_decay=Float64(plasticity_decay),mutation=Float64(mutation),
        mutation_probability=Float64(mutation_probability),developmental_partition=Float64(developmental_partition),packing_barrier=Float64(packing_barrier))
    all(x->isfinite(x)&&x>=0,values(pars)) && all(>(0),(density,store_capacity,energy_capacity,exchange_radius,half_saturation,response_time,adaptation_time,motor_cost,minimum_mass,division_mass)) || throw(ArgumentError("invalid cellular parameters"))
    0<=efficiency<=1 && 0<=mutation_probability<=1 && minimum_mass<division_mass/2 || throw(ArgumentError("invalid yield, mutation probability or size scales"))
    all(x->isfinite(x)&&x>=0,(biomass,energy,metabolites)) && biomass>=minimum_mass || throw(ArgumentError("invalid initial cell inventory"))
    2metabolites<=store_capacity*biomass || throw(ArgumentError("initial metabolites exceed cell capacity"))
    energy<=energy_capacity*biomass || throw(ArgumentError("initial energy exceeds cell capacity"))
    rng=MersenneTwister(seed)
    weights=.12f0.*randn(rng,Float32,h+CELLULAR_ACTIONS,2h+CELLULAR_SENSORS,n)
    # A recurrent developmental substrate with a slow inhibitory adaptation.
    # This prior makes temporal state available; it does not prescribe morphology.
    for i in 1:n,j in 1:h
        weights[j,j,i]+=1.2f0
        weights[j,mod1(j+1,h),i]+=isodd(j) ? -.7f0 : .7f0
    end
    for d in 1:3
        weights[h+d,2h+2+d,:].+=.7f0
        weights[h+d,2h+5+d,:].+=.7f0
    end
    weights[h+8,end,:].+=.5f0 # adhesion is available, never compulsory
    weights[h+9,end,:].+=.3f0
    weights[h+10,end,:].+=.5f0
    xyz=Float64.(positions)
    for d in 1:3; xyz[d,:].=d==3 && w.grid.size[3]==1 ? 0. : mod.(xyz[d,:],Float64(w.grid.extent[d])); end
    p=CellularEcology(collect(1:n),zeros(Int,n),collect(1:n),zeros(Int,n),xyz,
        fill(Float64(biomass),n),fill(Float64(metabolites),2,n),fill(Float64(energy),n),
        zeros(Float32,h,n),zeros(Float32,h,n),weights,zeros(Float32,size(weights)),
        zeros(Float32,CELLULAR_ACTIONS,n),zeros(n),resources,medium,CellBonds(),
        CellStencil(w.grid,xyz,exchange_radius),Int(seed),n+1,Int(capacity),0,0,0,
        0.,0.,n*Float64(energy),zeros(2),0.,0.,0.,0.,0.,nothing,pars)
    storage=Symbol(genome_storage)
    storage in (:dense,:indexed) || throw(ArgumentError("genome_storage must be dense or indexed"))
    storage===:indexed && _cellular_enable_indexed_storage!(p)
    validate(w,p)
    p
end
