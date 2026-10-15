# Initialization only; all numerical evolution is the native Interior composition.
using TOML
function native_constructed_garden(backend,cfg)
    haskey(cfg,"matrix_mechanics") && error("CPU MatrixDrag is not implemented by this native recipe")
    precision=get(cfg["native"],"energy_type","Float32")
    precision in ("Float32","Float64") || throw(ArgumentError("native.energy_type must be Float32 or Float64"))
    energy_type=precision=="Float64" ? Float64 : Float32
    shape=Tuple(Int.(cfg["grid"]["shape"]));extent=Tuple(Float64.(cfg["grid"]["extent"]))
    length(shape)==3 && all(>(0),shape) && all(>(0),extent) || error("positive3D grid required")
    all(==(shape[1]),shape) || error("observation helper currently requires a cubic grid")
    seed=Int(cfg["seed"]);options=Dict(Symbol(k)=>v for (k,v) in cfg["cellular"])
    n=Int(pop!(options,:founders));capacity=Int(options[:capacity]);0<n<=capacity || error("invalid founder count")
    grid=Grid(shape...;extent);volume=prod(grid.spacing);V=prod(shape)
    cpu=World(grid,CPU(),NamedTuple(),NamedTuple(),0.,seed)
    function pool(section,names)
        q=Dict(Symbol(k)=>v for (k,v) in section)
        q[:initial]=Float64.(q[:initial]).*volume
        for k in (:source_rates,:sink_rates);q[k]=Float64.(q[k])./V;end
        ResourcePool(cpu;names, q...)
    end
    food=pool(cfg["resources"],[:precursor_a,:precursor_b,:waste])
    side=ceil(Int,cbrt(n));xyz=zeros(3,n);rng=MersenneTwister(seed+31)
    sites=shuffle!(rng,collect(Iterators.product(0:side-1,0:side-1,0:side-1)))
    for (i,pos) in enumerate(sites)
        i>n && break
        for d in 1:3;xyz[d,i]=extent[d]*(pos[d]+.5+.025*(rand(rng)-.5))/side;end
    end
    cells=CellularEcology(cpu,xyz,food;seed,options...)
    feedcfg=deepcopy(cfg["feedstock"]);contrast=Float64(pop!(feedcfg,"niche_contrast",0.))
    0<=contrast<1 || error("raw-food niche contrast must lie in[0,1)")
    feedstock=pool(feedcfg,[:raw_a,:raw_b])
    if contrast>0
        for k in 1:2
            sign=k==1 ? 1. : -1.;mask=feedstock.source_masks[k]
            for index in CartesianIndices(mask)
                x=(index[1]-.5)/shape[1];z=(index[3]-.5)/shape[3]
                mask[index]=Float32(1+sign*contrast*sinpi(2x)*sinpi(2z))
            end
            feedstock.source_rates[k]=Float32(feedcfg["source_rates"][k]/sum(Float64,mask))
            feedstock.fields[k].*=mask
        end
    end
    habitat=CatalyticHabitat(cpu,cells,feedstock;seed,Dict(Symbol(k)=>v for (k,v) in cfg["catalysis"])...)
    matrix=SorptiveMatrix(cpu,habitat;Dict(Symbol(k)=>v for (k,v) in cfg["matrix"])...)
    q=cfg["native"];w=World(grid,backend,NamedTuple(),NamedTuple(),0.,seed)
    device=DeviceCellular(w,cells;capacity,support_strategy=Symbol(q["support_strategy"]),
        mechanical_substeps=Int(q["mechanical_substeps"]),radius_bound=q["radius_bound"],bin_width=q["bin_width"],energy_type)
    p=DeviceConstructiveEcology(w,device,habitat;matrix)
    dt=Float32(q["dt"]);0<dt<=Interior.stability_limit(w,p) || error("submitted Float32 dt exceeds coupled limit")
    # Waste and residue are counted too: this deliberately loose bound cannot
    # mistake a desired population ceiling for an ecological carrying capacity.
    initial=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
    funding=sum(food.source_remaining)+sum(feedstock.source_remaining)
    bound=floor(Int,(initial.material+funding)/Float64(device.parameters.minimum_mass))+1
    bound<=capacity || error("capacity $capacity below conservative funded live-count bound $bound")
    w,p,(initial=initial,remaining_source_budget=funding,funded_population_upper_bound=bound,
        nominal_source_duration=minimum(cfg["feedstock"]["source_budgets"][k]/cfg["feedstock"]["source_rates"][k]
            for k in 1:2 if cfg["feedstock"]["source_rates"][k]>0;init=Inf))
end
