"""Initialize device cells directly from positions and uniform resource concentrations.

This optional factory uses the canonical CPU parameter and seeded genome helpers,
but constructs no CPU spatial stencil. Only founder-sized biological arrays are
staged on the host. Resource concentrations are converted to voxel masses;
source/sink rates and finite budgets are total rates/masses, with uniform masks.
The axis strategy avoids capacity × physical-support storage. Transient device
resource initialization buffers are released after import.
"""
function device_cellular_initialize(w::World,positions;capacity=4096,hidden=8,seed=w.seed,
        biomass=1.0,energy=1.0,metabolites=.12,initial=(.2,.2,0.0),
        potentials=(4.,4.,0.),diffusion=(.08,.08,.02),source_rates=(0.,0.,0.),
        sink_rates=(0.,0.,0.),source_remaining=(0.,0.,0.),sink_remaining=(0.,0.,0.),
        support_strategy=:axis,radius_bound=2f0,bin_width=radius_bound/2,
        mechanical_substeps=16,kwargs...)
    size(positions,1)==3 || throw(ArgumentError("positions must have three rows"))
    n=size(positions,2);h=Int(hidden);xyz=Float64.(Array(positions))
    0<=n<=capacity && capacity>0 || throw(ArgumentError("invalid founder capacity"))
    h>=2 && all(isfinite,xyz) || throw(ArgumentError("invalid founder state"))
    pars=cellular_parameters(;kwargs...)
    cellular_validate_founder_inventory(pars,biomass,energy,metabolites)
    all(v->length(v)==3,(initial,potentials,diffusion,source_rates,sink_rates,source_remaining,sink_remaining)) || throw(ArgumentError("three resource species required"))
    all(x->isfinite(x)&&x>=0,Iterators.flatten((initial,potentials,diffusion,source_rates,sink_rates,source_remaining,sink_remaining))) || throw(ArgumentError("invalid resource initialization"))
    nd=w.grid.size[3]==1 ? 2 : 3;volume=prod(w.grid.spacing[d] for d in 1:nd)
    fields=map(initial) do concentration
        a=KernelAbstractions.zeros(w.backend,Float32,w.grid.size...)
        mass=Float32(concentration*volume)
        isfinite(mass) || throw(ArgumentError("initial concentration exceeds Float32 mass range"))
        fill!(a,mass);a
    end
    masks=ntuple(3) do _
        a=KernelAbstractions.zeros(w.backend,Float32,w.grid.size...);fill!(a,1f0);a
    end
    for d in 1:3
        xyz[d,:].=d>nd ? 0. : mod.(xyz[d,:],w.grid.extent[d])
    end
    reference=(medium=nothing,resources=(medium=nothing,fields=collect(fields),potentials=potentials,
        diffusion=diffusion,source_rates=map(x->x/prod(w.grid.size),source_rates),
        sink_rates=map(x->x/prod(w.grid.size),sink_rates),source_masks=collect(masks),
        sink_masks=collect(masks),source_remaining=collect(source_remaining),sink_remaining=collect(sink_remaining),
        source_totals=ntuple(_->prod(w.grid.size),3),sink_totals=ntuple(_->prod(w.grid.size),3)),
        capacity=capacity,seed=seed,parameters=pars,stencil=(grid=w.grid,),
        ids=collect(1:n),parents=zeros(Int,n),lineages=collect(1:n),generations=zeros(Int,n),
        positions=xyz,biomass=fill(Float64(biomass),n),energy=fill(Float64(energy),n),
        metabolites=fill(Float64(metabolites),2,n),age=zeros(n),hidden=zeros(Float32,h,n),
        adaptation=zeros(Float32,h,n),last_actions=zeros(Float32,CELLULAR_ACTIONS,n),
        weights=cellular_founder_weights(seed,h,n),plastic_marks=zeros(Float32,h+CELLULAR_ACTIONS,2h+CELLULAR_SENSORS,n),
        bonds=(endpoints=Tuple{Int,Int}[],rest_lengths=Float64[]),next_id=n+1)
    device_cellular_import(w,reference;capacity,radius_bound,bin_width,mechanical_substeps,seed,support_strategy)
end
