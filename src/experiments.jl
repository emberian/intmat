"""Descriptive process interfaces, not proofs of stability or conservation.

`writes` lists all modified world fields; `publishes` marks fields replaced by an
owned state or exclusive driver. `owned`, `drives`, `memory`, and `law` describe
process semantics. Custom processes should implement `process_ports`.
"""
Base.@kwdef struct ProcessPorts
    reads::Tuple=()
    writes::Tuple=()
    publishes::Tuple=()
    owned::Tuple=()
    drives::String="unspecified"
    memory::String="unspecified"
    law::String="unspecified; inspect implementation"
end
process_ports(p)=ProcessPorts()
process_ports(::Excitable)=ProcessPorts(reads=(:activity,:recovery,:deposit),writes=(:activity,:recovery),publishes=(:activity,:recovery),memory="world activity and recovery",law="explicit FitzHugh–Nagumo reaction and diffusion")
process_ports(::NeuralField)=ProcessPorts(reads=(:nutrient,:deposit,:pressure),writes=(:activity,),publishes=(:activity,),owned=(:latent1,:latent2),memory="two recurrent latent channels",law="fixed recurrent tanh target with explicit relaxation and diffusion")
process_ports(::Deposition)=ProcessPorts(reads=(:activity,:nutrient,:deposit),writes=(:nutrient,:deposit),memory="world material inventory",law="bounded nutrient/deposit exchange")
process_ports(p::Diffusion)=ProcessPorts(reads=(p.field,:deposit),writes=(p.field,),law="periodic symmetric flux diffusion")
process_ports(::WaveForcing)=ProcessPorts(reads=(:deposit,),writes=(:pressure,:vx,:vy,:vz,:activity),publishes=(:pressure,:vx,:vy,:vz),drives="standing sinusoidal pressure and velocity",law="prescribed wave plus additive activity forcing")
process_ports(::ThermalWriter)=ProcessPorts(reads=(:nutrient,:deposit),writes=(:nutrient,:deposit),owned=(:temperature,:dose),drives="parametric Gaussian heat source and shutter",memory="temperature and irreversible accumulated exposure",law="explicit heat diffusion/cooling and bounded material conversion")
process_ports(::AcousticField)=ProcessPorts(reads=(:deposit,:activity),writes=(:pressure,:activity),publishes=(:pressure,),owned=(:pressure,:previous,:envelope),drives="localized acoustic source",memory="two wave time levels and intensity envelope",law="damped heterogeneous scalar wave; additive activity forcing")
process_ports(::DielectricGrowth)=ProcessPorts(reads=(:nutrient,:deposit),writes=(:nutrient,:deposit,:activity),publishes=(:activity,),owned=(:potential,:occupied,:frontier),drives="fixed electrical boundary values",memory="irreversible occupied cells and event clock",law="quasistatic potential solve and stochastic frontier growth")
process_ports(::ParticleCarrier)=ProcessPorts(reads=(:nutrient,:deposit,:vx,:vy,:vz),writes=(:nutrient,:deposit),owned=(:positions,:hidden,:storage,:ids),memory="particle latent channels, inventory and stable IDs",law="particle motion and conservative material exchange")


process_ports(::PhaseSeparation)=ProcessPorts(reads=(:deposit,:nutrient,:activity),writes=(:nutrient,:activity),owned=(:phi,:chemical_potential,:interface),memory="conserved signed composition",law="explicit periodic Cahn–Hilliard flux with mineral pinning, permeability and activity relaxation")
process_ports(::BondNetwork)=ProcessPorts(owned=(:position,:velocity,:bonds,:damage,:peak_strain),drives="prescribed grip displacement and node load",memory="irreversible bond damage and mechanical energy ledgers",law="central-force springs, velocity Verlet, exponential damping and tensile failure")
process_ports(::Ecology)=ProcessPorts(reads=(:nutrient,:deposit,:vx,:vy,:vz),writes=(:nutrient,:deposit),owned=(:positions,:hidden,:storage,:energy,:weights,:ids,:lineages),drives="nutrient uptake funds abstract energy",memory="heritable recurrent controllers, generations and inventory ledgers",law="cell-local conservative exchange, bounded births and deaths with independent energy accounting")

"""Ordered, named process composition with a fixed Float32 step size.

Construction preflights every process. Two publishers of the same field require
`allow_overwrite=true`: tuple order then explicitly determines the last writer.
Additive updates and material exchanges may share fields without this opt-in.
Traits describe interfaces; they cannot establish that a composition is sound.
"""
struct Experiment{W,P,N}
    world::W
    processes::P
    names::N
    dt::Float32
end
function Experiment(w::World, processes::NamedTuple;dt,allow_overwrite=false)
    _experiment(w,Tuple(values(processes)),Tuple(keys(processes)),dt,allow_overwrite)
end
function Experiment(w::World, processes::Tuple;dt,allow_overwrite=false)
    if all(p->p isa Pair,processes)
        names=Tuple(Symbol(first(p)) for p in processes)
        ps=Tuple(last(p) for p in processes)
    else
        ps=processes
        names=ntuple(i->Symbol("process_",i),length(ps))
    end
    _experiment(w,ps,names,dt,allow_overwrite)
end
function _experiment(w,ps,names,dt,allow_overwrite)
    length(unique(names))==length(names) || throw(ArgumentError("process names must be unique"))
    isfinite(dt) && dt>0 && isfinite(Float32(dt)) && Float32(dt)>0 || throw(ArgumentError("dt must be positive and finite in Float32"))
    publishers=Dict{Symbol,Symbol}()
    for (name,p) in zip(names,ps)
        check_dt(w,p,dt)
        for field in process_ports(p).publishes
            if haskey(publishers,field) && !allow_overwrite
                throw(ArgumentError("$(publishers[field]) and $name both publish $field; specify allow_overwrite=true to accept this order"))
            end
            publishers[field]=name
        end
    end
    Experiment(w,ps,names,Float32(dt))
end

"""Explicit mapping from simulation state to render fields; never advances time."""
struct OpticalMap{F}
    map::F
end
(m::OpticalMap)(w::World,processes)=m.map(w,processes)
(m::OpticalMap)(e::Experiment)=m(e.world,e.processes)

struct ExperimentReport{M,S}
    metadata::M
    samples::S
end
"""Material inventory held outside the world; custom carriers should specialize this hook."""
owned_material_mass(p)=0.0
owned_material_mass(p::ParticleCarrier)=sum(Float64,p.storage;init=0.0)
owned_material_mass(p::Ecology)=sum(Float64,p.storage;init=0.0)+p.waste

"""Small host diagnostics; callback cadence does not change the integration schedule.

`material_mass` includes world inventory and declared process inventories. Finiteness
checks world fields and direct floating array fields of every process. Custom nested
state requires a custom diagnostics callback.
"""
function experiment_diagnostics(w::World,processes)
    fields=map(a->Array(a),values(w.fields))
    grid_mass=sum(Float64,Array(w.fields.nutrient))+sum(Float64,Array(w.fields.deposit))
    owned_mass=sum(owned_material_mass(p) for p in processes; init=0.0)
    owned_arrays=(getfield(p,n) for p in processes for n in fieldnames(typeof(p))
                  if getfield(p,n) isa AbstractArray{<:AbstractFloat})
    (material_mass=grid_mass+owned_mass,grid_mass=grid_mass,owned_mass=owned_mass,
     finite=all(a->all(isfinite,a),fields) && all(a->all(isfinite,Array(a)),owned_arrays),
     activity_max=maximum(Array(w.fields.activity)),
     nutrient_min=minimum(Array(w.fields.nutrient)),
     deposit_mass=sum(Float64,Array(w.fields.deposit)))
end

"""
    run!(experiment; steps, capture_every=1, observe=nothing,
         diagnostics=experiment_diagnostics)

Advance exactly `steps` fixed steps. Capture initial state, cadence boundaries,
and the final state (once). Callbacks receive `(world, processes, step)` for
`observe` and `(world, processes)` for diagnostics; they must be read-only.
Their returned values are deep-copied into samples so later updates cannot alter
past observations. `observe` can invoke `record_state` at these boundaries.
No adaptive stepping, automatic process switching, or optical blending occurs.
A repeated run reports step indices relative to that run and keeps world time.
"""
function run!(e::Experiment;steps::Integer,capture_every::Integer=1,
              observe=nothing,diagnostics=experiment_diagnostics)
    steps>=0 || throw(ArgumentError("steps must be nonnegative"))
    capture_every>0 || throw(ArgumentError("capture_every must be positive"))
    # Revalidate before even the initial observation (process state may have changed).
    for p in e.processes; check_dt(e.world,p,e.dt); end
    start=e.world.time
    samples=NamedTuple[]
    function capture(k)
        d=diagnostics===nothing ? nothing : diagnostics(e.world,e.processes)
        o=observe===nothing ? nothing : observe(e.world,e.processes,k)
        push!(samples,(step=k,time=e.world.time,diagnostics=deepcopy(d),observation=deepcopy(o)))
    end
    capture(0)
    for k in 1:steps
        step!(e.world,e.processes,e.dt)
        (k%capture_every==0 || k==steps) && capture(k)
    end
    metadata=(seed=e.world.seed,dt=e.dt,steps=Int(steps),start_time=start,
              end_time=e.world.time,process_order=e.names,
              process_types=map(p->string(typeof(p)),e.processes),
              ports=map(process_ports,e.processes))
    ExperimentReport(metadata,samples)
end

# Arena processes share the same composition and inventory interfaces.
process_ports(::ImplicitPhaseSeparation)=ProcessPorts(reads=(:deposit,:nutrient,:activity),writes=(:nutrient,:activity),owned=(:phi,:chemical_potential,:interface,:chemical_forcing),memory="conserved composition and external chemical-potential bias",law="constant-mobility semi-implicit periodic Cahn–Hilliard with the finite-difference Laplacian spectrum")
process_ports(::ActiveMedium)=ProcessPorts(reads=(:nutrient,:deposit),writes=(:nutrient,),owned=(:permeability,:signal,:signal_source,:remodeling,:forcing),memory="diffusing signals and decaying remodeling bias",law="conservative face transport through active permeability; drives phase chemical potential")
process_ports(::NeuralEcology)=ProcessPorts(reads=(:nutrient,:deposit),writes=(:nutrient,:deposit),owned=(:positions,:hidden,:storage,:energy,:weights,:plastic_marks,:ids,:lineages),drives="chemical uptake funds abstract energy",memory="inherited recurrent genotype and separate local plastic memory",law="funded neural actions and material exchange in a sensed compartment environment")
owned_material_mass(p::NeuralEcology)=sum(Float64,p.storage;init=0.0)+p.waste
