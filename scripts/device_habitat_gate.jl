# Callable correctness gate: caller owns backend setup, exclusive GPU window,
# runtime logging and output. All per-step simulation stays in shared Julia code.
function device_habitat_3d_fixture(backend;size=16,support_strategy=:axis)
    grid=Grid(size,size,size;extent=(4,4,4));volume=prod(grid.spacing)
    host=World(grid;seed=911);fill!(host.fields.nutrient,0)
    food=ResourcePool(host;potentials=(4,4,0),initial=(.02volume,.02volume,0),
        diffusion=(.001,.001,.001))
    positions=Float64[1 3;2 2;1 3]
    cells=CellularEcology(host,positions,food;capacity=32,biomass=1.85,energy=1.,
        exchange_radius=.18,density=128.,mutation=.01,minimum_mass=.25)
    cells.weights[Base.size(cells.hidden,1)+10,end,:].=2f0
    raw=ResourcePool(host;names=(:raw_a,:raw_b),potentials=(6,6),
        initial=(4volume,4volume),diffusion=(.001,.001))
    habitat=Interior.CatalyticHabitat(host,cells,raw;secretion_rate=.2,regulation=false,investment_bias=1.,mutation=.01)
    matrix=Interior.SorptiveMatrix(host,habitat;capacity_per_enzyme=6.,adsorption_rate=4.,desorption_rate=.12,transport_barrier=20.)
    world=World(grid;backend,seed=911);fill!(world.fields.nutrient,0)
    device=Interior.DeviceCellular(world,cells;capacity=32,radius_bound=.3,mechanical_substeps=16,support_strategy)
    world,Interior.DeviceConstructiveEcology(world,device,habitat;matrix)
end
function device_habitat_3d_gate(backend;size=16,steps=20,support_strategy=:axis)
    1<=steps<=100 || throw(ArgumentError("gate requires1–100steps"))
    cw,cpu=device_habitat_3d_fixture(CPU();size,support_strategy)
    gw,gpu=device_habitat_3d_fixture(backend;size,support_strategy)
    initial=Interior.device_habitat_audit(Interior.device_habitat_snapshot(gpu))
    for _ in 1:steps
        advance!(cw,cpu,.001);advance!(gw,gpu,.001)
    end
    cs=Interior.device_habitat_snapshot(cpu);gs=Interior.device_habitat_snapshot(gpu)
    final=Interior.device_habitat_audit(gs);errors=Dict{String,Float64}()
    for (label,a,b) in (("cells",cs.cells.state,gs.cells.state),("environment",cs.cells.environment,gs.cells.environment),
            ("habitat",cs.arrays,gs.arrays),("feedstock",cs.feedstock,gs.feedstock))
        for key in keys(a)
            x=a[key];y=b[key]
            if eltype(x)<:AbstractFloat
                errors[label*"/"*string(key)]=isempty(x) ? 0. : maximum(abs,Float64.(x).-Float64.(y))
            elseif key===:enzyme_bits
                errors[label*"/enzyme"]=maximum(abs,Float64.(reinterpret(Float32,x)).-Float64.(reinterpret(Float32,y)))
            elseif label=="cells" && key in (:active,:ids,:parents,:lineages,:generations)
                x==y || error("CPU/device identity mismatch $key")
            end
        end
    end
    material_error=final.material-initial.material-final.supplied+initial.supplied+final.removed-initial.removed
    potential_error=final.potential-initial.potential+final.work-initial.work+final.heat-initial.heat-
        final.supplied_potential+initial.supplied_potential+final.removed_potential-initial.removed_potential
    maximum(values(errors))<1e-4 || error("3D habitat device parity failed")
    abs(material_error)<2e-4 && abs(potential_error)<5e-4 || error("3D habitat conservation failed")
    final.population>initial.population && all(>(0),final.secreted) && all(>(0),final.solubilized) && all(>(0),final.bound) ||
        error("3D gate did not exercise birth/construction/release/sorption")
    checkpoint=Interior.device_habitat_checkpoint(gw,gpu)
    rw=World(gw.grid;backend,seed=911);fill!(rw.fields.nutrient,0)
    restored=Interior.device_habitat_restore(rw,checkpoint)
    Interior.device_habitat_snapshot(restored)==gs || error("3D checkpoint restore is not exact")
    (;size,steps,support_strategy,initial,final,material_error,potential_error,array_errors=errors,
        checkpoint_bitwise=true,source=Interior.DEVICE_HABITAT_SOURCE_PIN,base_source=Interior.DEVICE_CELLULAR_SOURCE_PIN)
end
