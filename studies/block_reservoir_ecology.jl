"""Experimental spatial lumping: fine-grid conservative transport and coarse
cell/catalyst/sorption reservoirs. No production process or checkpoint format.
Fine pools are the owned free/raw inventories; coarse pool arrays are temporary
representatives. Cell physiology, fair allocation, lifecycle and chemistry reuse
existing kernels. Sub-block competition, gradients, material localization and
reaction kinetics are deliberately approximated, not claimed equivalent.
"""
module BlockReservoirEcology
using Interior, Random, Statistics
include(joinpath(@__DIR__,"..","scripts","device_habitat_pilot_fixture.jl"))
include("nested_blocks.jl")
export block_fixture, block_step!, block_audit, block_storage_bytes, block_comparison,
    restrict_mass!, reconcile_mass!

# Construct a fresh carrier over the same exclusively borrowed arrays. This
# isolated adapter disables coarse diffusion; the fine transport stage owns it.
function without_diffusion(p)
    c=p.cells;q=merge(c.parameters,(diffusion=(0f0,0f0,0f0),))
    c2=Interior.DeviceCellular(c.state,c.environment,c.backend,q,c.geometry,c.attempted_steps,
        Interior.device_cellular_contract(c.state,c.environment,q,c.geometry))
    h=p.habitat;feed=merge(h.feedstock,(parameters=merge(h.feedstock.parameters,(diffusion=(0f0,0f0,0f0),)),))
    argument=Interior.device_habitat_argument(c2,h.arrays,h.parameters)
    ownership=(cells=c2,arrays=h.arrays,feedstock=feed,parameters=h.parameters,argument,
        array_shapes=map(size,h.arrays),feed_shapes=map(size,feed.arrays))
    h2=Interior.DeviceCatalyticHabitat(c2,h.arrays,feed,h.parameters,argument,ownership)
    Interior.DeviceConstructiveEcology(c2,h2)
end

struct BlockEcology{W,P,G,E,Q,B,F}
    world::W
    coarse::P
    fine_geometry::G
    free::E
    raw::E
    fine_parameters::Q
    blocks::B
    before_free::Matrix{Float32}
    before_raw::Matrix{Float32}
    body::F
end

block_indices(fine,coarse)=NestedBlocks(fine,coarse)

function restrict_mass!(coarse,fine,blocks)
    for species in axes(coarse,2),i in axes(coarse,1)
        amount=sum(j->fine[j,species],blocks[i])
        coarse[i,species]=Interior.dh_float_down(amount)
    end
    coarse
end

"""Apply the actual coarse mass change to fine donors, retaining unresolved
fine residuals. Withdrawals are proportional to available fine mass; additions
are uniform per volume. Float64 transfer workspace exposes rounding separately.
"""
function reconcile_mass!(fine,before,after,blocks)
    for species in axes(before,2),i in axes(before,1)
        ids=blocks[i];total=sum(j->fine[j,species],ids)
        delta=Float64(after[i,species])-Float64(before[i,species])
        target=total+delta
        target>=-64eps(Float64)*max(1.,total)||error("coarse withdrawal exceeded fine stock")
        if delta<0 && total>0
            scale=max(0.,target)/total
            for j in ids;fine[j,species]*=scale;end
        else
            addition=delta/length(ids)
            for j in ids;fine[j,species]+=addition;end
        end
    end
    fine
end

function block_fixture(;fine=32,coarse=16,founders=8,capacity=128,seed=911)
    fine>=coarse||error("fine lattice must resolve coarse")
    w,p=device_habitat_pilot_fixture(CPU();size=coarse,founders,capacity,seed,support_strategy=:materialized)
    q=p.cells.parameters;blocks=block_indices(fine,coarse);p=without_diffusion(p)
    makepool()=(fields=zeros(Float64,fine^3,3),scratch=zeros(Float64,fine^3,3),
        permeability=ones(Float32,fine^3),status=p.cells.environment.status)
    free=makepool();raw=makepool()
    for (source,dst) in ((p.cells.environment.fields,free.fields),(p.habitat.feedstock.arrays.fields,raw.fields))
        for s in 1:3,i in eachindex(blocks),j in blocks[i];dst[j,s]=Float64(source[i,s])/length(blocks[i]);end
    end
    for i in eachindex(blocks),j in blocks[i];free.permeability[j]=p.cells.environment.permeability[i];end
    g=merge(p.cells.geometry,(dims=(Int32(fine),Int32(fine),Int32(fine)),
        spacing=(4f0/fine,4f0/fine,4f0/fine),volume=(4f0/fine)^3))
    candidates=ntuple(d->Int32(min(fine,2ceil(Int,q.exchange_radius/g.spacing[d])+2)),3)
    g=merge(g,(candidates=candidates,support_strategy=Int32(1)))
    C=capacity;A=maximum(candidates);cs=p.cells.state
    body_state=(active=cs.active,position=cs.position,mass=cs.mass,metabolites=cs.metabolites,
        axis_indices=zeros(Int32,C,A,3),axis_distance2=zeros(Float32,C,A,3),
        support_normalization=zeros(Float32,C),support_fallback=zeros(Int32,C),support_count=zeros(Int32,C))
    body_environment=(occupancy=zeros(Float32,fine^3),permeability=free.permeability,status=p.cells.environment.status)
    b=BlockEcology(w,p,g,free,raw,q,blocks,zeros(Float32,coarse^3,3),zeros(Float32,coarse^3,3),
        (state=body_state,environment=body_environment))
    update_fine_porosity!(b)
    b
end

function update_fine_porosity!(b)
    s=b.body.state;e=b.body.environment;g=b.fine_geometry;q=b.fine_parameters
    fill!(e.occupancy,0)
    Interior.dc_support!(CPU(),128)(s,e,q,g;ndrange=length(s.active))
    Interior.dc_porosity_scatter!(CPU(),128)(s,e,g;ndrange=length(s.active))
    Interior.dc_porosity_commit!(CPU(),128)(e,q;ndrange=length(e.permeability))
    Interior.KernelAbstractions.synchronize(CPU())
    for i in eachindex(b.blocks),j in b.blocks[i]
        e.permeability[j]*=b.coarse.habitat.arrays.permeability[i]
    end
    nothing
end

function block_step!(b,dt)
    c=b.coarse.cells;e=c.environment
    dt=Float32(dt)
    isfinite(dt)&&dt>0||error("positive representable timestep required")
    limit=.45/(maximum(b.fine_parameters.diffusion)*sum(inv(Float64(h)^2) for h in b.fine_geometry.spacing))
    dt<=limit||error("fine diffusion CFL exceeded")
    for pool in (b.free,b.raw)
        Interior.dc_diffuse!(CPU(),128)(pool,b.fine_parameters,b.fine_geometry,Float32(dt);ndrange=size(pool.fields))
        Interior.KernelAbstractions.synchronize(CPU())
        copyto!(pool.fields,pool.scratch)
    end
    restrict_mass!(e.fields,b.free.fields,b.blocks)
    restrict_mass!(b.coarse.habitat.feedstock.arrays.fields,b.raw.fields,b.blocks)
    copyto!(b.before_free,e.fields);copyto!(b.before_raw,b.coarse.habitat.feedstock.arrays.fields)
    Interior.device_cellular_step!(b.world,c,dt;anabolism=b.coarse.habitat)
    reconcile_mass!(b.free.fields,b.before_free,e.fields,b.blocks)
    reconcile_mass!(b.raw.fields,b.before_raw,b.coarse.habitat.feedstock.arrays.fields,b.blocks)
    update_fine_porosity!(b)
    b
end

function block_audit(b)
    snapshot=Interior.device_habitat_snapshot(b.coarse);a=Interior.device_habitat_audit(snapshot)
    coarsefree=snapshot.cells.environment.fields;coarseraw=snapshot.feedstock.fields
    free_delta=sum(b.free.fields)-sum(Float64,coarsefree)
    raw_delta=sum(b.raw.fields)-sum(Float64,coarseraw)
    potential_delta=sum((sum(b.free.fields[:,k])-sum(Float64,coarsefree[:,k]))*b.fine_parameters.potentials[k] for k in 1:3)+
        sum((sum(b.raw.fields[:,k])-sum(Float64,coarseraw[:,k]))*b.coarse.habitat.feedstock.parameters.potentials[k] for k in 1:3)
    merge(a,(material=a.material+free_delta+raw_delta,potential=a.potential+potential_delta,
        minimum_fine=min(minimum(b.free.fields),minimum(b.raw.fields))))
end
function block_storage_bytes(b)
    support=sum(sizeof,getproperty(b.coarse.cells.state,k) for k in (:sites,:support,:support_count))
    (support=support,fine_pools=sum(sizeof,values(b.free))+sum(sizeof,values(b.raw)),
        block_descriptor=sizeof(b.blocks),transfer_workspace=sizeof(b.before_free)+sizeof(b.before_raw),
        fine_body_support=sum(sizeof,getproperty(b.body.state,k) for k in (:axis_indices,:axis_distance2,:support_normalization,:support_fallback,:support_count)),
        fine_body_occupancy=sizeof(b.body.environment.occupancy))
end

function field_error(x,y,volume)
    delta=Float64.(x).-Float64.(y)
    (relative_l2=sqrt(sum(abs2,delta)/max(sum(abs2,Float64.(y)),eps())),
        maximum_concentration_error=maximum(abs,delta)/volume,
        mass_l1_error=sum(abs,delta),
        per_species_relative_l2=[sqrt(sum(abs2,delta[:,k])/max(sum(abs2,Float64.(y[:,k])),eps())) for k in axes(x,2)])
end
function prolong_mass(coarse,blocks)
    fine=zeros(Float64,sum(length,blocks),size(coarse,2))
    for k in axes(coarse,2),i in eachindex(blocks),j in blocks[i]
        fine[j,k]=Float64(coarse[i,k])/length(blocks[i])
    end
    fine
end

function block_comparison(;fine=32,coarse=16,founders=8,capacity=128,steps=20,dt=.001,seed=911)
    # Tiny signature warmup, excluded from target timings.
    warm=block_fixture(;fine=8,coarse=8,founders=2,capacity=32,seed);block_step!(warm,dt);warm=nothing;GC.gc()
    b=block_fixture(;fine,coarse,founders,capacity,seed);initial=block_audit(b)
    blocked_seconds=@elapsed for _ in 1:steps;block_step!(b,dt);end
    final=block_audit(b)
    w,p=device_habitat_pilot_fixture(CPU();size=fine,founders,capacity,seed,support_strategy=:materialized)
    base_initial=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
    baseline_seconds=@elapsed for _ in 1:steps;Interior.device_cellular_step!(w,p.cells,dt;anabolism=p.habitat);end
    base_snapshot=Interior.device_habitat_snapshot(p);base_final=Interior.device_habitat_audit(base_snapshot)
    enzyme=reinterpret(Float32,b.coarse.habitat.arrays.enzyme_bits)
    (recipe=(fine,coarse,founders,capacity,steps,dt,seed),blocked_seconds,baseline_seconds,
        complete_step_speedup=baseline_seconds/blocked_seconds,
        storage=block_storage_bytes(b),baseline_support_bytes=sum(sizeof,getproperty(p.cells.state,k) for k in (:sites,:support,:support_count)),
        blocked_initial=initial,blocked_final=final,baseline_initial=base_initial,baseline_final=base_final,
        blocked_material_error=final.material-initial.material,
        blocked_potential_error=final.potential+final.work+final.heat-initial.potential,
        baseline_material_error=base_final.material-base_initial.material,
        baseline_potential_error=base_final.potential+base_final.work+base_final.heat-base_initial.potential,
        identities_equal=b.coarse.cells.state.ids==p.cells.state.ids&&b.coarse.cells.state.active==p.cells.state.active,
        free_field_error=field_error(b.free.fields,base_snapshot.cells.environment.fields,(4/fine)^3),
        raw_field_error=field_error(b.raw.fields,base_snapshot.feedstock.fields,(4/fine)^3),
        enzyme_field_error=field_error(prolong_mass(enzyme,b.blocks),reinterpret(Float32,base_snapshot.arrays.enzyme_bits),(4/fine)^3),
        bound_field_error=field_error(prolong_mass(b.coarse.habitat.arrays.bound,b.blocks),base_snapshot.arrays.bound,(4/fine)^3),
        permeability_maximum_error=maximum(abs.(b.free.permeability.-base_snapshot.cells.environment.permeability)),
        free_mass_difference=vec(sum(b.free.fields;dims=1).-sum(Float64.(base_snapshot.cells.environment.fields);dims=1)),
        raw_mass_difference=vec(sum(b.raw.fields;dims=1).-sum(Float64.(base_snapshot.feedstock.fields);dims=1)))
end
end
