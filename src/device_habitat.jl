"""Device-resident paid catalyst and sorptive habitat for explicit new runs.

The cell carrier retains its declared Float32 law. Conservative local transfer
rounding uses Float64 scalar arithmetic, so this first helper supports CPU and
AMDGPU ROCBackend; Metal is deliberately rejected. No per-cell state crosses to
host during a step. Catalyst values use Float32 bit patterns in UInt32 storage
so bounded recipient additions can use ordinary integer compare-and-swap.
"""
mutable struct DeviceCatalyticHabitat{C,A,F,Q,K,O}
    cells::C
    arrays::A
    feedstock::F
    parameters::Q
    argument::K
    ownership::O
end
struct DeviceConstructiveEcology{C,H}
    cells::C
    habitat::H
end
const DeviceCatalyticArgument=NamedTuple{(:requests,:last_secretion,:enzyme_bits,:support_descriptor,
    :ledger,:ledger_low,:inert,:inert_low,:q,:status)}

function device_habitat_arrays(backend,C,V,h)
    z(T,d...)=KernelAbstractions.zeros(backend,T,d...)
    (weights=z(Float32,C,h+10,3),expression=z(Float32,C,3),requests=z(Float32,C,2),
     last_secretion=z(Float32,C,2),enzyme_bits=z(UInt32,V,2),bound=z(Float32,V,2),
     permeability=z(Float32,V),release=z(Float32,V,2),
     cell_ledger=z(Float32,C,10),cell_ledger_low=z(Float32,C,10),
     cell_inert=z(Float32,C),cell_inert_low=z(Float32,C),
     grid_ledger=z(Float32,V,12),grid_ledger_low=z(Float32,V,12),
     grid_inert=z(Float32,V),grid_inert_low=z(Float32,V),observation=z(Float32,16))
end
function device_habitat_feed_arrays(backend,V)
    z(T,d...)=KernelAbstractions.zeros(backend,T,d...)
    (fields=z(Float32,V,3),scratch=z(Float32,V,3),permeability=z(Float32,V),
     source_masks=z(Float32,V,3),sink_masks=z(Float32,V,3),
     resource_ledger=z(Float32,3,6),resource_ledger_low=z(Float32,3,6),quota=z(Float32,3,2))
end
function device_habitat_argument(c,a,q)
    (requests=a.requests,last_secretion=a.last_secretion,enzyme_bits=a.enzyme_bits,
     support_descriptor=device_support_descriptor(c.state,c.geometry),
     ledger=a.cell_ledger,ledger_low=a.cell_ledger_low,inert=a.cell_inert,inert_low=a.cell_inert_low,
     q=(secretion_work=q.secretion_work,),status=c.environment.status)
end

@inline dh_fail!(status,code)=KernelAbstractions.@atomic status[1]|=Int32(code)
@inline function dh_add!(hi,lo,i,k,value)
    first=Float32(value);second=Float32(Float64(value)-Float64(first))
    a,b=dc_pair_add(hi[i,k],lo[i,k],first);a,b=dc_pair_add(a,b,second)
    hi[i,k]=a;lo[i,k]=b
    nothing
end
@inline function dh_inert_add!(hi,lo,i,value)
    first=Float32(value);second=Float32(Float64(value)-Float64(first))
    a,b=dc_pair_add(hi[i],lo[i],first);a,b=dc_pair_add(a,b,second)
    hi[i]=a;lo[i]=b
    nothing
end
@inline function dh_float_down(value)
    result=Float32(value)
    Float64(result)>value && (result=prevfloat(result))
    max(0f0,result)
end
# Withdraw at least the exact delivered amount; any extra representable donor
# debit becomes explicitly owned inert mass (or heat for an energy debit).
@inline function dh_take_up(before::Float64,amount::Float64)
    after=max(0.0,before-amount)
    before-after<amount && after>0 && (after=prevfloat(after))
    after,before-after
end
@inline function dh_take_up(before::Float32,amount::Float64)
    after=Float32(max(0.,Float64(before)-amount))
    Float64(before)-Float64(after)<amount && after>0 && (after=prevfloat(after))
    after,Float64(before)-Float64(after)
end
# Withdraw no more than a request, used for kinetic source release.
@inline function dh_take_down(before::Float32,amount::Float64)
    after=Float32(max(0.,Float64(before)-amount))
    Float64(before)-Float64(after)>amount && (after=nextfloat(after))
    after,Float64(before)-Float64(after)
end
@inline function dh_add_down(before::Float32,amount::Float64)
    after=Float32(Float64(before)+amount)
    Float64(after)-Float64(before)>amount && (after=prevfloat(after))
    after,Float64(after)-Float64(before)
end

# Linearizable bounded deposit. Unlike floating atomic add, a successful update
# cannot credit more mass than this donor requested. Unrepresentable export is
# left in the donor. A bounded contention failure latches an incomplete step.
@inline function dh_catalyst_deposit!(bits,status,site,k,want::Float64)
    want<=0 && return 0.
    old=KernelAbstractions.@atomic :monotonic bits[site,k]
    for attempt in 1:4096
        before=reinterpret(Float32,old)
        if !(isfinite(before)&&before>=0)
            dh_fail!(status,8);return 0.
        end
        after,gained=dh_add_down(before,want)
        gained<=0 && return 0.
        if !isfinite(after)
            dh_fail!(status,8);return 0.
        end
        result=KernelAbstractions.@atomicreplace :monotonic bits[site,k] old=>reinterpret(UInt32,after)
        result.success && return gained
        old=result.old
    end
    dh_fail!(status,32);0.
end
@inline function dh_export!(arg,i,k,amount::Float64)
    actual=0.
    for (site,weight) in device_support_edges(arg.support_descriptor,i)
        want=min(max(0.,amount-actual),amount*Float64(weight))
        actual+=dh_catalyst_deposit!(arg.enzyme_bits,arg.status,site,k,want)
        arg.status[1]!=0 && break
    end
    actual
end

@inline device_anabolic_planned_work(planned,arg::DeviceCatalyticArgument,i)=
    planned+arg.q.secretion_work*(arg.requests[i,1]+arg.requests[i,2])

@inline function device_anabolic_build!(arg::DeviceCatalyticArgument,i,a,b,mass,energy,
        growth_request,growth_work,potential_a,potential_b)
    first=Float64(arg.requests[i,1]);second=Float64(arg.requests[i,2]);growth=Float64(growth_request)
    total=growth+first+second;workplan=Float64(growth_work)*growth+Float64(arg.q.secretion_work)*(first+second)
    available=2min(Float64(a),Float64(b))
    scale=min(1.,total>0 ? available/total : 1.,workplan>0 ? Float64(energy)/workplan : 1.)
    growth*=scale;first*=scale;second*=scale
    dh_add!(arg.ledger,arg.ledger_low,i,1,Float64(growth_request))
    dh_add!(arg.ledger,arg.ledger_low,i,2,Float64(arg.requests[i,1]))
    dh_add!(arg.ledger,arg.ledger_low,i,3,Float64(arg.requests[i,2]))
    dh_add!(arg.ledger,arg.ledger_low,i,4,growth)
    dh_add!(arg.ledger,arg.ledger_low,i,5,first)
    dh_add!(arg.ledger,arg.ledger_low,i,6,second)
    # Growth can receive only its representable biomass increment. Debit donor
    # endpoints conservatively, retaining the excess as zero-potential residue.
    nextmass,build=dh_add_down(mass,growth)
    nexta,taken_a=dh_take_up(a,build/2);nextb,taken_b=dh_take_up(b,build/2)
    nextenergy,paid=dh_take_up(energy,build*Float64(growth_work))
    inert=taken_a+taken_b-build;heat=taken_a*potential_a+taken_b*potential_b
    work=build*Float64(growth_work);heat+=max(0.,paid-work)
    a=nexta;b=nextb;mass=nextmass;energy=nextenergy
    for k in 1:2
        requested=k==1 ? first : second
        amount=min(requested,2min(Float64(a),Float64(b)))
        arg.q.secretion_work>0 && (amount=min(amount,Float64(energy)/Float64(arg.q.secretion_work)))
        actual=dh_export!(arg,i,k,max(0.,amount))
        a,used_a=dh_take_up(a,actual/2);b,used_b=dh_take_up(b,actual/2)
        cost=actual*Float64(arg.q.secretion_work);energy,spent=dh_take_up(energy,cost)
        inert+=used_a+used_b-actual;heat+=used_a*potential_a+used_b*potential_b+max(0.,spent-cost);work+=cost
        arg.last_secretion[i,k]=Float32(actual)
        dh_add!(arg.ledger,arg.ledger_low,i,6+k,actual)
    end
    if !(isfinite(inert)&&inert>=0&&isfinite(work)&&isfinite(heat)&&a>=0&&b>=0&&energy>=0)
        dh_fail!(arg.status,64)
    end
    dh_inert_add!(arg.inert,arg.inert_low,i,max(0.,inert))
    # Base ledgers receive Float32 scalars. Preserve their conversion corrections
    # once in this owner's compensated ledgers, not as a second work/heat charge.
    work32=Float32(work);heat32=Float32(heat)
    dh_add!(arg.ledger,arg.ledger_low,i,9,work-Float64(work32))
    dh_add!(arg.ledger,arg.ledger_low,i,10,heat-Float64(heat32))
    a,b,mass,energy,work32,heat32,Float32(build)
end

include("device_habitat_sorption.jl")
include("device_habitat_kernels.jl")

function DeviceConstructiveEcology(w::World,c::DeviceCellular,reference::CatalyticHabitat;matrix=nothing)
    (w.backend isa CPU || nameof(typeof(w.backend))===:ROCBackend) ||
        throw(ArgumentError("device habitat currently requires CPU or AMD ROCBackend; scalar Float64 rounding is not Metal compatible"))
    device_cellular_validate(w,c)
    C=Int(c.geometry.capacity);V=prod(Int.(c.geometry.dims));h=Int(c.geometry.hidden);n=length(reference.ids)
    n<=C && size(reference.weights)==(3,h+10,n) || throw(DimensionMismatch("habitat readout/capacity mismatch"))
    hostids=Array(c.state.ids);hostactive=Array(c.state.active)
    hostids[1:n]==reference.ids && count(==(1),hostactive)==n || throw(ArgumentError("habitat initial identities differ from device cells"))
    f=reference.substrates
    length(f.fields)==2 && all(a->size(a)==w.grid.size,f.fields) || throw(DimensionMismatch("two matching feedstock fields required"))
    matrix===nothing || matrix.habitat===reference || throw(ArgumentError("matrix belongs to another habitat"))
    mq=matrix===nothing ? (capacity_per_enzyme=0.,adsorption_rate=0.,desorption_rate=0.,cross_binding=0.,transport_barrier=0.) : matrix.parameters
    q=merge(map(x->x isa Bool ? x : Float32(x),reference.parameters),map(Float32,mq),
        (seed=UInt32(reference.seed),food_potentials=c.parameters.potentials))
    fq=(potentials=(Float32(f.potentials[1]),Float32(f.potentials[2]),0f0),
        diffusion=(Float32(f.diffusion[1]),Float32(f.diffusion[2]),0f0),
        source_rates=(Float32(f.source_rates[1]),Float32(f.source_rates[2]),0f0),
        sink_rates=(Float32(f.sink_rates[1]),Float32(f.sink_rates[2]),0f0),
        source_totals=(Float32(sum(Float64,f.source_masks[1])),Float32(sum(Float64,f.source_masks[2])),0f0),
        sink_totals=(Float32(sum(Float64,f.sink_masks[1])),Float32(sum(Float64,f.sink_masks[2])),0f0))
    all(k->fq.potentials[k]>=c.parameters.potentials[k],1:2) || throw(ArgumentError("feedstock conversion must be downhill in device precision"))
    all(x->isfinite(x)&&x>=0,Iterators.flatten(v isa Tuple ? v : (v,) for v in values(fq))) || throw(ArgumentError("invalid device feedstock parameters"))
    dh_validate_parameters(q,fq)
    a=device_habitat_arrays(w.backend,C,V,h);owned=device_habitat_feed_arrays(w.backend,V)
    weights=zeros(Float32,C,h+10,3);weights[1:n,:,:].=permutedims(reference.weights,(3,2,1));copyto!(a.weights,weights)
    enzyme=hcat((vec(map(dh_float_down,x)) for x in reference.enzymes)...)
    bound=matrix===nothing ? zeros(Float32,V,2) : hcat((vec(map(dh_float_down,x)) for x in matrix.bound)...)
    copyto!(a.enzyme_bits,copy(reinterpret(UInt32,enzyme)));copyto!(a.bound,bound)
    inert=zeros(Float64,V);heat=zeros(Float64,V)
    inert[1]+=reference.inert_residual+(matrix===nothing ? 0. : matrix.inert_residual)
    for k in 1:2
        inert.+=vec(reference.enzymes[k]).-Float64.(enzyme[:,k])
        if matrix!==nothing
            loss=vec(matrix.bound[k]).-Float64.(bound[:,k]);inert.+=loss;heat.+=loss.*q.food_potentials[k]
        end
    end
    copyto!(a.grid_inert,Float32.(inert));copyto!(a.grid_inert_low,Float32.(inert.-Float64.(Float32.(inert))))
    gh=zeros(Float32,V,12);gl=similar(gh);fill!(gl,0);gh[:,11].=Float32.(heat);gl[:,11].=Float32.(heat.-Float64.(gh[:,11]))
    copyto!(a.grid_ledger,gh);copyto!(a.grid_ledger_low,gl)
    copyto!(owned.fields,hcat(vec(f.fields[1]),vec(f.fields[2]),zeros(Float32,V)))
    copyto!(owned.source_masks,hcat(vec(f.source_masks[1]),vec(f.source_masks[2]),zeros(Float32,V)))
    copyto!(owned.sink_masks,hcat(vec(f.sink_masks[1]),vec(f.sink_masks[2]),zeros(Float32,V)))
    copyto!(owned.permeability,vec(f.permeability))
    ledger=zeros(Float32,3,6);low=similar(ledger);fill!(low,0)
    for (column,values) in ((1,f.source_remaining),(2,f.sink_remaining))
        ledger[1:2,column].=Float32.(values);low[1:2,column].=Float32.(values.-Float64.(ledger[1:2,column]))
    end
    copyto!(owned.resource_ledger,ledger);copyto!(owned.resource_ledger_low,low)
    feed=(arrays=owned,parameters=fq,argument=merge(owned,(status=c.environment.status,)))
    argument=device_habitat_argument(c,a,q)
    ownership=(cells=c,arrays=a,feedstock=feed,parameters=q,argument=argument,array_shapes=map(size,a),feed_shapes=map(size,owned))
    habitat=DeviceCatalyticHabitat(c,a,feed,q,argument,ownership)
    result=DeviceConstructiveEcology(c,habitat)
    validate(w,result)
    device_cellular_launch!(c,dh_initial_collapse!,a,c.environment,q;ndrange=V)
    device_cellular_launch!(c,dh_permeability!,a,c.environment,c.geometry,q;ndrange=V)
    device_cellular_status!(c);result
end
function validate(w::World,p::DeviceConstructiveEcology)
    c=p.cells;h=p.habitat;o=h.ownership
    h.cells===c && o.cells===c && h.arrays===o.arrays && h.feedstock===o.feedstock && h.parameters===o.parameters && h.argument===o.argument ||
        throw(ArgumentError("device habitat ownership/roles were replaced"))
    map(size,h.arrays)==o.array_shapes && map(size,h.feedstock.arrays)==o.feed_shapes || throw(ArgumentError("device habitat array shapes changed"))
    # The same scratch is leased sequentially by food and feedstock exchange.
    # Only feed diffusion scratch may overlap grants; catalysts, bound material,
    # masks, fields and all other owned habitat roles must stay disjoint.
    device_resource_workspace_nonalias(device_resource_workspace(c),h.arrays,h.feedstock.arrays;allowed_transients=(:scratch,))
    device_cellular_validate(w,c);device_cellular_status!(c);nothing
end
function stability_limit(w,p::DeviceConstructiveEcology)
    c=p.cells;fq=p.habitat.feedstock.parameters;D=maximum(fq.diffusion)
    S=sum(inv(Float64(c.geometry.spacing[d])^2) for d in 1:Int(c.geometry.nd))
    min(device_cellular_timestep_limit(c),D>0 ? .45/(D*S) : Inf)
end
function advance!(w::World,p::DeviceConstructiveEcology,dt::Real)
    check_dt(w,p,dt)
    device_cellular_step!(w,p.cells,dt;anabolism=p.habitat)
    p
end
process_children(p::DeviceConstructiveEcology)=(p.cells,p.habitat)
composition_references(p::DeviceCatalyticHabitat)=(p.cells,)
process_ports(::DeviceConstructiveEcology)=ProcessPorts(reads=(:resources,:feedstock,:enzymes,:bound_food),
    writes=(:cellular_state,:resources,:feedstock,:enzymes,:bound_food,:permeability),
    memory="resident inherited construction readout and paid catalytic/sorptive habitat",
    law="shared anabolic budget; downhill release; closed-form finite-capacity sorption")
device_anabolic_argument(h::DeviceCatalyticHabitat)=h.argument
function device_anabolic_begin!(h::DeviceCatalyticHabitat,w,p,dt)
    f=h.feedstock;g=p.geometry;V=prod(Int.(g.dims))
    device_resource_phase!(p,f.argument,f.parameters,g,dt)
    device_cellular_launch!(p,dc_diffuse!,f.argument,f.parameters,g,dt;ndrange=(V,3))
    device_cellular_launch!(p,dc_commit_resources!,f.argument;ndrange=(V,3))
    nothing
end
function device_anabolic_prepare!(h::DeviceCatalyticHabitat,w,p,dt)
    device_cellular_launch!(p,dh_prepare!,h.arrays,h.feedstock.arrays,p.state,device_support_descriptor(p),p.environment,h.parameters,p.parameters,p.geometry,dt;ndrange=Int(p.geometry.capacity))
    nothing
end
function device_anabolic_inherit!(h::DeviceCatalyticHabitat,w,p,step)
    device_cellular_launch!(p,dh_inherit!,h.arrays,p.state,p.environment,h.parameters,p.geometry,step;ndrange=Int(p.geometry.capacity))
    nothing
end
function device_anabolic_finish!(h::DeviceCatalyticHabitat,w,p,dt)
    V=prod(Int.(p.geometry.dims))
    device_cellular_launch!(p,dh_react_sorb!,h.arrays,h.feedstock.arrays,p.environment,h.parameters,h.feedstock.parameters,p.geometry,dt;ndrange=V)
    device_cellular_launch!(p,dh_permeability!,h.arrays,p.environment,p.geometry,h.parameters;ndrange=V)
    device_cellular_launch!(p,dh_validate_complete!,h.arrays,h.feedstock.arrays,p.environment,p.geometry;ndrange=max(V,Int(p.geometry.capacity)))
    nothing
end

"""Explicit full transfer for audits/checkpoints, never called by a simulation step."""
function device_habitat_snapshot(p::DeviceConstructiveEcology;check_status=true)
    check_status && device_cellular_status!(p.cells)
    h=p.habitat
    (cells=device_cellular_snapshot(p.cells;check_status=false),arrays=map(Array,h.arrays),
        feedstock=map(Array,h.feedstock.arrays),parameters=h.parameters,feed_parameters=h.feedstock.parameters)
end
function device_habitat_audit(x)
    base=device_cellular_audit(x.cells);a=x.arrays;feed=x.feedstock;q=x.parameters;fq=x.feed_parameters
    cellledger=Float64.(a.cell_ledger).+Float64.(a.cell_ledger_low)
    gridledger=Float64.(a.grid_ledger).+Float64.(a.grid_ledger_low)
    external=Float64.(feed.resource_ledger).+Float64.(feed.resource_ledger_low)
    enzyme=reinterpret(Float32,a.enzyme_bits)
    extra_inert=sum(Float64,a.cell_inert)+sum(Float64,a.cell_inert_low)+sum(Float64,a.grid_inert)+sum(Float64,a.grid_inert_low)
    material=base.material+sum(Float64,feed.fields)+sum(Float64,enzyme)+sum(Float64,a.bound)+extra_inert
    potential=base.potential+sum(sum(Float64,feed.fields[:,k])*fq.potentials[k] for k in 1:3)+
        sum(sum(Float64,a.bound[:,k])*q.food_potentials[k] for k in 1:2)
    (material=material,potential=potential,work=base.work+sum(cellledger[:,9]),
        heat=base.heat+sum(cellledger[:,10])+sum(gridledger[:,11]),
        supplied=base.supplied+sum(external[:,3]),removed=base.removed+sum(external[:,4]),
        supplied_potential=base.supplied_potential+sum(external[k,3]*fq.potentials[k] for k in 1:3),
        removed_potential=base.removed_potential+sum(external[k,4]*fq.potentials[k] for k in 1:3),
        population=base.population,completed_steps=base.completed_steps,
        secreted=vec(sum(cellledger[:,7:8];dims=1)),solubilized=vec(sum(gridledger[:,1:2];dims=1)),
        bound=vec(sum(Float64.(a.bound);dims=1)),enzyme=vec(sum(Float64.(enzyme);dims=1)),inert=extra_inert,
        finite=base.finite&&all(a->!(eltype(a)<:AbstractFloat)||all(isfinite,a),values(a))&&
            all(a->!(eltype(a)<:AbstractFloat)||all(isfinite,a),values(feed))&&all(isfinite,enzyme))
end

function dh_validate_parameters(q,fq)
    expected=(:secretion_rate,:secretion_work,:allocation_exponent,:catalytic_rate,:half_saturation,
        :decay_rate,:investment_bias,:allocation_bias_scale,:mutation,:mutation_probability,:regulation,
        :capacity_per_enzyme,:adsorption_rate,:desorption_rate,:cross_binding,:transport_barrier,:seed,:food_potentials)
    keys(q)==expected && keys(fq)==(:potentials,:diffusion,:source_rates,:sink_rates,:source_totals,:sink_totals) ||
        throw(ArgumentError("device habitat parameter roles differ"))
    all(k->k===:regulation ? q[k] isa Bool : k===:seed ? q[k] isa UInt32 :
        k===:food_potentials ? q[k] isa NTuple{3,Float32} : q[k] isa Float32,keys(q)) &&
        all(v->v isa NTuple{3,Float32},values(fq)) || throw(ArgumentError("device habitat parameter types differ"))
    all(isfinite,Iterators.flatten(v isa Tuple ? v : (v,) for v in values(q))) || throw(ArgumentError("nonfinite device habitat parameters"))
    nonnegative=(:secretion_rate,:secretion_work,:allocation_exponent,:catalytic_rate,:half_saturation,
        :decay_rate,:allocation_bias_scale,:mutation,:mutation_probability,:capacity_per_enzyme,
        :adsorption_rate,:desorption_rate,:cross_binding,:transport_barrier)
    all(k->getproperty(q,k)>=0,nonnegative) && q.half_saturation>0 && q.allocation_exponent>=1 &&
        q.mutation_probability<=1 && q.cross_binding<=1 || throw(ArgumentError("invalid device habitat parameters"))
    all(x->isfinite(x)&&x>=0,Iterators.flatten(v isa Tuple ? v : (v,) for v in values(fq))) &&
        all(k->fq.potentials[k]>=q.food_potentials[k],1:2) || throw(ArgumentError("invalid/downhill feedstock contract"))
    nothing
end
const DEVICE_HABITAT_SOURCE_PIN=(
    helper=bytes2hex(SHA.sha256(read(joinpath(@__DIR__,"device_habitat.jl")))),
    kernels=bytes2hex(SHA.sha256(read(joinpath(@__DIR__,"device_habitat_kernels.jl")))),
    sorption=bytes2hex(SHA.sha256(read(joinpath(@__DIR__,"device_habitat_sorption.jl")))))
function device_habitat_checkpoint(w::World,p::DeviceConstructiveEcology)
    validate(w,p);h=p.habitat;c=p.cells
    device_cellular_launch!(c,dh_validate_complete!,h.arrays,h.feedstock.arrays,c.environment,c.geometry;
        ndrange=max(Int(c.geometry.capacity),prod(Int.(c.geometry.dims))))
    device_cellular_status!(c)
    arrays=map(Array,h.arrays);feedstock=map(Array,h.feedstock.arrays)
    dh_validate_host_arrays(arrays,feedstock,h.parameters,h.feedstock.parameters)
    (format="Interior.DeviceConstructiveEcology",version=1,julia_version=VERSION,source=DEVICE_HABITAT_SOURCE_PIN,
        cells=device_cellular_checkpoint(w,c),arrays=arrays,feedstock=feedstock,
        parameters=h.parameters,feed_parameters=h.feedstock.parameters)
end
# Strict host payload audit belongs to checkpoint/restore, never the queued step.
function dh_validate_host_arrays(a,feed,q,fq)
    for group in (a,feed),array in values(group)
        array isa AbstractArray || throw(ArgumentError("habitat payload role is not an array"))
        eltype(array)<:AbstractFloat && !all(isfinite,array) && throw(ArgumentError("nonfinite habitat checkpoint array"))
    end
    for (hi,lo) in ((a.cell_ledger,a.cell_ledger_low),(a.grid_ledger,a.grid_ledger_low),
            (a.cell_inert,a.cell_inert_low),(a.grid_inert,a.grid_inert_low),
            (feed.resource_ledger,feed.resource_ledger_low))
        all(i->Float32(Float64(hi[i])+Float64(lo[i]))==hi[i],eachindex(hi)) ||
            throw(ArgumentError("nonnormalized compensated habitat pair"))
    end
    for key in (:expression,:requests,:last_secretion,:bound,:permeability,:release)
        all(>=(0),getproperty(a,key)) || throw(ArgumentError("negative habitat checkpoint inventory/readout"))
    end
    all(<=(1),a.expression) && all(<=(1),a.permeability) || throw(ArgumentError("invalid habitat checkpoint readout/permeability"))
    enzyme=reinterpret(Float32,a.enzyme_bits)
    all(x->isfinite(x)&&x>=0,enzyme) || throw(ArgumentError("invalid catalyst bits"))
    for k in 1:2,i in axes(enzyme,1)
        capacity=dh_float_down(Float64(q.capacity_per_enzyme)*((1-Float64(q.cross_binding))*enzyme[i,k]+Float64(q.cross_binding)*enzyme[i,3-k]))
        a.bound[i,k]<=capacity || throw(ArgumentError("bound food exceeds catalyst capacity"))
    end
    for (hi,lo,signed) in ((a.cell_ledger,a.cell_ledger_low,9),(a.grid_ledger,a.grid_ledger_low,13))
        for k in axes(hi,2),i in axes(hi,1)
            k<signed && Float64(hi[i,k])+Float64(lo[i,k])<0 && throw(ArgumentError("negative compensated habitat ledger"))
        end
    end
    for (hi,lo) in ((a.cell_inert,a.cell_inert_low),(a.grid_inert,a.grid_inert_low))
        all(i->Float64(hi[i])+Float64(lo[i])>=0,eachindex(hi)) || throw(ArgumentError("negative compensated retained material"))
    end
    for key in (:fields,:scratch,:source_masks,:sink_masks,:permeability,:quota)
        all(>=(0),getproperty(feed,key)) || throw(ArgumentError("negative feedstock checkpoint role"))
    end
    all(<=(1),feed.quota) && all(<=(1),feed.permeability) || throw(ArgumentError("invalid feedstock quota/permeability"))
    all(i->Float64(feed.resource_ledger[i])+Float64(feed.resource_ledger_low[i])>=0,eachindex(feed.resource_ledger)) ||
        throw(ArgumentError("negative compensated feedstock budget/counter"))
    for (mask,totals) in ((feed.source_masks,fq.source_totals),(feed.sink_masks,fq.sink_totals)),k in 1:3
        Float32(sum(Float64,view(mask,:,k)))==totals[k] || throw(ArgumentError("feedstock mask total differs from quota contract"))
    end
    for key in (:fields,:scratch,:source_masks,:sink_masks)
        all(iszero,view(getproperty(feed,key),:,3)) || throw(ArgumentError("unused third feedstock species is nonzero"))
    end
    all(iszero,view(feed.resource_ledger,3,:)) && all(iszero,view(feed.resource_ledger_low,3,:)) &&
        all(iszero,view(feed.quota,3,:)) && all(values(fq)) do value
            value isa Tuple ? iszero(value[3]) : true
        end || throw(ArgumentError("unused third feedstock species contract differs"))
    nothing
end
function device_habitat_restore(w::World,x)
    x.format=="Interior.DeviceConstructiveEcology" && x.version==1 && x.julia_version==VERSION ||
        throw(ArgumentError("unsupported device habitat checkpoint/runtime"))
    x.source==DEVICE_HABITAT_SOURCE_PIN || throw(ArgumentError("device habitat source pin differs"))
    (w.backend isa CPU || nameof(typeof(w.backend))===:ROCBackend) || throw(ArgumentError("device habitat requires CPU or AMD"))
    dh_validate_parameters(x.parameters,x.feed_parameters)
    c=device_cellular_restore(w,x.cells);g=c.geometry
    a=device_habitat_arrays(w.backend,Int(g.capacity),prod(Int.(g.dims)),Int(g.hidden))
    owned=device_habitat_feed_arrays(w.backend,prod(Int.(g.dims)))
    for (target,source) in ((a,x.arrays),(owned,x.feedstock))
        keys(target)==keys(source) || throw(ArgumentError("habitat checkpoint array roles differ"))
        for key in keys(target)
            dst=getproperty(target,key);src=getproperty(source,key)
            size(dst)==size(src)&&eltype(dst)==eltype(src) || throw(DimensionMismatch("habitat checkpoint array shape/type differs"))
            copyto!(dst,src)
        end
    end
    x.parameters.food_potentials==c.parameters.potentials || throw(ArgumentError("habitat food potentials differ from cells"))
    dh_validate_host_arrays(x.arrays,x.feedstock,x.parameters,x.feed_parameters)
    feed=(arrays=owned,parameters=x.feed_parameters,argument=merge(owned,(status=c.environment.status,)))
    argument=device_habitat_argument(c,a,x.parameters)
    ownership=(cells=c,arrays=a,feedstock=feed,parameters=x.parameters,argument=argument,
        array_shapes=map(size,a),feed_shapes=map(size,owned))
    h=DeviceCatalyticHabitat(c,a,feed,x.parameters,argument,ownership)
    p=DeviceConstructiveEcology(c,h);validate(w,p)
    device_cellular_launch!(c,dh_validate_complete!,a,owned,c.environment,g;ndrange=max(Int(g.capacity),prod(Int.(g.dims))))
    device_cellular_status!(c);p
end

# Generic observations are explicit synchronization/download boundaries. Arena's
# portable observation/checkpoint adapter is separate from the native protocol.
arena_process_supported(::DeviceConstructiveEcology)=false
arena_process_supported(::DeviceCatalyticHabitat)=false
function owned_material_mass(h::DeviceCatalyticHabitat)
    device_cellular_status!(h.cells)
    a=h.arrays;feed=h.feedstock.arrays
    # CPU/AMD support the Float64 accumulator required by strict inventory.
    # Only scalar reduction results cross the backend boundary.
    sum(Float64,feed.fields)+sum(x->Float64(reinterpret(Float32,x)),a.enzyme_bits)+sum(Float64,a.bound)+
        sum(Float64,a.cell_inert)+sum(Float64,a.cell_inert_low)+sum(Float64,a.grid_inert)+sum(Float64,a.grid_inert_low)
end
owned_material_mass(p::DeviceConstructiveEcology)=owned_material_mass(p.cells)+owned_material_mass(p.habitat)
function process_state_finite(p::DeviceConstructiveEcology)
    h=p.habitat
    process_state_finite(p.cells) &&
        all(a->!(eltype(a)<:AbstractFloat)||all(isfinite,a),values(h.arrays)) &&
        all(a->!(eltype(a)<:AbstractFloat)||all(isfinite,a),values(h.feedstock.arrays)) &&
        all(x->isfinite(reinterpret(Float32,x)),h.arrays.enzyme_bits)
end

"""Detached generic observation with named cell, habitat and feedstock roles.
This full snapshot is explicit; scalar diagnostics do not call this method.
Use device_habitat_checkpoint/device_habitat_restore for strict restart.
"""
function record_process_fields(p::DeviceConstructiveEcology)
    c=record_process_fields(p.cells);h=p.habitat
    a=map(Array,h.arrays);feed=map(Array,h.feedstock.arrays)
    (;schema="interior-device-constructive-observation",schema_version=1,
        catalyst_storage="habitat_enzyme_bits contains UInt32 bit patterns of Float32 catalyst inventory",
        habitat_parameter_names=Tuple(keys(h.parameters)),habitat_parameters=Tuple(values(h.parameters)),
        feedstock_parameter_names=Tuple(keys(h.feedstock.parameters)),feedstock_parameters=Tuple(values(h.feedstock.parameters)),
        habitat_source_names=Tuple(keys(DEVICE_HABITAT_SOURCE_PIN)),habitat_source=Tuple(values(DEVICE_HABITAT_SOURCE_PIN)),
        (Symbol("cells_",name)=>value for (name,value) in pairs(c))...,
        (Symbol("habitat_",name)=>array for (name,array) in pairs(a))...,
        (Symbol("feedstock_",name)=>array for (name,array) in pairs(feed))...)
end
