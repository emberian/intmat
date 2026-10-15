#!/usr/bin/env julia
# Bounded native trajectory, complete checkpoints and same-backend restart probe.
include("device_cellular_setup.jl")
include("device_cellular_population.jl")
using Serialization
side=parse(Int,get(ENV,"DEVICE_CELLS_SIDE","8"));steps=parse(Int,get(ENV,"DEVICE_STEPS","1000"))
batch=parse(Int,get(ENV,"DEVICE_BATCH","8"));capture=parse(Int,get(ENV,"DEVICE_CAPTURE","50"))
restart_steps=parse(Int,get(ENV,"DEVICE_RESTART_STEPS","20"));mechanical_substeps=parse(Int,get(ENV,"DEVICE_MECHANICAL_SUBSTEPS","16"))
min(side,steps,batch,capture,restart_steps,mechanical_substeps)>0 || error("positive run limits required")
out=length(ARGS)>0 ? ARGS[1] : "outputs/device-cellular-trajectory"
ispath(out) && error("fresh output required")
mkpath(dirname(abspath(out)));mkdir(out);mkdir(joinpath(out,"frames"));mkdir(joinpath(out,"source"))
for file in ("device_cellular.jl","device_cellular_support.jl","device_cellular_initialization.jl","cellular_initialization.jl","device_cellular_kernels.jl","device_cellular_checkpoint.jl","device_cellular_protocol.jl")
    cp(joinpath(@__DIR__,"..","src",file),joinpath(out,"source",file))
end
for file in ("run_device_cellular.jl","device_cellular_setup.jl","device_cellular_population.jl")
    cp(joinpath(@__DIR__,file),joinpath(out,"source",file))
end
function capture_frame!(out,p,step)
    # Caller has just checked status. Transfer selected fields, then store only
    # active rows. Original slot indices preserve the meaning of bond slots.
    active=Array(p.state.active);slots=findall(==(1),active)
    names=(:ids,:parents,:lineages,:generations,:position,:radii,:mass,:energy,:metabolites,:hidden,:actions,:bonds,:bond_ids,:rest)
    source=NamedTuple{names}(map(key->Array(getproperty(p.state,key)),names))
    selected=map(a->ndims(a)==1 ? a[slots] : a[slots,:],source)
    fields=Array(p.environment.fields);permeability=Array(p.environment.permeability)
    frame=merge((step=step,slot_indices=Int32.(slots),),selected,
        (resource_fields=fields,permeability=permeability,geometry=p.geometry))
    serialize(joinpath(out,"frames","frame-$(lpad(step,7,'0')).jls"),frame)
    downloaded=sizeof(active)+sum(sizeof,values(source))+sizeof(fields)+sizeof(permeability)
    stored=sum(sizeof(v) for v in values(frame) if v isa AbstractArray)
    downloaded,stored
end
function run_trajectory!(w,p,out,steps,batch,capture)
    observations=Any[];capture_bytes=0;capture_download_bytes=0;delta=.01f0
    elapsed=@elapsed for step in 1:steps
        Interior.device_cellular_step!(w,p,delta;check_status=false);w.time+=Float64(delta)
        if step%batch==0 || step%capture==0 || step==steps
            push!(observations,(step=step,stats=Interior.device_cellular_statistics(p)))
        end
        if step%capture==0 || step==steps
            downloaded,stored=capture_frame!(out,p,step)
            capture_download_bytes+=downloaded;capture_bytes+=stored
        end
    end
    observations,capture_bytes,capture_download_bytes,elapsed
end
w,p=device_cellular_population(side,backend;mechanical_substeps)
initial=Interior.device_cellular_audit(Interior.device_cellular_snapshot(p))
remaining_sources=sum(Float64.(Array(p.environment.resource_ledger)[:,1]).+Float64.(Array(p.environment.resource_ledger_low)[:,1]))
population_material_bound=(initial.material+remaining_sources)/p.parameters.minimum_mass
headroom_forecast=(initial_population=initial.population,capacity=Int(p.geometry.capacity),
    initial_material=initial.material,remaining_external_material=remaining_sources,
    minimum_viable_mass=Float64(p.parameters.minimum_mass),conservative_live_count_bound=population_material_bound,
    capacity_covers_material_bound=Int(p.geometry.capacity)>1.01population_material_bound)
get(ENV,"DEVICE_REQUIRE_HEADROOM","0")=="1" && !headroom_forecast.capacity_covers_material_bound &&
    error("capacity does not cover conservative material forecast: $headroom_forecast")
serialize(joinpath(out,"initial.jls"),Interior.device_cellular_checkpoint(w,p))
# Compile on the real state, then restore exact initial state before timed trajectory.
Interior.device_cellular_step!(w,p,.01);Interior.device_cellular_statistics(p)
# Verify the queued-failure boundary on this backend before the sustained run.
copyto!(p.environment.status,Int32[8])
latched_before=Interior.device_cellular_snapshot(p;check_status=false)
for _ in 1:3;Interior.device_cellular_step!(w,p,.01;check_status=false);end
latched_after=Interior.device_cellular_snapshot(p;check_status=false)
latch_arrays_frozen=latched_before.state==latched_after.state && latched_before.environment==latched_after.environment
latch_arrays_frozen || error("queued failure boundary mutated device arrays")
p=nothing;GC.gc()
p=Interior.device_cellular_restore(w,deserialize(joinpath(out,"initial.jls")))
observations,capture_bytes,capture_download_bytes,elapsed=run_trajectory!(w,p,out,steps,batch,capture)
final=Interior.device_cellular_audit(Interior.device_cellular_snapshot(p))
serialize(joinpath(out,"final.jls"),Interior.device_cellular_checkpoint(w,p))
# Continue once, restore on the same backend, repeat; compare all persistent arrays.
for _ in 1:restart_steps;Interior.device_cellular_step!(w,p,.01;check_status=false);w.time+=Float64(.01f0);end
expected=Interior.device_cellular_snapshot(p)
p=nothing;GC.gc()
saved=deserialize(joinpath(out,"final.jls"))
p=Interior.device_cellular_restore(w,saved)
restore_exact=Interior.device_cellular_snapshot(p)==saved.snapshot
for _ in 1:restart_steps;Interior.device_cellular_step!(w,p,.01;check_status=false);w.time+=Float64(.01f0);end
actual=Interior.device_cellular_snapshot(p)
function restart_comparison(expected,actual)
    errors=Dict{String,Float64}();integer_mismatches=String[]
    for group in (:state,:environment),key in keys(getproperty(expected,group))
        a=getproperty(getproperty(expected,group),key);b=getproperty(getproperty(actual,group),key)
        if eltype(a)<:AbstractFloat
            errors["$(group).$(key)"]=isempty(a) ? 0. : maximum(abs.(Float64.(a).-Float64.(b)))
        else
            a==b || push!(integer_mismatches,"$(group).$(key)")
        end
    end
    (exact_integer_arrays=isempty(integer_mismatches),integer_mismatches=integer_mismatches,maximum_float_error=maximum(values(errors)),array_errors=errors,
        attempted_steps_equal=expected.attempted_steps==actual.attempted_steps)
end
restart=restart_comparison(expected,actual)
# Explicit same-backend engineering replay tolerances for this short probe.
# Bin linked-list order is an implementation cache, not biological identity.
allowed_integer_differences=Set(["state.next","environment.heads"])
identity_ok=all(k->k in allowed_integer_differences,restart.integer_mismatches)
physical_keys=("state.position","state.mass","state.metabolites","state.energy","state.hidden","state.adaptation","state.weights","state.marks","environment.fields")
physical_error=maximum(restart.array_errors[k] for k in physical_keys)
restart_pass=restore_exact && identity_ok && restart.attempted_steps_equal && physical_error<=1e-4 && restart.maximum_float_error<=1e-3

result=(runtime=device_cellular_runtime(),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,backend=kind,amd_optlevel=get(ENV,"DEVICE_AMD_OPTLEVEL","2"),
    headroom_forecast=headroom_forecast,latch_arrays_frozen=latch_arrays_frozen,capacity=Int(p.geometry.capacity),side=side,steps=steps,dt=Float64(.01f0),batch=batch,capture=capture,capture_bytes=capture_bytes,capture_download_bytes=capture_download_bytes,
    seconds=elapsed,full_steps_per_second=steps/elapsed,initial=initial,final=final,restart_steps=restart_steps,restore_all_arrays_exact=restore_exact,restart=restart,restart_pass=restart_pass,
    restart_tolerances=(physical_arrays=1e-4,all_float_arrays=1e-3,allowed_integer_differences=collect(allowed_integer_differences)),
    material_error=final.material-initial.material-final.supplied+initial.supplied+final.removed-initial.removed,
    potential_error=final.potential-initial.potential-final.supplied_potential+initial.supplied_potential+final.removed_potential-initial.removed_potential+final.work-initial.work+final.heat-initial.heat,
    timing_boundary="full trajectory including compact observations, selected-field downloads and frame serialization; excludes compilation and complete checkpoint/audit I/O")
open(joinpath(out,"trajectory.json"),"w") do io;JSON3.pretty(io,result);end
open(joinpath(out,"observations.json"),"w") do io;JSON3.pretty(io,observations);end
println(result)
restart_pass || error("restart continuation exceeded declared tolerances; inspect trajectory.json")
