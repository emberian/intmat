#!/usr/bin/env julia
# Backend must be preloaded before this include for the pinned AMD runtime.
include("device_cellular_setup.jl")
include(joinpath(@__DIR__,"..","studies","multirate_resolved_mechanics.jl"))
using .MultirateResolvedMechanics
using Statistics

function mr_arrays(x)
    (x.cells.state,x.cells.environment,x.arrays,x.feedstock)
end
function mr_difference(x,y;exact=false,atol=2e-5)
    worst=0.;discrete=true
    for (group,(a,b)) in enumerate(zip(mr_arrays(x),mr_arrays(y))),key in keys(a)
        # Atomic linked-list insertion order and cached radius bits are not physical state.
        if !exact && ((group==1 && key==:next)||(group==2 && key==:heads))
            continue
        end
        aa=a[key];bb=b[key]
        if !exact && group==2 && key==:metadata
            aa=aa[setdiff(eachindex(aa),[7])];bb=bb[setdiff(eachindex(bb),[7])]
        end
        size(aa)==size(bb)||error("shape mismatch $key")
        if eltype(aa)<:AbstractFloat
            err=isempty(aa) ? 0. : maximum(abs.(Float64.(aa).-Float64.(bb)))
            worst=max(worst,err)
            (exact ? aa==bb : err<=atol)||error("array mismatch $key: $err")
        else
            discrete&=aa==bb
            aa==bb||error("discrete array mismatch $key")
        end
    end
    (maximum_float_error=worst,discrete_equal=discrete,exact,excluded_nonphysical_caches=exact ? String[] : ["state.next","environment.heads","environment.metadata[7] radius bits"])
end
function mr_prepare!(r,mode)
    s=r.process.cells.state
    # Explicit bounded staging, outside all timings; never device scalar access.
    if mode==:lifecycle
        mass=Array(s.mass);energy=Array(s.energy);mass[1]=.1f0;mass[2]=2f0;energy[2]=2.
        copyto!(s.mass,mass);copyto!(s.energy,energy)
    elseif mode==:contact
        position=Array(s.position);bonds=Array(s.bonds);ids=Array(s.ids);bid=Array(s.bond_ids);rest=Array(s.rest)
        weights=Array(s.weights);hidden=Array(s.hidden);h=Int(r.process.cells.geometry.hidden)
        position[1,:].=Float32[1,1,1];position[2,:].=Float32[1.32,1,1]
        for (i,j) in ((1,2),(2,1))
            bonds[i,1]=j;bid[i,1]=ids[j];rest[i,1]=.3f0
            weights[i,:,h+8].=0;weights[i,end,h+8]=2
        end
        hidden[1,1]=.5f0;hidden[2,1]=-.5f0
        for (a,b) in ((s.position,position),(s.bonds,bonds),(s.bond_ids,bid),(s.rest,rest),(s.weights,weights),(s.hidden,hidden));copyto!(a,b);end
    end
    r
end
function mr_advance!(r,fast_steps;native=false,cadence=8)
    ratio=native ? 1 : r.carrier.ratio
    fast_steps%ratio==0||error("duration misaligned")
    p=r.process.cells
    for k in 1:fast_steps÷ratio
        if native
            Interior.device_cellular_step!(r.world,p,r.carrier.fast_dt;anabolism=r.process.habitat,check_status=false)
        else
            resolved_step!(r;check_status=false)
        end
        (k*ratio)%cadence==0&&Interior.device_cellular_status!(p)
    end
    Interior.device_cellular_status!(p)
    r
end
mr_snapshot(r)=Interior.device_habitat_snapshot(r.process)
function mr_residual(a,b)
    (material=a.material-b.material-a.supplied+b.supplied+a.removed-b.removed,
        potential=a.potential-b.potential+a.work-b.work+a.heat-b.heat-a.supplied_potential+b.supplied_potential+a.removed_potential-b.removed_potential)
end
function mr_gates(backend)
    rows=NamedTuple[]
    for mode in (:ordinary,:lifecycle,:contact)
        native=mr_prepare!(resolved_fixture(;size=8,founders=2,capacity=32,ratio=1,fast_controller=true,backend,energy_type=Float64),mode)
        split=mr_prepare!(resolved_fixture(;size=8,founders=2,capacity=32,ratio=1,fast_controller=true,backend,energy_type=Float64),mode)
        eltype(native.process.cells.state.energy)===Float64||error("energy precision was ignored")
        mr_advance!(native,8;native=true);mr_advance!(split,8)
        push!(rows,(mode,ratio=1,result=mr_difference(mr_snapshot(split),mr_snapshot(native);exact=true)))
        native=nothing;split=nothing;GC.gc()
    end
    for mode in (:ordinary,:lifecycle,:contact)
        cpu=mr_prepare!(resolved_fixture(;size=8,founders=2,capacity=32,ratio=4,fast_controller=true,backend=CPU(),energy_type=Float64),mode)
        gpu=mr_prepare!(resolved_fixture(;size=8,founders=2,capacity=32,ratio=4,fast_controller=true,backend,energy_type=Float64),mode)
        initial=Interior.device_habitat_audit(mr_snapshot(gpu))
        nfast=mode==:lifecycle ? 4 : 8
        mr_advance!(cpu,nfast);mr_advance!(gpu,nfast)
        snap=mr_snapshot(gpu);audit=Interior.device_habitat_audit(snap);residual=mr_residual(audit,initial)
        audit.finite&&abs(residual.material)<1e-4&&abs(residual.potential)<1e-4||error("gate closure failed")
        if mode==:lifecycle
            alive=findall(==(1),snap.cells.state.active)
            sort(snap.cells.state.ids[alive])==[2,3]||error("lifecycle gate missing death/division")
            child=only(filter(i->snap.cells.state.ids[i]==3,alive))
            snap.cells.state.parents[child]==2&&all(iszero,snap.arrays.last_secretion[child,:])||error("inheritance mismatch")
        elseif mode==:contact
            snap.cells.state.bonds[1,1]==2&&snap.cells.state.bonds[2,1]==1||error("lost contact")
            h=Int(gpu.process.cells.geometry.hidden)
            abs(snap.cells.state.inputs[1,h+1])>0||error("cached neighbour communication absent")
        end
        push!(rows,(mode,ratio=4,result=mr_difference(snap,mr_snapshot(cpu)),residual))
        cpu=nothing;gpu=nothing;snap=nothing;GC.gc()
    end
    println("MULTIRATE_GATES_COMPLETE");flush(stdout)
    rows
end
function multirate_gate_main()
    out=only(ARGS);ispath(out)&&error("fresh output required")
    size=parse(Int,get(ENV,"MULTIRATE_SIZE","64"));founders=parse(Int,get(ENV,"MULTIRATE_FOUNDERS","256"))
    capacity=parse(Int,get(ENV,"MULTIRATE_CAPACITY","4096"));steps=parse(Int,get(ENV,"MULTIRATE_FAST_STEPS","256"))
    pairs=parse(Int,get(ENV,"MULTIRATE_PAIRS","3"));dt=Float32(.001)
    size in (8,16,32,64,128)&&0<founders<=512&&capacity>=2founders&&0<steps<=512&&steps%8==0&&0<pairs<=3||error("bounded aligned fixture required")
    mkpath(dirname(abspath(out)));mkdir(out);started=time_ns()
    gates=mr_gates(backend)
    make(ratio)=resolved_fixture(;size,founders,capacity,ratio,fast_dt=dt,fast_controller=true,backend,energy_type=Float64)
    for ratio in (1,2,4)
        warm=make(ratio);mr_advance!(warm,8;native=ratio==1);warm=nothing;GC.gc()
    end
    rows=NamedTuple[];errors=NamedTuple[]
    for pair in 1:pairs
        snapshots=Dict{Int,Any}()
        for ratio in (isodd(pair) ? (1,2,4) : (4,2,1))
            r=make(ratio);initial=Interior.device_habitat_audit(mr_snapshot(r))
            bound=floor(Int,initial.material/r.process.cells.parameters.minimum_mass)+1
            capacity>=bound||error("capacity below funded headroom: $bound")
            GC.gc();seconds=@elapsed mr_advance!(r,steps;native=ratio==1)
            snap=mr_snapshot(r);audit=Interior.device_habitat_audit(snap)
            snap.cells.environment.metadata[5]==0||error("blocked birth in timing fixture")
            residual=mr_residual(audit,initial);audit.finite||error("nonfinite measured state")
            push!(rows,(pair,ratio,seconds,physical_duration=steps*Float64(dt),
                simulated_time_per_wall_second=steps*Float64(dt)/seconds,
                neural_updates=steps,mechanical_substep=Float64(dt)/16,
                biological_updates=steps÷ratio,funded_population_upper_bound=bound,
                energy_type=string(eltype(r.process.cells.state.energy)),audit,residual))
            snapshots[ratio]=snap;r=nothing;GC.gc()
            println("MULTIRATE_SAMPLE pair=",pair," ratio=",ratio," seconds=",seconds);flush(stdout)
        end
        ref=snapshots[1]
        for ratio in (2,4)
            x=snapshots[ratio];volume=(4/size)^3
            push!(errors,(pair,ratio,
                identities_equal=x.cells.state.ids==ref.cells.state.ids&&x.cells.state.active==ref.cells.state.active,
                free=MultirateResolvedMechanics.Shared.field_error(x.cells.environment.fields,ref.cells.environment.fields,volume),
                bound=MultirateResolvedMechanics.Shared.field_error(x.arrays.bound,ref.arrays.bound,volume),
                enzyme=MultirateResolvedMechanics.Shared.field_error(reinterpret(Float32,x.arrays.enzyme_bits),reinterpret(Float32,ref.arrays.enzyme_bits),volume)))
        end
        snapshots=nothing;GC.gc()
    end
    result=(source_commit=get(ENV,"DEVICE_SOURCE_COMMIT","unspecified"),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,
        habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,runtime=device_cellular_runtime(),
        driver_sha256=bytes2hex(Interior.SHA.sha256(read(@__FILE__))),
        cached_controller_source_sha256=MultirateResolvedMechanics.CACHED_CONTROLLER_SOURCE_SHA256,
        recipe=(size=size,founders=founders,capacity=capacity,fast_steps=steps,fast_dt=Float64(dt),pairs=pairs,energy_type="Float64"),
        timing_boundary="full coupled trajectory; identical batch8 fast-clock status checks; excludes construction and final snapshots; target warmups precede timings",
        checkpoint_support=false,gates,rows,errors,study_seconds=(time_ns()-started)/1e9)
    open(io->JSON3.pretty(io,result),joinpath(out,"gate.json"),"w")
    println("MULTIRATE_GATE_COMPLETE seconds=",result.study_seconds);flush(stdout)
end
if abspath(PROGRAM_FILE)==@__FILE__ || get(ENV,"MULTIRATE_INCLUDE_MAIN","0")=="1"
    multirate_gate_main()
end
