#!/usr/bin/env julia
# Native fused material benchmark; optional preserved-kernel comparison.
using Interior
using Interior: JSON3
isdefined(Interior,:RecurrentMaterial) || Base.include(Interior,joinpath(@__DIR__,"..","src","recurrent_material.jl"))
using Interior: RecurrentMaterial
include(joinpath(@__DIR__,"recurrent_runtime.jl"))
Base.include(Interior,joinpath(@__DIR__,"..","test","recurrent_material_reference.jl"))
function runbatch(w,p,n,kernel=Interior.recurrent_material_kernel!)
    for _ in 1:n
        if kernel===Interior.recurrent_material_kernel!
            advance!(w,p,.05)
            continue
        end
        Interior.check_dt(w,p,.05)
        kernel(w.backend,256)(p.scratch,p.state,w.scratch.nutrient,w.scratch.deposit,
            w.fields.nutrient,w.fields.deposit,w.fields.activity,w.grid.spacing,.05f0,
            p.coupling,p.feedback,p.writing,p.erosion,p.contrast;ndrange=w.grid.size)
        p.state,p.scratch=p.scratch,p.state
        oldn=w.fields.nutrient;oldd=w.fields.deposit
        w.fields=merge(w.fields,(nutrient=w.scratch.nutrient,deposit=w.scratch.deposit))
        w.scratch=merge(w.scratch,(nutrient=oldn,deposit=oldd))
    end
    Interior.KernelAbstractions.synchronize(w.backend)
end
function benchmark_material(backend,side,steps,kernel)
    w=World(Grid(side,side,side;extent=Float32.((side,side,side)));backend,seed=31)
    p=RecurrentMaterial(w;coupling=.28,feedback=2,writing=.5,erosion=.035,contrast=1.1)
    runbatch(w,p,20,kernel)
    bytes=@allocated runbatch(w,p,20,kernel)
    times=Float64[]
    for start in 1:20:steps
        push!(times,@elapsed runbatch(w,p,min(20,steps-start+1),kernel))
    end
    elapsed=sum(times)
    (seconds=elapsed,voxel_steps_per_second=side^3*steps/elapsed,
     host_bytes_per_step=bytes/20,batch_seconds=times,
     finite=all(isfinite,Array(p.state)),material_mass=sum(Float64,Array(w.fields.nutrient))+sum(Float64,Array(w.fields.deposit)))
end
function benchmark_main()
    length(ARGS)>=4 || throw(ArgumentError("backend side steps output.json [--reference]"))
    kind=ARGS[1];side=parse(Int,ARGS[2]);steps=parse(Int,ARGS[3]);out=ARGS[4]
    side>0 && steps>0 || throw(ArgumentError("positive side/steps required"))
    ispath(out) && throw(ArgumentError("benchmark output exists; use a fresh destination"))
    backend=recurrent_backend(kind)
    # Heterogeneous material/resources, all eight channels and multiple steps.
    a=World(Grid(8,6,4;extent=(8,6,4));seed=31);b=World(a.grid;backend,seed=31)
    for I in CartesianIndices(a.fields.nutrient)
        a.fields.nutrient[I]=.1f0+.7f0*isodd(sum(Tuple(I)))
        a.fields.deposit[I]=.2f0+.4f0*iseven(sum(Tuple(I)))
    end
    copyto!(b.fields.nutrient,a.fields.nutrient);copyto!(b.fields.deposit,a.fields.deposit)
    x=RecurrentMaterial(a);y=RecurrentMaterial(b)
    runbatch(a,x,101,Interior.recurrent_material_reference_kernel!);runbatch(b,y,101)
    errors=(latent=maximum(abs.(x.state.-Array(y.state))),nutrient=maximum(abs.(a.fields.nutrient.-Array(b.fields.nutrient))),deposit=maximum(abs.(a.fields.deposit.-Array(b.fields.deposit))))
    maximum(values(errors))<2e-5 || error("CPU reference/GPU mismatch $errors")
    optimized=benchmark_material(backend,side,steps,Interior.recurrent_material_kernel!)
    reference="--reference" in ARGS ? benchmark_material(backend,side,steps,Interior.recurrent_material_reference_kernel!) : nothing
    result=(backend=kind,side=side,steps=steps,dt=.05,optimized=optimized,reference=reference,
        reference_over_optimized=reference===nothing ? nothing : reference.seconds/optimized.seconds,
        cpu_reference_errors=errors,resident_array_bytes=side^3*128,
        model_source_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"..","src","recurrent_material.jl")))))
    mkpath(dirname(out));temporary,io=mktemp(dirname(out))
    try
        JSON3.pretty(io,result);close(io);mv(temporary,out;force=false)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
    println(result)
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    isempty(ARGS) && throw(ArgumentError("backend side steps output.json [--reference]"))
    prepare_recurrent_backend!(ARGS[1])
    Base.invokelatest(benchmark_main)
end
