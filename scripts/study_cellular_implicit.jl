#!/usr/bin/env julia
# Bounded, opt-in mechanics study. Does not replace a process or resume a run.
using Interior,Serialization,JSON3,SHA,Statistics,LinearAlgebra
BLAS.set_num_threads(1)
solver_path=joinpath(@__DIR__,"..","studies","cellular_implicit_mechanics.jl")
solver_bytes=read(solver_path)
Base.include_string(Interior,String(copy(solver_bytes)),solver_path)
pair_path=joinpath(@__DIR__,"..","studies","cellular_pair_proximal.jl")
pair_bytes=read(pair_path)
Base.include_string(Interior,String(copy(pair_bytes)),pair_path)
coordinate_path=joinpath(@__DIR__,"..","studies","cellular_coordinate_proximal.jl")
coordinate_bytes=read(coordinate_path)
Base.include_string(Interior,String(copy(coordinate_bytes)),coordinate_path)
coordinate(args...;kwargs...)=Interior.cellular_implicit_mechanics!(args...;kwargs...,
    symmetric_motor="--midpoint" in ARGS,relaxer=Interior._CellularCoordinateCache(midpoint="--midpoint" in ARGS))
radial(args...;kwargs...)=Interior.cellular_implicit_mechanics!(args...;kwargs...,relaxer=Interior._CellularRadialCache())
reference_path=joinpath(@__DIR__,"..","studies","cellular_capped_reference.jl")
reference_bytes=read(reference_path)
Base.include_string(Interior,String(copy(reference_bytes)),reference_path)
reference_spec=Interior.install_cellular_capped_reference!()
fine_reference=reference_spec.method
function option(name,default=nothing)
    i=findfirst(==(name),ARGS);i===nothing ? default : ARGS[i+1]
end
fixture_path=abspath(option("--fixture"));out=abspath(option("--output"))
ispath(out)&&!isempty(readdir(out))&&error("study output must be fresh and empty")
mkpath(out)
write(joinpath(out,"worker.pid"),string(getpid())*"\n")
write(joinpath(out,"solver-source.jl"),solver_bytes)
write(joinpath(out,"pair-source.jl"),pair_bytes)
write(joinpath(out,"coordinate-source.jl"),coordinate_bytes)
write(joinpath(out,"reference-generator.jl"),reference_bytes)
write(joinpath(out,"fine-reference-source.jl"),reference_spec.source)
write(joinpath(out,"driver-source.jl"),read(@__FILE__))
fixture=open(deserialize,fixture_path)
grid=fixture.grid;ids=fixture.ids;radii=fixture.radii;q=fixture.parameters
motility=hasproperty(fixture,:motility) ? fixture.motility : Float64.(fixture.last_actions[1:3,:]).*q.speed
adhesion=hasproperty(fixture,:adhesion) ? fixture.adhesion : Float64.((fixture.last_actions[8,:].+1)./2)
kwargs=(;motility,adhesion,repulsion=q.repulsion,stiffness=q.stiffness,
    break_strain=q.break_strain,max_degree=q.max_degree,motor_cost=q.motor_cost)
duration=parse(Float64,option("--duration",string(hasproperty(fixture,:dt) ? fixture.dt : .1)))
quick="--quick" in ARGS;samples=quick ? 1 : 3
reference_splits=parse(Int,option("--reference-splits",quick ? "128" : "512"))
duration>0&&reference_splits>0||error("positive duration and reference-splits required")
function run(method,splits)
    x=copy(fixture.positions);b=deepcopy(fixture.bonds);energy=copy(fixture.energy)
    substeps=0;iterations=0;cg=0;evaluations=0;retries=0;ledger=0.;work=0.;loss=0.
    timing=@timed for _ in 1:1
        result=if method===explicit
            method(grid,ids,x,radii,b,duration;energy,kwargs...)
        elseif method===fine_reference
            method(grid,ids,x,radii,b,duration;energy,kwargs...,max_step=duration/splits)
        else
            method(grid,ids,x,radii,b,duration;energy,kwargs...,max_step=duration/splits,max_evaluations=256)
        end
        substeps+=result.substeps;ledger=max(ledger,abs(result.energy_residual))
        work+=sum(result.spent_motor_work);loss+=result.relaxation_energy_loss
        if hasproperty(result,:nonlinear_iterations)
            iterations+=result.nonlinear_iterations;cg+=result.cg_iterations;evaluations+=result.evaluations;retries+=result.retries
        end
    end
    (;x,b,energy,seconds=timing.time,bytes=timing.bytes,gc_seconds=timing.gctime,
      substeps,iterations,cg,evaluations,retries,ledger,work,loss)
end
explicit=Interior.cellular_mechanics!;implicit="--coordinate" in ARGS ? coordinate : "--pair" in ARGS ? radial : Interior.cellular_implicit_mechanics!
method_name="--coordinate" in ARGS ? ("--midpoint" in ARGS ? "coordinate-midpoint" : "coordinate") : "--pair" in ARGS ? "radial" : "implicit"
function position_difference(x,y)
    delta=x-y;nd=Interior._cellular_ndim(grid);maximum_distance=0.
    for i in axes(delta,2)
        distance2=0.
        for d in 1:nd
            L=Float64(grid.extent[d]);delta[d,i]=mod(delta[d,i]+L/2,L)-L/2;distance2+=delta[d,i]^2
        end
        maximum_distance=max(maximum_distance,sqrt(distance2))
    end
    (;rms=sqrt(sum(abs2,delta)/(nd*length(ids))),maximum=maximum_distance,
      center=ntuple(d->mean(view(delta,d,:)),nd))
end
# Compile measurement wrappers and each numerical branch before any samples.
println("Warming adaptive reference");flush(stdout)
run(fine_reference,reference_splits)
reference=run(fine_reference,reference_splits)
reference_uncertainty=nothing
if "--verify-reference" in ARGS
    println("Verifying capped reference at ",2reference_splits," substeps");flush(stdout)
    finer=run(fine_reference,2reference_splits)
    reference_uncertainty=position_difference(reference.x,finer.x)
    open(io->serialize(io,reference),joinpath(out,"reference-coarse.jls"),"w")
    reference=finer;reference_splits*=2
end
open(io->serialize(io,reference),joinpath(out,"reference.jls"),"w")
rows=Any[]
for (name,method,splits) in (("adaptive",explicit,1),(method_name,implicit,1),(method_name,implicit,2),(method_name,implicit,4),(method_name,implicit,8))
    println("Warming ",name," splits=",splits);flush(stdout)
    try
        run(method,splits)
    catch error
        failure=(;name,splits,status="failed",error=sprint(showerror,error))
        push!(rows,failure);println(JSON3.write(failure));flush(stdout)
        continue
    end
    println("Measuring ",name," splits=",splits);flush(stdout)
    measurements=Any[];last=nothing
    for _ in 1:samples
        GC.gc()
        result=try
            run(method,splits)
        catch error
            failure=(;name,splits,status="failed",error=sprint(showerror,error))
            push!(rows,failure);println(JSON3.write(failure));flush(stdout)
            break
        end
        last=result
        push!(measurements,(;result.seconds,result.bytes,result.gc_seconds,result.substeps,
            result.iterations,result.cg,result.evaluations,result.retries,result.ledger,result.work,result.loss))
    end
    last===nothing&&continue
    difference=position_difference(last.x,reference.x)
    a=Set(last.b.endpoints);b=Set(reference.b.endpoints);both=union(a,b)
    max_overlap=0.
    for (i,j) in cellular_pairs(grid,last.x,radii)
        distance=Interior._cellular_length(Interior._cellular_delta(grid,last.x,i,j))
        max_overlap=max(max_overlap,(radii[i]+radii[j]-distance)/(radii[i]+radii[j]))
    end
    row=(;name,splits,median_seconds=median(m.seconds for m in measurements),measurements,
        position_rms=difference.rms,position_max=difference.maximum,center_difference=difference.center,
        bond_jaccard=isempty(both) ? 1. : length(intersect(a,b))/length(both),
        max_overlap_fraction=max_overlap,energy_l1_difference=sum(abs,last.energy-reference.energy))
    push!(rows,row);println(JSON3.write(row));flush(stdout)
    open(io->serialize(io,last),joinpath(out,name*"-"*string(splits)*".jls"),"w")
end
metadata=(;scope="mechanics-only fixed bodies; not complete ecology steps",
    refinement="one fixed outer interval and funding horizon; cap mechanical h by duration/splits",duration,reference_splits,reference_uncertainty,samples,
    fixture=fixture_path,fixture_sha256=bytes2hex(sha256(read(fixture_path))),
    solver_sha256=bytes2hex(sha256(solver_bytes)),pair_sha256=bytes2hex(sha256(pair_bytes)),coordinate_sha256=bytes2hex(sha256(coordinate_bytes)),
    reference_sha256=bytes2hex(sha256(reference_spec.source)),reference_contacts_sha256=reference_spec.contacts_sha256,
    runtime=(;julia=string(VERSION),threads=Threads.nthreads(),blas_threads=BLAS.get_num_threads(),machine=Sys.MACHINE,kernel=Sys.KERNEL,finished_at=time()),
    source=Interior.arena_sources(),population=length(ids),grid=grid.size,extent=grid.extent,rows)
write(joinpath(out,"results.json"),JSON3.write(metadata))
