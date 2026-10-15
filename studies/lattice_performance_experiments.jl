# Diagnostic variants only. Never modifies the frozen ecological law or old runs.
using Interior,SHA,Profile,Serialization,Statistics
using Interior:JSON3
include("lattice_ecology.jl")
function same_state(a,b)
    a.tick==b.tick||return false
    all(name->name==:weights||getfield(a,name)==getfield(b,name),fieldnames(typeof(a)))||return false
    a.weights==permutedims(b.weights,(3,1,2))
end
function variant_state(w,p,mod)
    q=mod.LatticeEcology((name==:weights ? permutedims(p.weights,(2,3,1)) : deepcopy(getfield(p,name)) for name in fieldnames(typeof(p)))...)
    v=deepcopy(w);Interior.validate(v,q);@assert same_state(p,q);v,q
end
function timed_window(w,p,n)
    elapsed=@elapsed for _ in 1:n;step!(w,p,.1);end
    elapsed
end
function performance_main(output;bounds=false)
    output=abspath(output);Threads.nthreads()==1||error("CPU1 only");ispath(output)&&error("fresh output directory required");mkpath(output)
    source=read(joinpath(@__DIR__,"lattice_ecology.jl"),String)
    variant=replace(source,"module LatticeEcologyStudy"=>"module LatticeEcologyPackedStudy",
        "weights=lattice_genotypes(;loci,seed=genotype_seed)"=>"weights=permutedims(lattice_genotypes(;loci,seed=genotype_seed),(2,3,1))",
        "size(p.weights)==(1<<p.loci,INPUTS,OUTPUTS)"=>"size(p.weights)==(INPUTS,OUTPUTS,1<<p.loci)",
        "bank[code,col,output]"=>"bank[col,output,code]")
    if bounds
        variant=replace(variant,"module LatticeEcologyPackedStudy"=>"module LatticeEcologyPackedBoundsStudy",
            "for col in 1:INPUTS;value=muladd"=>"@inbounds for col in 1:INPUTS;value=muladd")
    end
    write(joinpath(output,"driver.jl"),read(@__FILE__,String))
    write(joinpath(output,"baseline.jl"),source);vpath=joinpath(output,"packed.jl");write(vpath,variant)
    Base.include(Main,vpath);mod=Base.invokelatest(getproperty,Main,bounds ? :LatticeEcologyPackedBoundsStudy : :LatticeEcologyPackedStudy)
    # Invokelatest enters after the generated module's methods have been defined.
    Base.invokelatest(performance_runs,output,mod,source,variant)
end
function performance_runs(output,mod,source,variant)
    rows=Any[]
    for side in (32,64)
        w=LatticeEcologyStudy.lattice_world(Grid(side,side,1;extent=(side,side,1)))
        p=LatticeEcologyStudy.LatticeEcology(w;occupancy=.5)
        for _ in 1:64;step!(w,p,.1);end
        v,q=variant_state(w,p,mod)
        # Compile both window paths and variant kernels using throwaway copies.
        for (x,y) in (deepcopy((w,p)),deepcopy((v,q)));timed_window(x,y,2);end
        # Every measured window starts from these exact same warmed states.
        checks=Bool[];measurements=Any[]
        for repeat in 1:3
            for arm in ("baseline","packed","packed","baseline")
                x,y=deepcopy(arm=="baseline" ? (w,p) : (v,q));GC.gc()
                elapsed=timed_window(x,y,100)
                push!(measurements,(repeat=repeat,arm=arm,elapsed_seconds=elapsed,voxel_steps_per_second=side^2*100/elapsed))
            end
        end
        a,b=deepcopy((w,p));c,d=deepcopy((v,q))
        allocated_base=@allocated timed_window(a,b,100)
        allocated_packed=@allocated timed_window(c,d,100)
        exact=same_state(b,d)&&a.time==c.time;@assert exact
        serialize(joinpath(output,"state-$(side)-baseline.jls"),(world=a,process=b))
        serialize(joinpath(output,"state-$(side)-packed.jls"),(world=c,process=d))
        for (arm,state) in (("baseline",(w,p)),("packed",(v,q)))
            x,y=deepcopy(state);Profile.clear();Profile.init(n=2_000_000,delay=.0005)
            @profile timed_window(x,y,200)
            open(joinpath(output,"profile-$(side)-$(arm).txt"),"w") do io
                Profile.print(io;format=:flat,sortedby=:count,mincount=3)
            end
            serialize(joinpath(output,"profile-$(side)-$(arm).jls"),Profile.retrieve())
        end
        push!(rows,(side=side,warmup_steps=64,window_steps=100,exact_all_state_arrays=exact,
            population_initial=LatticeEcologyStudy.lattice_statistics(w,p).population,
            population_final=LatticeEcologyStudy.lattice_statistics(a,b).population,
            baseline_bytes_per_step=allocated_base/100,packed_bytes_per_step=allocated_packed/100,
            measurements=measurements))
        open(joinpath(output,"results.json"),"w") do io
            JSON3.write(io,(julia=string(VERSION),threads=Threads.nthreads(),baseline_sha256=bytes2hex(sha256(source)),
                packed_sha256=bytes2hex(sha256(variant)),driver_sha256=bytes2hex(sha256(read(@__FILE__))),
                shared_sources=Interior.arena_sources(),rows=rows))
        end
        println("completed side=",side," exact=",exact)
    end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    length(ARGS) in (1,2)||error("fresh output directory [--bounds]");performance_main(first(ARGS);bounds=length(ARGS)==2&&ARGS[2]=="--bounds")
end
