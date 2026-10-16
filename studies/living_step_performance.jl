#!/usr/bin/env julia
#= Exact full-step comparison against the pre-shortcut chemistry, CPU only.
Run JULIA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 julia --project=. --compiled-modules=existing studies/living_step_performance.jl [output.json].
The pinned Git baseline supplies only dh_react_sorb!; all other code is current.
=#
using Interior, KernelAbstractions, TOML, Statistics
using Interior: JSON3
include("native_living_construction.jl")
const N=NativeLivingConstruction
const BASELINE="9c9f1a51f277f66949da00fdb63c7a20c2fcaed5"
source=haskey(ENV,"LIVING_BASELINE_KERNEL") ? read(ENV["LIVING_BASELINE_KERNEL"],String) : read(`git show $(BASELINE):src/device_habitat_kernels.jl`,String)
body=split(split(source,"@kernel function dh_react_sorb!")[2],"@kernel function dh_validate_complete!")[1]
Base.include_string(Interior,"@kernel function dh_react_sorb_baseline!"*body,"baseline-reaction.jl")
@eval Interior begin
    const living_probe_baseline=Ref(false)
    function device_anabolic_finish!(h::DeviceCatalyticHabitat,w,p,dt)
        V=prod(Int.(p.geometry.dims))
        reaction=living_probe_baseline[] ? dh_react_sorb_baseline! : dh_react_sorb!
        device_cellular_launch!(p,reaction,h.arrays,h.feedstock.arrays,p.environment,h.parameters,h.feedstock.parameters,p.geometry,dt;ndrange=V)
        device_cellular_launch!(p,dh_permeability!,h.arrays,p.environment,p.geometry,h.parameters;ndrange=V)
        device_cellular_launch!(p,dh_validate_complete!,h.arrays,h.feedstock.arrays,p.environment,p.geometry;ndrange=max(V,Int(p.geometry.capacity)))
        nothing
    end
end
function arrays(x)
    (;cells=map(copy,x.p.cells.state),environment=map(copy,x.p.cells.environment),
      habitat=map(copy,x.p.habitat.arrays),feed=map(copy,x.p.habitat.feedstock.arrays),
      private=map(copy,x.h.state),events=map(copy,x.h.events))
end
function compare(a,b)
    count=0
    for group in keys(a),key in keys(a[group])
        isequal(a[group][key],b[group][key])||error("full-state mismatch $group.$key")
        count+=length(a[group][key])
    end
    count
end
function advance(x,baseline,n)
    Interior.living_probe_baseline[]=baseline
    N.L.advance!(x.w,x.p,x.h,n*.02;dt=.02)
end
function run(out)
    rows=Any[]
    for resolution in (16,32,64),founders in (24,46)
        cfg=TOML.parsefile("experiments/native/living-patch-recruitment-16.toml")
        cfg["grid"]["shape"]=fill(resolution,3);cfg["cellular"]["founders"]=founders
        # Extra founders are an explicitly separate dense-contact performance fixture.
        recipe=tempname()*".toml";open(io->TOML.print(io,cfg),recipe,"w")
        println((;stage="fixture",resolution,founders));flush(stdout)
        f=N.fixture(recipe);h=N.sidecar(f.p;dt=.02,horizon=2.,act_time=100.,inheritance_fraction=.5)
        original=(w=f.w,p=f.p,h=h)
        println("warming living fixture");flush(stdout)
        advance(original,false,220)
        # The second fixture spreads catalyst to every voxel (paid law unchanged
        # thereafter), making the empty shortcut unavailable as a negative control.
        for coverage in (:native,:full)
            start=deepcopy(original)
            if coverage==:full
                bits=start.p.habitat.arrays.enzyme_bits
                reinterpret(Float32,bits) .= max.(reinterpret(Float32,bits),1f-7)
            end
            println((;stage="parity",resolution,founders,coverage));flush(stdout)
            a=deepcopy(start);b=deepcopy(start)
            advance(a,true,2);advance(b,false,2)
            checked=compare(arrays(a),arrays(b))
            times=Dict(true=>Float64[],false=>Float64[])
            for repeat in 1:5
                # Matched states and alternating execution order control temporal drift.
                a=deepcopy(start);b=deepcopy(start)
                for baseline in (isodd(repeat) ? (true,false) : (false,true))
                    x=baseline ? a : b
                    GC.gc() # exclude clone-triggered collection from step timing
                    push!(times[baseline],@elapsed advance(x,baseline,30))
                end
                checked+=compare(arrays(a),arrays(b))
            end
            occupied=count(i->start.p.habitat.arrays.enzyme_bits[i,1]!=0||start.p.habitat.arrays.enzyme_bits[i,2]!=0||start.p.habitat.arrays.bound[i,1]!=0||start.p.habitat.arrays.bound[i,2]!=0,axes(start.p.habitat.arrays.bound,1))
            origin=N.C.A.audit(b.p)
            N.frame(b.w,b.p,b.h,origin) # warm observer separately
            observer=@elapsed N.frame(b.w,b.p,b.h,origin)
            push!(rows,(;resolution,founders,coverage,occupied_fraction=occupied/resolution^3,
                baseline_step_seconds=median(times[true])/30,optimized_step_seconds=median(times[false])/30,
                speedup=median(times[true])/median(times[false]),checked_array_elements=checked,observer_seconds=observer,
                population=count(==(1),b.p.cells.state.active),baseline_samples=times[true],optimized_samples=times[false]))
            open(io->JSON3.write(io,(;baseline=BASELINE,kernel_sha256=bytes2hex(N.sha256(read(joinpath(@__DIR__,"..","src","device_habitat_kernels.jl")))),threads=Threads.nthreads(),rows)),out,"w")
            println(last(rows));flush(stdout)
        end
    end
    open(io->JSON3.write(io,(;baseline=BASELINE,kernel_sha256=bytes2hex(N.sha256(read(joinpath(@__DIR__,"..","src","device_habitat_kernels.jl")))),threads=Threads.nthreads(),rows)),out,"w")
end
abspath(PROGRAM_FILE)==(@__FILE__) && run(length(ARGS)>0 ? ARGS[1] : "/tmp/living-step-performance.json")
