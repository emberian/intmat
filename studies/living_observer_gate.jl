#!/usr/bin/env julia
# Exact observer parity and complete capture/serialization timing; no simulation benchmark.
using Interior, KernelAbstractions, TOML, SHA, Statistics
using Interior: JSON3
include("native_living_construction.jl")
const N=NativeLivingConstruction
const BASELINE="1ff4acea0673b9697c0dd9e543f0c8638c0a570b"
source=read(`git show $BASELINE:studies/native_living_construction.jl`,String)
a=findfirst("function frame(",source).start;b=findnext("function birth_events",source,a).start
Base.include_string(N,replace(source[a:prevind(source,b)],"function frame("=>"function frame_reference("),"observer_reference")
results=Any[]
for side in parse.(Int, isempty(ARGS) ? ["16","32"] : ARGS)
    cfg=TOML.parsefile(joinpath(@__DIR__,"../experiments/native/living-patch-recruitment-16.toml"))
    cfg["grid"]["shape"]=[side,side,side]
    path=tempname()*".toml";open(io->TOML.print(io,cfg),path,"w")
    f=N.fixture(path);h=N.sidecar(f.p;dt=.02,horizon=2.,act_time=100.,inheritance_fraction=.5)
    origin=N.C.A.audit(f.p)
    # Nonzero private fields exercise selected transfer, reduction and learned norms.
    h.state.learned .= .125f0;h.state.requested .= .25;h.state.paid .= .5
    # Exercise native edge extraction on an ID-valid pair without changing mass/energy.
    state=f.p.cells.state;residents=findall(==(1),state.active)
    if length(residents)>=2 && size(state.bonds,2)>0
        i,j=residents[1:2];state.bonds[i,1]=j;state.bond_ids[i,1]=state.ids[j];state.rest[i,1]=.125f0
    end
    old=N.frame_reference(f.w,f.p,h,origin);new=N.frame(f.w,f.p,h,origin)
    encode(x)=JSON3.write(x)
    @assert isequal(old,new)
    @assert encode(old)==encode(new)
    digest=bytes2hex(sha256(Vector{UInt8}(codeunits(encode(new)))))
    # Mutate every physical and learner array after capture, including native edges.
    for group in (f.p.cells.state,f.p.cells.environment,f.p.habitat.arrays,f.p.habitat.feedstock.arrays,h.state)
        for array in values(group)
            array isa AbstractArray && fill!(array,zero(eltype(array)))
        end
    end
    @assert bytes2hex(sha256(Vector{UInt8}(codeunits(encode(new)))))==digest "captured frame aliases live state"
    # Recreate a valid static state, origin once, excluding initialization from timing.
    f=N.fixture(path);h=N.sidecar(f.p;dt=.02,horizon=2.,act_time=100.,inheritance_fraction=.5);origin=N.C.A.audit(f.p)
    rows=Dict{String,Any}()
    for serialize in (false,true)
        capture(fun)=serialize ? JSON3.write(fun(f.w,f.p,h,origin)) : fun(f.w,f.p,h,origin)
        capture(N.frame_reference);capture(N.frame)
        samples=Dict("reference"=>Any[],"optimized"=>Any[])
        for repeat in 1:9
            for (name,fun) in (isodd(repeat) ? (("reference",N.frame_reference),("optimized",N.frame)) : (("optimized",N.frame),("reference",N.frame_reference)))
                m=@timed capture(fun)
                push!(samples[name],(;seconds=m.time,bytes=m.bytes,gc_seconds=m.gctime))
            end
        end
        rows[serialize ? "frame_and_json" : "frame"]=Dict(k=>(;median_seconds=median(x.seconds for x in v),median_bytes=median(x.bytes for x in v),samples=v) for (k,v) in samples)
    end
    push!(results,(;side,exact_record=true,exact_json=true,mutation_independent=true,json_sha256=digest,measurements=rows))
    println(stderr,"observer gate completed ",side,"³")
end
JSON3.write(stdout,(;schema="living-observer-gate-v1",baseline=BASELINE,threads=Threads.nthreads(),results));println()
