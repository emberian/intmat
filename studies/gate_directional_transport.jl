# CPU-only scoped numerical/performance gate; no production edits or GPU use.
using Interior
using Interior: JSON3,SHA
out=only(ARGS);ispath(out)&&error("fresh output required")
mkdir(out)
source=joinpath(@__DIR__,"pair_transport_ecology.jl")
source_hash=bytes2hex(SHA.sha256(read(source)))
cp(source,joinpath(out,"study.jl"))
if get(ENV,"TRANSPORT_RUN_TESTS","1")=="1"
    include(joinpath(@__DIR__,"..","test","pair_transport_ecology.jl"))
else
    include(source);using .PairTransportEcology
end
size=parse(Int,get(ENV,"TRANSPORT_SIZE","24"));pairs=parse(Int,get(ENV,"TRANSPORT_PAIRS","3"))
size in (24,32,64)&&1<=pairs<=3||error("bounded CPU gate")
strategies=size==64 ? (:explicit,:lines) : (:explicit,:pairs,:lines)
runs=NamedTuple[]
for pair in 1:pairs
    x=PairTransportEcology.compare_transport(;size,founders=8,capacity=64,diffusion=20.,dt=.01,steps=4,cycles=(1,2,4),strategies)
    push!(runs,x);println("TRANSPORT_PAIR_COMPLETE ",pair);flush(stdout)
end
result=(source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
    study_sha256=source_hash,driver_sha256=bytes2hex(SHA.sha256(read(@__FILE__))),
    runtime=(julia=string(VERSION),threads=Threads.nthreads(),backend="CPU"),results=runs)
open(io->JSON3.pretty(io,result),joinpath(out,"results.json"),"w")
for x in runs,r in x.rows;println((r.mode,r.cycles,r.seconds,r.field_relative_l2,r.raw_relative_l2,r.bound_relative_l2));end
