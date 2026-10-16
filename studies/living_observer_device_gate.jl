#!/usr/bin/env julia
# Bounded device observation gate: no simulation or performance claims.
include(joinpath(@__DIR__,"../scripts/device_cellular_setup.jl"))
using KernelAbstractions, TOML, SHA
include("native_living_construction.jl")
const N=NativeLivingConstruction
length(ARGS)==2 || error("usage: living_observer_device_gate.jl REFERENCE_SOURCE OUTPUT_JSON")
source=read(ARGS[1],String)
a=findfirst("function frame(",source).start;b=findnext("function birth_events",source,a).start
Base.include_string(N,replace(source[a:prevind(source,b)],"function frame("=>"function frame_reference("),"observer_reference")
f=N.fixture(joinpath(@__DIR__,"../experiments/native/living-patch-recruitment-16.toml");backend)
h=N.sidecar(f.p;dt=.02,horizon=2.,act_time=100.,inheritance_fraction=.5)
fill!(h.state.learned,.125f0);fill!(h.state.paid,.5);fill!(h.state.requested,.25)
s=f.p.cells.state;residents=findall(==(1),Array(s.active));i,j=residents[1:2]
bonds=Array(s.bonds);ids=Array(s.ids);bond_ids=Array(s.bond_ids);rest=Array(s.rest)
bonds[i,1]=j;bond_ids[i,1]=ids[j];rest[i,1]=.125f0
copyto!(s.bonds,bonds);copyto!(s.bond_ids,bond_ids);copyto!(s.rest,rest)
Interior.device_cellular_status!(f.p.cells)
origin=N.C.A.audit(f.p)
groups=(f.p.cells.state,f.p.cells.environment,f.p.habitat.arrays,f.p.habitat.feedstock.arrays,h.state)
originals=[map(Array,g) for g in groups]
old=N.frame_reference(f.w,f.p,h,origin);new=N.frame(f.w,f.p,h,origin)
@assert typeof(old)==typeof(new)
@assert isequal(old,new)
@assert JSON3.write(old)==JSON3.write(new)
@assert !isempty(new.bonds)
for (group,original) in zip(groups,originals),key in keys(group)
    @assert isequal(Array(group[key]),original[key]) "observer changed source state"
end
digest=bytes2hex(sha256(Vector{UInt8}(codeunits(JSON3.write(new)))))
for group in groups,array in values(group);fill!(array,zero(eltype(array)));end
KernelAbstractions.synchronize(backend)
@assert bytes2hex(sha256(Vector{UInt8}(codeunits(JSON3.write(new)))))==digest "captured frame aliases source state"
result=(schema="living-observer-device-gate-v1",runtime=device_cellular_runtime(),side=16,
    exact_record=true,exact_json=true,record_types_equal=true,source_unchanged=true,mutation_independent=true,
    cells=length(new.cells),bonds=length(new.bonds),json_sha256=digest,
    observer_source_sha256=bytes2hex(sha256(read(joinpath(@__DIR__,"native_living_construction.jl")))),
    reference_source_sha256=bytes2hex(sha256(Vector{UInt8}(codeunits(source)))))
open(io->JSON3.write(io,result),ARGS[2],"w");println(JSON3.write(result))
