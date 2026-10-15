#!/usr/bin/env julia
# CPU-only matched full mechanical steps from actual repair-weave snapshots.
using Interior, LinearAlgebra, JSON3, TOML, SHA, Serialization
include(joinpath(@__DIR__,"..","test","fixtures","bond_network_reference.jl"))
function weave_fixture(root,step)
    report=JSON3.read(read(joinpath(root,"run.json"),String));c=report.params
    w=World(Grid(1,1,1;extent=(2,2,2));seed=911)
    p=BondNetwork(w;size=Tuple(Int.(c.dimensions)),extent=(1.6,1.35,.075),stiffness=5.,mass=c.node_mass,
        damping=7.,failure_strain=.17,weak_seam=.22,jitter=.10,stretch_rate=0.,shear_rate=0.)
    frame=JSON3.read(read(joinpath(root,"frames",lpad(step,5,'0'),"geometry.json"),String))
    firstframe=JSON3.read(read(joinpath(root,"frames","00000","geometry.json"),String))
    # JSON3 stores the exported N×3 / B×2 matrices in flat column-major order.
    matrix(values,T,width)=permutedims(reshape(T.(values),:,width))
    original=length(p.bond_ids)
    p.position=matrix(frame.positions,Float64,3);p.velocity=matrix(frame.velocities,Float64,3)
    p.reference=matrix(frame.reference,Float64,3);p.node_ids=Int.(frame.node_ids)
    p.bonds=matrix(frame.bonds,Int,2);p.bond_ids=Int.(frame.bond_ids)
    p.rest_length=Float64.(frame.rest_length);p.stiffness=Float64.(frame.stiffness)
    append!(p.threshold,fill(.22,length(p.bond_ids)-original))
    p.damage=Float64.(frame.damage);p.peak_strain=Float64.(frame.peak_strain);p.elapsed=Float64(frame.time)
    history=only(x for x in report.history if x.step==step)
    p.damage_energy=history.energy.damage;p.damping_energy=history.energy.damping;p.external_work=history.energy.external_work
    w.time=p.elapsed
    # Snapshot reconstruction is accepted only if every energy component matches.
    actual=mechanical_energy(p)
    all(isequal(getproperty(actual,k),Float64(getproperty(history.energy,k))) for k in keys(actual)) || error("Snapshot energy mismatch")
    w,p,matrix(firstframe.reference,Float64,3),c,original
end
function drive!(w,p,initial_reference,c,advance_fn)
    dt=Float64(c.dt);t=p.elapsed+dt;phase=2pi*t/c.cycle_time
    extension=.30*(1-cos(phase))/2+.063t/c.cycle_time;shear=.16sin(phase)
    for q in eachindex(p.grip)
        p.grip[q]==0 && continue
        angle=p.grip[q]*c.grip_twist*sin(phase)
        y,z=initial_reference[2,q],initial_reference[3,q]
        p.reference[2,q]=cos(angle)*y-sin(angle)*z;p.reference[3,q]=sin(angle)*y+cos(angle)*z
    end
    p.stretch_rate=extension/t;p.shear_rate=shear/t
    advance_fn(w,p,dt);w.time=p.elapsed
end
function run_steps!(w,p,reference,c,f,n)
    for _ in 1:n;drive!(w,p,reference,c,f);end
end
function state_sha(p)
    io=IOBuffer();serialize(io,p);bytes2hex(sha256(take!(io)))
end
function compare_case(root,step)
    w,p,reference,c,original=weave_fixture(root,step)
    for f in (BondNetworkReference.advance!,Interior.advance!)
        run_steps!(deepcopy(w),deepcopy(p),reference,c,f,3)
    end
    samples=Dict("reference"=>Any[],"candidate"=>Any[])
    for _ in 1:5,(name,f) in (("reference",BondNetworkReference.advance!),("candidate",Interior.advance!))
        a=deepcopy(p);v=deepcopy(w);GC.gc()
        t=@timed run_steps!(v,a,reference,c,f,20)
        push!(samples[name],(;seconds=t.time,bytes=t.bytes,gc_seconds=t.gctime,steps=20))
    end
    a=deepcopy(p);b=deepcopy(p)
    run_steps!(deepcopy(w),a,reference,c,BondNetworkReference.advance!,100)
    run_steps!(deepcopy(w),b,reference,c,Interior.advance!,100)
    sha_a=state_sha(a);sha_b=state_sha(b);sha_a==sha_b || error("Whole trajectory bit/alias mismatch")
    (;snapshot_step=step,nodes=length(p.node_ids),bonds=length(p.bond_ids),broken=count(==(1.),p.damage),
      repaired=length(p.bond_ids)-original,
      samples,gate=(steps=100,exact_serialized_state=true,reference_sha=sha_a,candidate_sha=sha_b),
      note="Timed steps include validation, tangent bound, grip motion, two force evaluations, damage, damping and work ledgers. Existing repaired bonds are retained; no new repair is inserted inside this fixed mechanical-step benchmark.")
end
length(ARGS)==2 || error("Usage: benchmark_bond_network.jl PILOT_OUTPUT RESULTS_JSON")
result=(;runtime=string(VERSION),threads=Threads.nthreads(),blas_threads=BLAS.get_num_threads(),
    source_sha=bytes2hex(sha256(read(joinpath(@__DIR__,"..","src","mechanics.jl")))),
    reference_sha=bytes2hex(sha256(read(joinpath(@__DIR__,"..","test","fixtures","bond_network_reference.jl")))),
    cases=[compare_case(ARGS[1],step) for step in (0,6000)])
open(io->JSON3.pretty(io,result),ARGS[2],"w");println(JSON3.write(result))
