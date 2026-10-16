#!/usr/bin/env julia
# Exact direct-kernel edge fixtures; positive dt is the public-step precondition.
include(joinpath(@__DIR__,"..","scripts","device_cellular_setup.jl"))
using KernelAbstractions
const EDGE_BASELINE="9c9f1a51f277f66949da00fdb63c7a20c2fcaed5"
source=haskey(ENV,"LIVING_BASELINE_KERNEL") ? read(ENV["LIVING_BASELINE_KERNEL"],String) : read(`git show $(EDGE_BASELINE):src/device_habitat_kernels.jl`,String)
body=split(split(source,"@kernel function dh_react_sorb!")[2],"@kernel function dh_validate_complete!")[1]
Base.include_string(Interior,"@kernel function dh_react_sorb_edge_baseline!"*body,"edge-baseline-reaction.jl")
function edge_host(case)
    a=(enzyme_bits=zeros(UInt32,1,2),bound=zeros(Float32,1,2),release=fill(7f0,1,2),
        grid_ledger=reshape(Float32.(1:12),1,12),grid_ledger_low=fill(1f-6,1,12),
        grid_inert=Float32[1],grid_inert_low=Float32[1f-6])
    feed=(fields=reshape([.25,.4,0.],1,3),)
    e=(fields=reshape([.1,.2,.3],1,3),status=Int32[0])
    if case==:bound_negative_zero;a.bound[1,1]=-0f0
    elseif case==:enzyme_negative_zero;a.enzyme_bits[1,1]=reinterpret(UInt32,-0f0)
    elseif case==:positive_bound;a.bound[1,1]=1f-4
    elseif case==:positive_enzyme;a.enzyme_bits[1,1]=reinterpret(UInt32,1f-4)
    elseif startswith(string(case),"raw_")
        feed.fields[1,1]=endswith(string(case),"nan") ? NaN : endswith(string(case),"inf") ? Inf : -0.
    elseif startswith(string(case),"food_")
        e.fields[1,1]=endswith(string(case),"nan") ? NaN : endswith(string(case),"inf") ? Inf : -0.
    elseif startswith(string(case),"waste_")
        e.fields[1,3]=endswith(string(case),"nan") ? NaN : endswith(string(case),"inf") ? Inf : -0.
    end
    (;a,feed,e)
end
function edge_upload(x)
    map(x) do array
        dst=KernelAbstractions.zeros(backend,eltype(array),size(array)...);copyto!(dst,array);dst
    end
end
function run_edges(out)
    cases=(:empty_compensated,:raw_nan,:raw_inf,:raw_negative_zero,:food_nan,:food_inf,:food_negative_zero,
        :waste_nan,:waste_inf,:waste_negative_zero,:bound_negative_zero,:enzyme_negative_zero,:positive_bound,:positive_enzyme)
    q=(catalytic_rate=6f0,half_saturation=.2f0,food_potentials=(4f0,4f0,0f0),decay_rate=.05f0,
        capacity_per_enzyme=6f0,cross_binding=0f0,adsorption_rate=4f0,desorption_rate=.12f0)
    fq=(potentials=(4f0,4f0,0f0),);g=(volume=1f0,);rows=Any[];checks=0
    for case in cases
        host=edge_host(case);states=Any[]
        for kernel in (Interior.dh_react_sorb_edge_baseline!,Interior.dh_react_sorb!)
            a=edge_upload(host.a);feed=edge_upload(host.feed);e=edge_upload(host.e)
            kernel(backend,128)(a,feed,e,q,fq,g,.02f0;ndrange=1);KernelAbstractions.synchronize(backend)
            push!(states,(a=map(Array,a),feed=map(Array,feed),e=map(Array,e)))
        end
        for group in (:a,:feed,:e),key in keys(states[1][group])
            isequal(states[1][group][key],states[2][group][key]) || error("edge mismatch $case $group.$key")
            checks+=1
        end
        normalized=!isequal(states[1].a.grid_ledger,host.a.grid_ledger)
        if case==:empty_compensated
            normalized || error("fixture failed to exercise compensated normalization")
            states[2].a.release==zeros(Float32,1,2)||error("release reset missing")
        end
        push!(rows,(;case,status=only(states[2].e.status),normalized,release=vec(states[2].a.release)))
    end
    # NaN release values are expected diagnostics; do not serialize non-JSON floats.
    publicrows=[(;case=r.case,status=r.status,normalized=r.normalized,release_finite=all(isfinite,r.release)) for r in rows]
    report=(;schema="living-chemistry-edge-gate-v1",runtime=device_cellular_runtime(),baseline=EDGE_BASELINE,
        kernel_sha256=bytes2hex(Interior.SHA.sha256(read(joinpath(@__DIR__,"..","src","device_habitat_kernels.jl")))),checks,cases=publicrows)
    open(io->JSON3.write(io,report),out,"w");println((;backend=kind,cases=length(rows),checks))
end
run_edges(ARGS[1])
