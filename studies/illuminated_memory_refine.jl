include("illuminated_memory.jl")
using LinearAlgebra
function refine(root,grid)
    r=JSON3.read(read(joinpath(root,"run.json"),String),Dict{String,Any});c=r["params"]
    solid=loadfield(root,last(r["frames"]),"solid",r["shape"][1]);previous=nothing;records=Any[]
    for ds in (2c["ray_ds"],c["ray_ds"],c["ray_ds"]/2)
        elapsed=@elapsed q,l,e=absorbed_power(solid;grid=grid,ds=ds,alpha=c["alpha"],power=c["power"],contrast=c["contrast"])
        if previous!==nothing
            pq,pl,pe=previous;errors=[norm(collect(a).-collect(b)) for (a,b) in zip(e,pe)]
            push!(records,(ds=ds,relative_dose=norm(q.-pq)/norm(pq),relative_absorption=abs(l.absorbed-pl.absorbed)/pl.absorbed,max_endpoint=maximum(errors),mean_endpoint=mean(errors),p95_endpoint=quantile(errors,.95),elapsed_seconds=elapsed,energy=l))
        end
        previous=(q,l,e)
    end
    out=joinpath(root,"native-ray-refinement.json");ispath(out)&&error("Refusing existing refinement")
    open(io->JSON3.pretty(io,(grid=grid,records=records,source_sha256=bytes2hex(sha256(read(@__FILE__))),producer_sha256=bytes2hex(sha256(read(joinpath(@__DIR__,"illuminated_memory.jl")))),solid_sha256=bytes2hex(sha256(reinterpret(UInt8,vec(solid)))))),out,"w")
    cp(@__FILE__,joinpath(root,"native-ray-refinement-source.jl"));println("REFINEMENT ",JSON3.write(records))
end
refine(ARGS[1],parse(Int,ARGS[2]))
