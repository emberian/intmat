#!/usr/bin/env julia
# Detached frame export needs only Julia stdlibs, not the simulation package.
using Serialization,Printf
length(ARGS)==2 || error("usage: export_device_habitat_frames.jl FRAME_DIR FRESH_OUTPUT")
source,out=ARGS;ispath(out)&&error("fresh output required");mkdir(out)
files=sort(filter(f->occursin(r"^frame-\d+\.jls$",f),readdir(source)))
open(joinpath(out,"frames.json"),"w") do io
    println(io,"{\"field_names\":[\"free_a\",\"free_b\",\"raw_a\",\"raw_b\",\"enzyme_a\",\"enzyme_b\",\"bound_a\",\"bound_b\",\"permeability\"],\"fields\":\"concentration block means; permeability dimensionless\",\"order\":\"Fortran\",\"frames\":[")
    for (i,file) in enumerate(files)
        x=deserialize(joinpath(source,file));name=splitext(file)[1]
        for (suffix,value) in (("fields.f32",x.fields),("positions.f32",x.position),("radii.f32",x.radii),
            ("hidden.f32",x.hidden),("actions.f32",x.actions),("energy.f32",x.energy),("metabolites.f32",x.metabolites))
            open(f->write(f,Float32.(value)),joinpath(out,name*"."*suffix),"w")
        end
        for (suffix,value) in (("ids.i32",x.ids),("lineages.i32",x.lineages),("trace_ids.i32",x.trace_ids),("bond_ids.i32",x.bond_ids))
            open(f->write(f,Int32.(value)),joinpath(out,name*"."*suffix),"w")
        end
        @printf(io,"%s{\"name\":\"%s\",\"time\":%.17g,\"coarse\":%d,\"population\":%d,\"trace_count\":%d,\"hidden\":%d,\"actions\":%d,\"degree\":%d,\"extent\":[%.9g,%.9g,%.9g]}",i==1 ? "" : ",\n",name,x.time,size(x.fields,1),length(x.ids),length(x.trace_ids),size(x.hidden,2),size(x.actions,2),size(x.bond_ids,2),x.extent...)
    end
    println(io,"\n]}")
end
