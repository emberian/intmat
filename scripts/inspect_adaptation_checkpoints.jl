#!/usr/bin/env julia
# Read-only extraction of preregistered counters retained in trusted local checkpoints.
using Interior, Serialization, JSON3
out=abspath(first(ARGS)); rows=Any[]
for (dir,_,files) in walkdir(out)
    "terminal.jls" in files || continue
    saved=deserialize(joinpath(dir,"terminal.jls"));a,phase,m,p,b=saved.state;c=p.controller
    groups=occursin("assay-",dir)||occursin("life-only-",dir)
    push!(rows,(path=relpath(dir,out),time=a.world.time,population=length(c.ids),
        blocked_births=c.blocked_births,births=c.births,deaths=c.deaths,
        group_stored_energy=groups ? [sum((Float64(c.energy[i]) for i in eachindex(c.ids) if c.lineages[i]==g);init=0.) for g in 1:2] : nothing,
        group_stored_material=groups ? [sum((Float64(c.storage[i]) for i in eachindex(c.ids) if c.lineages[i]==g);init=0.) for g in 1:2] : nothing,
        original_source_hashes=saved.sources))
end
JSON3.write(stdout,rows);println()
