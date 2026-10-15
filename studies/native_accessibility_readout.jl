#!/usr/bin/env julia
# Read-only actual-frame spatial accounting. No checkpoint restore or simulation.
using Serialization, SHA, Statistics, JSON3
length(ARGS)==3||error("MIRROR_ROOT NATIVE_FRAME_RELATIVE_PATH FRESH_JSON required")
root,relative,out=ARGS;ispath(out)&&error("fresh output required")
manifest=JSON3.read(read(joinpath(root,"run.json"),String));path=joinpath(root,relative)
entry=only(filter(v->v.file==relative,manifest.frames))
hash=bytes2hex(sha256(read(path)));hash==entry.sha256||error("recorded frame SHA differs")
x=deserialize(path);x.shape==x.native_shape||error("requires actual native fields")
dims=Int.(x.shape);extent=Float64.(x.extent);spacing=extent./dims;V=prod(dims);volume=prod(spacing)
field=reshape(x.fields,V,9);R=Float64(manifest.config.cellular.exchange_radius)
covered=falses(V);percell=Any[]
for i in eachindex(x.ids)
    pos=Float64.(x.position[i,:]);lo=floor.(Int,(pos.-R)./collect(spacing).+.5)
    span=ceil.(Int,2R./collect(spacing)).+3;total=0.;values=zeros(9)
    for dz in 0:span[3]-1,dy in 0:span[2]-1,dx in 0:span[1]-1
        index=mod1.(lo.+[dx,dy,dz],collect(dims));r2=0.
        for d in 1:3
            delta=abs((index[d]-.5)*spacing[d]-pos[d]);delta=min(delta,extent[d]-delta);r2+=delta^2
        end
        if r2<R^2
            weight=(1-r2/R^2)^2;site=index[1]+dims[1]*((index[2]-1)+dims[2]*(index[3]-1))
            covered[site]=true;total+=weight
            for k in 1:9;values[k]+=weight*Float64(field[site,k]);end
        end
    end
    values./=total
    push!(percell,(id=x.ids[i],lineage=x.lineages[i],mass=x.mass[i],energy=x.energy[i],
        free=values[1:2],raw=values[3:4],enzyme=values[5:6],bound=values[7:8],combined_permeability=values[9]))
end
enzyme=Float64.(field[:,5]).+field[:,6];constructed=enzyme.>0
raw=Float64.(field[:,3]).+field[:,4];free=Float64.(field[:,1]).+field[:,2];bound=Float64.(field[:,7]).+field[:,8]
cfg=manifest.config;source_total=sum(Float64.(cfg.feedstock.source_rates));eff=.65;potential=4.
release=zeros(2)
for k in 1:2,site in 1:V
    release[k]+=cfg.catalysis.catalytic_rate*Float64(field[site,4+k])*Float64(field[site,2+k])/(cfg.catalysis.half_saturation+Float64(field[site,2+k]))*volume
end
enzyme_mass=[sum(Float64,field[:,4+k])*volume for k in 1:2]
substrate_limited_share=[sum(Float64(field[site,4+k]) for site in 1:V if field[site,2+k]<cfg.catalysis.half_saturation)*volume/max(enzyme_mass[k],eps()) for k in 1:2]
summary=(schema="native-spatial-accessibility-readout-v1",time=x.time,population=length(x.ids),frame=relative,frame_sha256=hash,
    source_commit=manifest.source_commit,driver_sha256=bytes2hex(sha256(read(@__FILE__))),
    scope="actual full-resolution frame; support readout in Float64; no uptake, future behavior or cumulative lineage attribution inferred",
    voxel_spacing=spacing,support_radius=R,current_support_union_fraction=count(covered)/V,any_positive_enzyme_fraction=count(constructed)/V,
    total_raw=sum(raw)*volume,total_free=sum(free)*volume,total_bound=sum(bound)*volume,total_enzyme=sum(enzyme)*volume,
    raw_inside_current_support=sum(raw[covered])*volume,enzyme_inside_current_support=sum(enzyme[covered])*volume,
    raw_inside_any_positive_enzyme=sum(raw[constructed])*volume,
    local_raw_median=median([sum(v.raw) for v in percell]),local_free_median=median([sum(v.free) for v in percell]),
    local_bound_median=median([sum(v.bound) for v in percell]),local_combined_permeability_median=median([v.combined_permeability for v in percell]),
    mean_energy=mean(x.energy),basal_work_rate=Float64(cfg.cellular.basal_rate)*sum(Float64,x.mass),
    estimated_instantaneous_release=release,
    enzyme_mass,enzyme_weighted_substrate_saturation=release./(Float64(cfg.catalysis.catalytic_rate).*enzyme_mass),
    enzyme_mass_share_below_half_saturation=substrate_limited_share,
    optimistic_ongoing_work_at_current_support=source_total*eff*potential*count(covered)/V,
    optimistic_ongoing_work_at_any_positive_enzyme=source_total*eff*potential*count(constructed)/V,
    bound_note="ongoing-work bounds use the opposed masks' constant summed source density; ignore transport, unequal species, accessibility, all other costs, and kinetic insufficiency of tiny positive enzyme",
    cells=percell)
open(io->JSON3.pretty(io,summary),out,"w")
println((time=x.time,population=summary.population,support=summary.current_support_union_fraction,constructed=summary.any_positive_enzyme_fraction,
    raw=summary.total_raw,free=summary.total_free,release,basal=summary.basal_work_rate,local_raw=summary.local_raw_median,
    local_free=summary.local_free_median,source_work_upper=summary.optimistic_ongoing_work_at_current_support))
