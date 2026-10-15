#!/usr/bin/env julia
# Optical observation only: no simulation, spatial refinement or invented detail.
using JSON3, SHA
include(joinpath(@__DIR__,"..","src","render.jl"))
length(ARGS)==2 || error("ASSAY_DIRECTORY FRESH_OUTPUT_DIRECTORY")
source,out=ARGS
ispath(out) && error("fresh output required")
mkpath(out)
names=("none_sham","neural_sham")
frames=[JSON3.read(read(joinpath(source,name*"-final-field.json"),String)) for name in names]
shape=Tuple(Int.(frames[1].shape))
all(x->Tuple(Int.(x.shape))==shape,frames) || error("different shapes")
# One shared monotone transfer for both branches. Scale from the joint maximum,
# record it explicitly; neither image is independently normalized.
free_scale=maximum(maximum(Float64.(x.fields[2])) for x in frames)
bound_scale=maximum(maximum(Float64.(x.fields[8])) for x in frames)
free_scale>0 && bound_scale>0 || error("empty optical channels")
settings=RenderSettings(width=640,height=640,steps=128,orbit=.6,elevation=.42,
    absorption=1.6,emission=5.,scattering=.1,exposure=1.4)
records=[]
for (name,x) in zip(names,frames)
    free=reshape(Float32.(x.fields[2]),shape)./Float32(free_scale)
    bound=reshape(Float32.(x.fields[8]),shape)./Float32(bound_scale)
    # Warm absorption follows retained B. Emission is localized to that actual
    # retained field, with gain modulated by free B; uniform free food otherwise
    # obscures the spatial memory under a glowing box.
    fields=(density=bound,structure=bound,activity=bound.*sqrt.(free))
    file=joinpath(out,name*".ppm")
    seconds=@elapsed render_volume(fields,file;settings)
    push!(records,(;name,time=x.time,file=basename(file),seconds,
        source_sha256=bytes2hex(sha256(read(joinpath(source,name*"-final-field.json")))),
        image_sha256=bytes2hex(sha256(read(file)))))
end
manifest=(;kind="optical comparison of recorded native fields",shape,
    mapping="retained B sets warm absorption/density; cyan emission is normalized bound B times sqrt(normalized free B); shared scales; full periodic box observed as a cut specimen; no cell glyphs",
    free_scale,bound_scale,settings=(; (k=>getproperty(settings,k) for k in fieldnames(RenderSettings))...),
    renderer_sha256=bytes2hex(sha256(read(joinpath(@__DIR__,"..","src","render.jl")))),
    driver_sha256=bytes2hex(sha256(read(@__FILE__))),records)
open(io->JSON3.pretty(io,manifest),joinpath(out,"optical.json"),"w")
