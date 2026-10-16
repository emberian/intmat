# Complete resident sidecar adapter; numerical laws remain in the study module.
include("../studies/native_living_construction.jl")
const LivingLaw=NativeLivingConstruction
# Includes are literal and relative throughout this study dependency tree.
function living_source_pins(path=joinpath(@__DIR__,"../studies/native_living_construction.jl"),pins=Dict{String,String}())
    path=normpath(path);haskey(pins,path)&&return pins
    source=read(path,String);pins[path]=bytes2hex(sha256(source))
    for m in eachmatch(r"include\(\"([^\"]+)\"\)",source)
        living_source_pins(joinpath(dirname(path),m.captures[1]),pins)
    end
    pins
end
const LIVING_SOURCE_PINS=merge(living_source_pins(),Dict(normpath(@__FILE__)=>bytes2hex(sha256(read(@__FILE__)))))
function living_snapshot(h)
    (;schema="interior-living-sidecar-1",source=copy(LIVING_SOURCE_PINS),state=map(Array,h.state),
        events=map(Array,h.events),original_masks=Array(h.original_masks),options=h.options,
        steps=h.steps,memory=h.memory,inheritance_fraction=h.inheritance_fraction,
        birth_transfer_adjusted=h.birth_transfer_adjusted)
end
function living_restore(p,x)
    x.schema=="interior-living-sidecar-1"&&x.source==LIVING_SOURCE_PINS || error("living sidecar source mismatch")
    # Rebuild ownership explicitly: never deserialize a device array or reset a
    # learner through a fresh constructor. Every queued cohort and clock survives.
    own(a)=copyto!(similar(p.cells.state.inputs,eltype(a),size(a)),a)
    z=map(own,x.state);events=map(own,x.events)
    base=LivingLaw.C.ConstructionHabitat(p.habitat,z,p.cells.backend,x.options,own(x.original_masks),x.steps,LivingLaw.L.NativePolicyLaunch())
    LivingLaw.LivingHabitat(LivingLaw.D.CohortHabitat(LivingLaw.J.JointHabitat(base,x.memory)),events,x.inheritance_fraction,x.birth_transfer_adjusted)
end
living_rebind!(h,p)=(h.habitat=p.habitat;h)
