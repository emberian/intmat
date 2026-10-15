# Experimental full-physiology carrier. No production methods are replaced.
mutable struct RKCMechanicalPolicy
 splits::Int
 max_evaluations::Int
 force_evaluations::Int
 calls::Int
end
RKCMechanicalPolicy(;splits=16,max_evaluations=4096)=begin
 splits>0&&max_evaluations>0||throw(ArgumentError("positive RKC policy controls required"))
 RKCMechanicalPolicy(splits,max_evaluations,0,0)
end
function (policy::RKCMechanicalPolicy)(grid,ids,positions,radii,bonds,dt;drag_multiplier=nothing,kwargs...)
 drag_multiplier===nothing||throw(ArgumentError("experimental RKC currently supports unit mechanical drag only"))
 cache=_CellularRKCCache()
 result=cellular_rkc_funded_mechanics!(grid,ids,positions,radii,bonds,dt;kwargs...,
  relaxer=cache,max_step=Float64(dt)/policy.splits,max_evaluations=policy.max_evaluations)
 policy.force_evaluations+=cache.force_evaluations;policy.calls+=1
 result
end
function install_cellular_rkc_ecology!()
 root=dirname(@__DIR__);loaded=arena_sources()
 original=Dict{String,String}();generated=Dict{String,String}()
 for file in ("cellular.jl","constructive_cells.jl")
  text=read(joinpath(root,"src",file),String)
  bytes2hex(SHA.sha256(text))==loaded[file]||error("loaded source differs from disk: $file")
  original[file]=text
 end
 text=original["cellular.jl"]
 first=findfirst("function _ecology_step_storage!(",text);last=findfirst("# The ordinary law remains",text)
 first!==nothing&&last!==nothing||error("cellular step source markers changed")
 body=text[Base.first(first):Base.first(last)-1]
 header="function _ecology_step_storage!(w::World,p::CellularEcology,dt::Real,q,weights,marks,anabolism)"
 length(findall(header,body))==1&&length(findall("p.mechanical=cellular_mechanics!(",body))==1||error("cellular mechanics hook changed")
 body=replace(body,header=>"function _rkc_ecology_step_storage!(w::World,p::CellularEcology,dt::Real,q,weights,marks,anabolism,mechanical_policy)",
  "p.mechanical=cellular_mechanics!("=>"p.mechanical=mechanical_policy(")
 generated["cellular.jl"]=body
 Base.include_string(@__MODULE__,body,"cellular.jl:rkc-study")
 text=original["constructive_cells.jl"]
 first=findfirst("function advance!(w::World,p::ConstructiveEcology,dt::Real)",text);last=findfirst("anabolic_statistics(p::ConstructiveEcology)",text)
 first!==nothing&&last!==nothing||error("constructive step source markers changed")
 body=text[Base.first(first):Base.first(last)-1]
 body=replace(body,"function advance!(w::World,p::ConstructiveEcology,dt::Real)"=>"function _rkc_constructive_step!(w::World,p::ConstructiveEcology,dt::Real,mechanical_policy)",
  "ecology_step!(w,c,dt;anabolism=allocator)"=>"_rkc_ecology_step!(w,c,dt,mechanical_policy;anabolism=allocator)")
 generated["constructive_cells.jl"]=body
 Base.include_string(@__MODULE__,body,"constructive_cells.jl:rkc-study")
 (;original,generated)
end
function _rkc_ecology_step!(w,c,dt,policy;anabolism=nothing)
 check_dt(w,c,dt)
 _cellular_mechanical_drag(w,c,anabolism)===nothing||throw(ArgumentError("RKC study does not implement matrix drag"))
 issorted(c.ids)||reorder_ecology!(c,sortperm(c.ids))
 _rkc_ecology_step_storage!(w,c,dt,c.parameters,c.weights,c.plastic_marks,anabolism,policy)
end
const CELLULAR_RKC_ECOLOGY_SOURCE=install_cellular_rkc_ecology!()
struct RKCExperimentalEcology{P,M}
 model::P
 policy::M
end
function RKCExperimentalEcology(model;policy=RKCMechanicalPolicy())
 model isa Union{CellularEcology,ConstructiveEcology}||throw(ArgumentError("unsupported RKC study ecology"))
 model isa ConstructiveEcology&&model.mechanics!==nothing&&throw(ArgumentError("matrix drag is outside this study"))
 RKCExperimentalEcology(model,policy)
end
validate(w::World,p::RKCExperimentalEcology)=validate(w,p.model)
stability_limit(w,p::RKCExperimentalEcology)=stability_limit(w,p.model)
process_children(p::RKCExperimentalEcology)=(p.model,)
process_ports(p::RKCExperimentalEcology)=process_ports(p.model)
function advance!(w::World,p::RKCExperimentalEcology,dt::Real)
 check_dt(w,p,dt)
 p.model isa CellularEcology ? _rkc_ecology_step!(w,p.model,dt,p.policy) : _rkc_constructive_step!(w,p.model,dt,p.policy)
 p
end

# The generic Arena inventory/checkpoint contract does not know this unloaded
# study type. Use physical processes for diagnostics; ordinary checkpoints must
# reject it instead of claiming the src-only fingerprint describes its law.
rkc_physical_processes(processes)=map(p->p isa RKCExperimentalEcology ? p.model : p,processes)
arena_process_supported(::RKCExperimentalEcology)=false
