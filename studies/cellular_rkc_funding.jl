# Source-derived funded schedule with the shared production contact workspace.
# The old experimental wrapper remains untouched; no core method is replaced.
function _rkc_max_speed(velocity)
 largest=0.
 for i in axes(velocity,2);largest=max(largest,sqrt(sum(abs2,view(velocity,:,i))));end
 largest
end
function _rkc_affordable(spend,funds)
 funds===nothing&&return true
 for i in eachindex(spend);spend[i]<=funds[i]||return false;end
 true
end
function install_cellular_rkc_funding!()
 path=joinpath(@__DIR__,"cellular_implicit_mechanics.jl");source=read(path,String)
 marker="function cellular_implicit_mechanics!("
 location=findfirst(marker,source);location===nothing&&error("funded schedule marker changed")
 body=source[first(location):end]
 replacements=(
  "function cellular_implicit_mechanics!("=>"function cellular_rkc_funded_mechanics!(",
  "relaxer=cellular_proximal_relax,symmetric_motor=false"=>"relaxer=_CellularRKCCache(),symmetric_motor=true",
  "    index=Dict(id=>i for (i,id) in enumerate(ids))"=>"    bond_range==0||throw(ArgumentError(\"RKC cached funding currently requires zero bond_range\"))\n    index=Dict(id=>i for (i,id) in enumerate(ids))\n    contact_workspace=_CellularContactWorkspace(ids,working_bonds,index)",
  "Val(false),_study_potential_pairs(relaxer,grid,z,radii)).energy"=>"Val(false),_study_potential_pairs(relaxer,grid,z,radii),_cellular_bond_columns(contact_workspace)).energy",
  "refresh=_cellular_refresh_bonds!(grid,ids,x,radii,working_bonds,index,adhesion,adhesion_threshold,break_strain,max_degree,bond_range,stiffness)"=>"refresh=_cellular_refresh_workspace!(contact_workspace,grid,ids,x,radii,working_bonds,index,adhesion,adhesion_threshold,break_strain,max_degree,bond_range,stiffness;pairs=_study_potential_pairs(relaxer,grid,x,radii))",
  "speed=maximum(sqrt(sum(abs2,view(velocity,:,i))) for i in 1:n)"=>"speed=_rkc_max_speed(velocity)",
  "if funds===nothing||all(i->spend[i]<=funds[i],eachindex(spend))"=>"if _rkc_affordable(spend,funds)",
  "repulsion,stiffness,tolerance,max_iterations,max_cg,exact_curvature,"=>"repulsion,stiffness,tolerance,max_iterations,max_cg,exact_curvature,bond_columns=_cellular_bond_columns(contact_workspace),")
 for (old,new) in replacements
  expected=startswith(old,"refresh=") ? 2 : 1
  length(findall(old,body))==expected||error("funded source replacement changed: $old")
  body=replace(body,old=>new)
 end
 Base.include_string(@__MODULE__,body,path*":rkc-workspace")
 (source=source,generated=body)
end
const CELLULAR_RKC_FUNDING_SOURCE=install_cellular_rkc_funding!()
