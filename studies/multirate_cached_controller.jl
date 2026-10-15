# Source-derived operators: retain the native neural update, separate its
# support/request work, and hold only the eight chemical sensor channels.
using Interior: KernelAbstractions
using Interior.KernelAbstractions: @kernel, @index
using Interior: CELLULAR_SENSORS,CELLULAR_ACTIONS,dc_squash,dc_bond_valid,dc_distance,
    dc_bin,dc_search_range,dc_search_count,dc_neighbor_bin,dc_index,
    device_support_edges,device_support_descriptor
const CACHED_CONTROLLER_SOURCE=joinpath(pkgdir(Interior),"src","device_cellular_kernels.jl")
const CACHED_CONTROLLER_SOURCE_SHA256=bytes2hex(Interior.SHA.sha256(read(CACHED_CONTROLLER_SOURCE)))
function native_section(name,nextname)
    source=read(CACHED_CONTROLLER_SOURCE,String)
    a=findfirst("@kernel function $name",source)
    a===nothing&&error("native operator start changed")
    b=findnext("@kernel function $nextname",source,last(a)+1)
    b===nothing&&error("native operator end changed")
    source[first(a):prevind(source,first(b))]
end
function derived_controller_operators!()
    sense=native_section("dc_sense!","dc_controller!")
    initial="for j in 1:CELLULAR_SENSORS;s.inputs[i,base+j]=0;end"
    occursin(initial,sense)||error("native sense reset changed")
    sense=replace(sense,initial=>"for j in 9:CELLULAR_SENSORS;s.inputs[i,base+j]=0;end")
    a=findfirst("        concentration1=",sense);b=findfirst("        degree=0",sense)
    a===nothing||b===nothing ? error("native chemical sense markers changed") : nothing
    sense=sense[1:prevind(sense,first(a))]*sense[first(b):end]
    sense=replace(sense,"dc_sense!"=>"mr_cached_sense!")
    regulate=native_section("dc_regulate!","dc_allocate_sites!")
    room=findfirst("        room=",regulate);neural=findfirst("        alpha=",regulate)
    scatter=findfirst("\n        for (site,weight)",regulate)
    any(isnothing,(room,neural,scatter))&&error("native regulation markers changed")
    fast=regulate[1:prevind(regulate,first(room))]*regulate[first(neural):prevind(regulate,first(scatter))]*"\n    end\nend\n"
    request=regulate[1:prevind(regulate,first(neural))]*regulate[first(scatter):end]
    fast=replace(fast,"dc_regulate!"=>"mr_fast_neural!")
    request=replace(request,"dc_regulate!"=>"mr_request_only!")
    for text in (sense,fast,request);Core.eval(@__MODULE__,Meta.parse(text));end
end
derived_controller_operators!()
