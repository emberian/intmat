"""Bounded, observation-only recorder for actual per-resident program snapshots.

A mechanism supplies plain completed-state observations. This module never
advances physics, evaluates a policy, reconstructs a prediction, or fills a gap.
"""
module CellProgramObservatory
using Interior: JSON3
export Recorder, capture!, document, write_observatory
mutable struct Recorder
    metadata::Dict{String,Any}
    frames::Vector{Any}
    selected_ids::Vector{Int}
    max_frames::Int
    max_cells::Int
    max_bytes::Int
    bytes::Int
    stop_reason::Union{Nothing,String}
end
plain(x)=JSON3.read(JSON3.write(x),Dict{String,Any})
function Recorder(metadata;ids=Int[],max_frames=256,max_cells=64,max_bytes=8_000_000)
    all(x->x isa Integer && !(x isa Bool) && x>0,(max_frames,max_cells,max_bytes)) || throw(ArgumentError("positive integer recorder bounds required"))
    all(x->x isa Integer && !(x isa Bool) && x>0,ids) || throw(ArgumentError("positive resident IDs required"))
    length(unique(ids))==length(ids) && length(ids)<=max_cells || throw(ArgumentError("unique selected IDs within cell bound required"))
    m=plain(metadata);bytes=ncodeunits(JSON3.write(m))+2048
    bytes<max_bytes || throw(ArgumentError("metadata exceeds recorder budget"))
    Recorder(m,Any[],Int.(ids),max_frames,max_cells,max_bytes,bytes,nothing)
end
finite_number(x)=x isa Real && !(x isa Bool) && isfinite(x)
function check_frame(frame)
    finite_number(get(frame,"time",nothing)) && frame["time"]>=0 || throw(ArgumentError("finite nonnegative observation time required"))
    step=get(frame,"step",nothing);step isa Integer && !(step isa Bool) && step>=0 || throw(ArgumentError("nonnegative actual step required"))
    cells=get(frame,"cells",nothing);cells isa AbstractVector || throw(ArgumentError("actual cells array required"))
    ids=Int[]
    for cell in cells
        cell isa AbstractDict || throw(ArgumentError("plain cell records required"))
        id=get(cell,"id",nothing);id isa Integer && !(id isa Bool) && id>0 || throw(ArgumentError("positive resident ID required"))
        push!(ids,id)
        for key in ("raw","food","catalyst","fast_raw","slow_raw","raw_next","prediction_error","funded","requested")
            haskey(cell,key) || continue
            cell[key]===nothing && continue
            cell[key] isa AbstractVector && length(cell[key])==2 && all(finite_number,cell[key]) || throw(ArgumentError("$key must be two actual finite values or null"))
        end
        if haskey(cell,"position")
            cell["position"] isa AbstractVector && length(cell["position"])==3 && all(finite_number,cell["position"]) || throw(ArgumentError("position must contain three actual finite coordinates"))
        end
    end
    length(unique(ids))==length(ids) || throw(ArgumentError("duplicate resident ID in one capture"))
    nothing
end
function capture!(r::Recorder,observation)
    r.stop_reason===nothing || return false
    frame=plain(observation);check_frame(frame)
    if !isempty(r.frames)
        previous=last(r.frames)
        frame["time"]>=previous["time"] && frame["step"]>=previous["step"] || throw(ArgumentError("completed capture clocks must not move backward"))
    end
    if !isempty(r.selected_ids)
        frame["cells"]=[c for c in frame["cells"] if c["id"] in r.selected_ids]
    end
    reason=length(r.frames)>=r.max_frames ? "frame cap reached" : length(frame["cells"])>r.max_cells ? "cell cap reached; select IDs or use a larger declared bound" : nothing
    bytes=ncodeunits(JSON3.write(frame))+1
    reason===nothing && r.bytes+bytes>r.max_bytes && (reason="byte cap reached")
    if reason!==nothing;r.stop_reason=reason;return false;end
    push!(r.frames,frame);r.bytes+=bytes;true
end
function document(r::Recorder)
    (schema="interior-cell-program-observatory-1",metadata=r.metadata,
        recorder=(selected_ids=r.selected_ids,max_frames=r.max_frames,max_cells=r.max_cells,max_bytes=r.max_bytes,
            captures=length(r.frames),stop_reason=r.stop_reason,note="Actual completed captures only; missing fields and residents remain unavailable."),frames=r.frames)
end
function write_observatory(r::Recorder,path)
    ispath(path) && throw(ArgumentError("fresh observatory file required"))
    payload=JSON3.write(document(r));ncodeunits(payload)+1<=r.max_bytes || throw(ArgumentError("encoded document exceeds recorder budget"))
    open(path*".pending","w") do io;write(io,payload);write(io,'\n');end
    mv(path*".pending",path);path
end
end
