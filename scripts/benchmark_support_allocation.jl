#!/usr/bin/env julia
# Transport-only synchronous competition benchmark; no cellular process mutation.
include("benchmark_physical_support.jl")
function allocation_case(side,cells;repeats=3)
    grid=Grid(side,side,side;extent=(4,4,4));rng=MersenneTwister(93)
    positions=4rand(rng,3,cells);base=make_fields(grid)
    reference=CellStencil(grid,positions,.18);candidate=AxisSupport(grid,positions,.18)
    a=deepcopy(base);b=deepcopy(base);workspace=SupportAllocationScratch(candidate)
    requests=.02rand(rng,3,cells);requests[3,:].=0
    reset!(pool)=(foreach(k->copyto!(pool.fields[k],base.fields[k]),1:3);fill!(pool.withdrawn,0))
    ga=allocate_cell_resources!(a,reference,requests)
    gb=allocate_support!(b,candidate,requests,workspace)
    baseline=Float64[];optimized=Float64[];grant_error=0.;field_error=0.;ledger_error=0.
    for pass in 1:repeats
        reset!(a);reset!(b)
        if isodd(pass)
            ta=@timed allocate_cell_resources!(a,reference,requests);ga=ta.value
            tb=@timed allocate_support!(b,candidate,requests,workspace);gb=tb.value
        else
            tb=@timed allocate_support!(b,candidate,requests,workspace);gb=tb.value
            ta=@timed allocate_cell_resources!(a,reference,requests);ga=ta.value
        end
        push!(baseline,ta.time);push!(optimized,tb.time)
        grant_error=max(grant_error,maximum(abs.(ga.-gb)))
        field_error=max(field_error,maximum(maximum(abs.(x.-y)) for (x,y) in zip(a.fields,b.fields)))
        ledger_error=max(ledger_error,abs(sum(gb)-sum(b.withdrawn)))
    end
    grant_error<1e-8 && ledger_error<1e-9 || error("fair allocation mismatch")
    (side=side,cells=cells,extent=4.,radius=.18,baseline_seconds=baseline,candidate_seconds=optimized,
     paired_speedup=median(baseline./optimized),max_grant_error=grant_error,max_field_error=field_error,
     grant_ledger_error=ledger_error,actual_withdrawn=sum(b.withdrawn),
     reference_stencil_logical_bytes=Base.summarysize(reference),candidate_geometry_bytes=support_storage_bytes(candidate),
     candidate_workspace_bytes=sizeof(workspace.largest)+sizeof(workspace.scaled)+sizeof(workspace.actual)+sizeof(workspace.active_sites))
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    out=isempty(ARGS) ? "outputs/physical-support-20261015-v1" : ARGS[1];mkpath(out);rows=Any[]
    allocation_case(8,8;repeats=1)
    for side in (64,128,256)
        r=allocation_case(side,2000);push!(rows,r);println(JSON3.write(r));flush(stdout)
        open(io->JSON3.pretty(io,rows),joinpath(out,"allocation-results.json"),"w");GC.gc()
    end
end
