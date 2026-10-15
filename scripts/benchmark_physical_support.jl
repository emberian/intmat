#!/usr/bin/env julia
# Bounded transport-only CPU measurements; no device or arena implementation edits.
using Interior,Random,Statistics,SHA
using Interior: JSON3
include(joinpath(@__DIR__,"..","studies","physical_support_prototype.jl"))
using .PhysicalSupportPrototype

function support_forecasts(;extent=4.,radius=.18)
    rows=Any[]
    for side in (64,128,256),cells in (2000,20000),capacity_factor in (1,2)
        h=extent/side;a=min(side,2ceil(Int,radius/h)+2);box=a^3
        expected=4pi/3*(radius/h)^3;capacity=capacity_factor*cells
        push!(rows,(side=side,extent=extent,radius=radius,cells=cells,capacity=capacity,
            candidate_box=box,expected_live_edges_per_cell=expected,
            device_dense_support_bytes=8capacity*box,
            device_axis_geometry_bytes=capacity*(24a+12),
            cpu_edge_arrays_bytes_estimate=40cells*expected,
            device_base_grid_bytes=84side^3,
            axis_candidate_work_per_rebuild=cells*box,
            expected_live_edges=cells*expected))
    end
    rows
end
function make_fields(grid)
    pool=ResourcePool(World(grid;fields=());initial=(0,0,0))
    n=grid.size[1];volume=prod(Float64.(grid.extent)./grid.size)
    for species in 1:3,I in CartesianIndices(pool.fields[species])
        x,y,z=Tuple(I)./n
        pool.fields[species][I]=Float32(volume*(.4-.1species+.04sinpi(2x)+.02cospi(2y)+.01sinpi(4z)))
    end
    pool
end
function sample_times(f;repeats=3)
    rows=Float64[];allocated=Int[]
    for _ in 1:repeats
        result=@timed f();push!(rows,result.time);push!(allocated,result.bytes)
    end
    (seconds=median(rows),samples=rows,allocated_bytes=median(allocated))
end
function support_case(side,cells;extent=4.,radius=.18,reference=true,sampling=true)
    grid=Grid(side,side,side;extent=(extent,extent,extent))
    positions=extent.*rand(MersenneTwister(71),3,cells)
    measured=@timed AxisSupport(grid,positions,radius);candidate=measured.value
    build=(seconds=measured.time,allocated_bytes=measured.bytes)
    ref=nothing;refbuild=nothing
    if reference
        r=@timed CellStencil(grid,positions,radius);ref=r.value
        refbuild=(seconds=r.time,allocated_bytes=r.bytes,logical_bytes=Base.summarysize(ref),edges=length(ref.sites))
    end
    results=nothing
    if sampling
        pool=make_fields(grid);concentrations=zeros(3,cells);gradients=zeros(3,3,cells)
        sample_support!(concentrations,gradients,pool.fields,candidate)
        candidate_time=sample_times(()->sample_support!(concentrations,gradients,pool.fields,candidate))
        reference_time=nothing;errors=nothing
        if reference
            expected=sample_concentrations(pool,ref);expected_grad=sample_concentration_gradients(pool,ref)
            errors=(concentration=maximum(abs.(concentrations.-expected)),gradient=maximum(abs.(gradients.-expected_grad)))
            errors.concentration<1e-11 && errors.gradient<1e-10 || error("physical sampling mismatch $errors")
            reference_time=sample_times(()->(sample_concentrations(pool,ref),sample_concentration_gradients(pool,ref)))
        end
        results=(candidate=candidate_time,reference=reference_time,errors=errors,
            sampling_speedup=reference_time===nothing ? nothing : reference_time.seconds/candidate_time.seconds)
    end
    (side=side,extent=extent,radius=radius,cells=cells,prototype_precision="Float64",
     candidate_build=build,reference_build=refbuild,candidate_storage_bytes=support_storage_bytes(candidate),
     candidate_count=sum(Int64,candidate.edge_count),candidate_slots=prod(candidate.counts),sampling=results)
end
function main()
    out=isempty(ARGS) ? "outputs/physical-support-20261015" : ARGS[1];mkpath(out)
    # Precompile exact types using a tiny independent fixture.
    support_case(8,8;radius=.6)
    rows=Any[]
    for side in (64,128,256)
        row=support_case(side,2000);push!(rows,row)
        open(io->JSON3.pretty(io,row),joinpath(out,"side-$side-cells-2000.json"),"w")
        println(JSON3.write(row));flush(stdout);GC.gc()
    end
    # Large-population descriptors and sampling remain bounded without edge arrays.
    for side in (64,128,256)
        row=support_case(side,20000;reference=false);push!(rows,row)
        open(io->JSON3.pretty(io,row),joinpath(out,"side-$side-cells-20000.json"),"w")
        println(JSON3.write(row));flush(stdout);GC.gc()
    end
    sources=Dict(p=>bytes2hex(sha256(read(joinpath(@__DIR__,"..",p)))) for p in
        ("src/cellular_transport.jl","src/device_cellular.jl","src/device_cellular_kernels.jl","studies/physical_support_prototype.jl","scripts/benchmark_physical_support.jl"))
    report=(julia=string(VERSION),threads=Threads.nthreads(),cpu=Sys.CPU_NAME,sources=sources,
        boundary="physical-support construction and concentration+gradient sampling only; no arena/GPU speedup claim",
        measurements=rows,forecasts=vcat(support_forecasts(),support_forecasts(;extent=16.,radius=1.)))
    open(io->JSON3.pretty(io,report),joinpath(out,"results.json"),"w")
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
