# Bounded CPU-only characterization. No whole-process/GPU throughput claim.
using Interior, Random, Statistics, JSON3, SHA
include(joinpath(@__DIR__,"radial_moment_tiles.jl"))
include(joinpath(@__DIR__,"physical_support_prototype.jl"))
using .RadialMomentTiles
using .PhysicalSupportPrototype
function axis_query(st,field,i)
    result=Ref(0.)
    foreach_support_edge(st,i) do site,weight,x,y,z
        result[]+=weight*Float64(field[site])
    end
    result[]
end
function timed(f;repeats=3)
    f();median([@elapsed f() for _ in 1:repeats])
end
function characterize()
    # Warm the actual field/moment/plan/query specializations outside samples.
    warmgrid=Grid(8,8,8;extent=(4.,4.,4.));warmfield=ones(Float32,warmgrid.size...)
    p=MomentPyramid(warmfield,warmgrid,.18);plan=compile_support(p,view(reshape([.1,.1,.1],3,1),:,1));sample(p,plan)
    st=AxisSupport(warmgrid,reshape([.1,.1,.1],3,1),.18);axis_query(st,warmfield,1)
    rng=MersenneTwister(20261015);positions=rand(rng,3,24).*4
    positions[:,1].=[.001,3.999,.001]
    rows=Any[]
    for resolution in (64,128,256,512)
        grid=Grid(resolution,resolution,resolution;extent=(4.,4.,4.));field=Array{Float32}(undef,grid.size...)
        for I in CartesianIndices(field)
            x,y,z=Tuple(I);field[I]=Float32(1+.2sin(2pi*x/resolution)+.1cos(2pi*y/resolution)+.05sin(4pi*z/resolution))
        end
        axis=nothing
        axis_build=@elapsed axis=AxisSupport(grid,positions,.18)
        direct_values=[axis_query(axis,field,i) for i in axes(positions,2)]
        axis_query_seconds=let a=axis,f=field,n=size(positions,2);timed(()->sum(axis_query(a,f,i) for i in 1:n));end
        for width in (resolution==512 ? (4,8) : (2,4,8))
            p=nothing;GC.gc()
            build=@elapsed p=MomentPyramid(field,grid,.18;leaf_width=width)
            plans=SupportPlan[]
            plan_build=@elapsed for i in axes(positions,2);push!(plans,compile_support(p,view(positions,:,i)));end
            actual=[sample(p,plan) for plan in plans]
            query=let active_p=p,active_plans=plans;timed(()->sum(sample(active_p,plan) for plan in active_plans));end
            stats=[sample(p,view(positions,:,i)) for i in axes(positions,2)]
            query_allocations=let active_p=p,active_plans=plans;@allocated sum(sample(active_p,plan) for plan in active_plans);end
            push!(rows,(grid=resolution,cells=size(positions,2),radius=.18,extent=4.,leaf_width=width,
                field_bytes=sizeof(field),moment_bytes=storage_bytes(p),plan_bytes=sum(storage_bytes,plans),
                axis_bytes=support_storage_bytes(axis),axis_build_seconds=axis_build,
                moment_build_seconds=build,plan_build_seconds=plan_build,
                axis_query_seconds=axis_query_seconds,moment_plan_query_seconds=query,
                moment_query_allocations=query_allocations,maximum_error=maximum(abs.(actual.-direct_values)),
                accepted_edges=sum(axis.edge_count),direct_boundary_tests=sum(x.direct_voxels for x in stats),
                interior_blocks=sum(x.full_nodes for x in stats),plan_blocks=sum(length(x.blocks) for x in plans),
                complete_single_field_epoch_seconds=build+plan_build+query,
                baseline_single_field_epoch_seconds=axis_build+axis_query_seconds))
        end
    end
    (;schema="radial-moment-cpu-probe-v1",julia=string(VERSION),threads=Threads.nthreads(),
      notes=("Moment/plan/query timings include actual CPU work; no GPU or coupled speedup claim.",
             "Single scalar field epoch; repeated coupled fields require their own refresh or a justified shared transform.",
             "Grid is refined at fixed extent4/radius.18, never shrinking physical support.",
             "Float64 stored moments change reduction order; maximum_error compares existing Float64 axis support."),
      source=(study=bytes2hex(sha256(read(joinpath(@__DIR__,"radial_moment_tiles.jl")))),
              probe=bytes2hex(sha256(read(@__FILE__))),reference=bytes2hex(sha256(read(joinpath(@__DIR__,"physical_support_prototype.jl"))))),
      forecasts=[merge((grid=n,leaf_width=w),memory_forecast((n,n,n);leaf_width=w)) for n in (256,512),w in (2,4,8)],cases=rows)
end
if !isempty(ARGS) && ispath(only(ARGS));error("output exists; choose a fresh file");end
result=characterize()
if isempty(ARGS)
    println(JSON3.write(result))
else
    open(only(ARGS),"w") do io;JSON3.write(io,result);end
end
