# CPU scalar-field epoch comparison, unchanged extent/radius. No GPU launches.
using Interior, Random, Statistics, JSON3, SHA
include(joinpath(@__DIR__,"radial_moment_tiles.jl"))
include(joinpath(@__DIR__,"physical_support_prototype.jl"))
using .RadialMomentTiles, .PhysicalSupportPrototype
function gather_axis(st,field,i)
    result=Ref(0.)
    foreach_support_edge(st,i) do site,weight,x,y,z
        result[]+=weight*Float64(field[site])
    end
    result[]
end
function median_time(f)
    f();median([@elapsed f() for _ in 1:3])
end
function characterize()
    rng=MersenneTwister(20261016);rows=Any[]
    for resolution in (16,128,256),n in (resolution==16 ? (4,) : (512,4096))
        grid=Grid(resolution,resolution,resolution;extent=(4.,4.,4.))
        field=Array{Float32}(undef,grid.size...)
        for I in CartesianIndices(field)
            x,y,z=Tuple(I);field[I]=Float32(1+.2sin(2pi*x/resolution)+.1cos(2pi*y/resolution)+.05sin(4pi*z/resolution))
        end
        positions=rand(rng,3,n).*4;positions[:,1].=[.001,3.999,.001]
        radius=resolution==16 ? 1.5 : .18
        geometry=MomentGeometry(grid,radius;leaf_width=4)
        axis=nothing;axis_build=@elapsed axis=AxisSupport(grid,positions,radius)
        plans=SupportPlan[]
        plan_build=@elapsed for i in 1:n;push!(plans,compile_support(geometry,view(positions,:,i)));end
        lazy=nothing;lazy_build=@elapsed lazy=LazyMomentPyramid(field,geometry,plans)
        dense=nothing;dense_build=@elapsed dense=MomentPyramid(field,grid,radius;leaf_width=4)
        lazy_values=zeros(n);axis_values=zeros(n)
        axis_query=let a=axis,f=field,v=axis_values,count=n;median_time(()->begin for i in 1:count;v[i]=gather_axis(a,f,i);end;end);end
        lazy_query=let p=lazy,ps=plans,v=lazy_values,count=n;median_time(()->begin for i in 1:count;v[i]=sample(p,ps[i]);end;end);end
        # Warm compilation case omitted. Real epoch refresh costs include finite scans
        # of dense fields, and only used-leaf reads for lazy fields.
        lazy_refresh=let p=lazy;median_time(()->refresh!(p));end
        dense_refresh=let p=dense;median_time(()->refresh!(p));end
        dense_query=let p=dense,ps=plans;median_time(()->sum(sample(p,plan) for plan in ps));end
        if resolution!=16
            used=count(!iszero,lazy.index[1]);leaf_voxels=sum(prod(hi.-lo.+1) for I in CartesianIndices(lazy.index[1]) if lazy.index[1][I]!=0 for (lo,hi) in (RadialMomentTiles.bounds(lazy,1,Tuple(I)...),))
            push!(rows,(grid=resolution,cells=n,leaf_width=4,radius=.18,extent=4.,
                axis_build_seconds=axis_build,plan_build_seconds=plan_build,lazy_build_seconds=lazy_build,
                dense_build_seconds=dense_build,lazy_refresh_seconds=lazy_refresh,dense_refresh_seconds=dense_refresh,
                axis_query_seconds=axis_query,lazy_query_seconds=lazy_query,dense_query_seconds=dense_query,
                axis_epoch_seconds=axis_build+axis_query,lazy_epoch_seconds=plan_build+lazy_build+lazy_query,
                lazy_held_geometry_epoch_seconds=lazy_refresh+lazy_query,
                dense_held_geometry_epoch_seconds=dense_refresh+dense_query,
                moment_nodes=size(lazy.moments,2),leaf_tiles_used=used,leaf_voxels_read=leaf_voxels,
                lazy_bytes=storage_bytes(lazy),dense_bytes=storage_bytes(dense),plan_bytes=sum(storage_bytes,plans),
                axis_bytes=support_storage_bytes(axis),maximum_error=maximum(abs.(lazy_values.-axis_values)),
                dense_lazy_exact=all(sample(dense,plans[i])==lazy_values[i] for i in 1:n)))
        end
        println(stderr,"done grid=",resolution," N=",n);flush(stderr)
    end
    (;schema="lazy-radial-cpu-epoch-v1",julia=string(VERSION),threads=Threads.nthreads(),
        notes=["One scalar field only; distinct or mutated fields require separate refreshes.",
        "Held geometry timings exclude plan rebuild; moving cells require new plans and new lazy closure.",
        "Lazy storage uses dense Int32 node maps plus compact14-moment values for used subtrees.",
        "Moment regrouping is not the production bitwise Float32 law; no GPU/coupled speedup claim."],
        study_sha256=bytes2hex(sha256(read(joinpath(@__DIR__,"radial_moment_tiles.jl")))),cases=rows)
end
length(ARGS)==1 || error("one fresh output path required")
ispath(only(ARGS)) && error("output exists")
result=characterize()
open(only(ARGS),"w") do io;JSON3.write(io,result);end
