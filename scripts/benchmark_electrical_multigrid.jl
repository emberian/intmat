#!/usr/bin/env julia
# Same model/operator and convergence tolerance; compare complete cold + warm solves.
using Interior, SHA, Statistics
using Interior: JSON3, KernelAbstractions
include(joinpath(@__DIR__,"..","studies","branch_solver_prototype.jl"))
include(joinpath(@__DIR__,"..","studies","branch_metal_multigrid.jl"))
const MG_BACKEND=isempty(ARGS) ? "cpu" : ARGS[1]
if MG_BACKEND=="amd"
    @eval using AMDGPU
    const mg_upload=AMDGPU.ROCArray
elseif MG_BACKEND=="cpu"
    const mg_upload=identity
else
    error("backend must be cpu or amd")
end
function mg_timed(f,levels)
    backend=KernelAbstractions.get_backend(first(levels).u)
    KernelAbstractions.synchronize(backend)
    started=time_ns(); stats=f(levels)
    KernelAbstractions.synchronize(backend)
    (seconds=(time_ns()-started)*1e-9,cycles=stats.cycles,residual=stats.residual)
end
function mg_gather!(p,levels,upload=mg_upload)
    sites=findall(p.frontier);indices=upload(Int32[LinearIndices(p.potential)[I] for I in sites])
    values=similar(first(levels).u,Float32,length(sites))
    backend=KernelAbstractions.get_backend(values)
    branch_gather_kernel!(backend,256)(values,first(levels).u,indices;ndrange=length(sites))
    for (I,v) in zip(sites,Array(values));p.potential[I]=v;end
end
function mg_case(n;events=3,tol=1.6e-5,workgroup=256)
    mask=specimen(n);p=process(n,mask)
    a=device_levels(hierarchy(copy(p.potential),copy(mask)),mg_upload)
    GC.gc()
    b=device_levels(hierarchy(copy(p.potential),copy(mask)),mg_upload)
    GC.gc()
    baselinegroup=Tuple(parse.(Int,split(get(ENV,"MG_BASE_GROUP","8,8,4"),",")))
    baseline=l->device_solve!(l;tol,maxcycles=500,workgroup=baselinegroup)
    variant="layout"
    rbgroup=Tuple(parse.(Int,split(get(ENV,"MG_RBGROUP","32,4,2"),",")))
    densegroup=Tuple(parse.(Int,split(get(ENV,"MG_DENSEGROUP","32,4,2"),",")))
    rbgroup==densegroup || error("canonical solver uses one shared workgroup shape")
    workspace=get(ENV,"MG_FUSED_NORM","0")=="1" ? DeviceMGNorm(first(b)) : nothing
    candidate=l->device_solve!(l;tol,maxcycles=500,workgroup=rbgroup,norm_workspace=workspace)
    # Warm this exact grid/reduction specialization before measuring cold physics.
    baseline(a); candidate(b)
    copyto!(first(a).u,p.potential); copyto!(first(b).u,p.potential)
    cold_a=mg_timed(baseline,a);cold_b=mg_timed(candidate,b)
    initialdiff=maximum(abs.(Array(first(a).u).-Array(first(b).u)))
    initialdiff<=1e-5 || error("cold potential mismatch $initialdiff")
    warms=Any[]
    for event in 1:events
        mg_gather!(p,a);I=select_site(p,(1.,1.,1.));Interior.electrical_mark!(p,I);p.events+=1
        device_mark!(a,I);device_mark!(b,I)
        # Alternate order to reduce systematic thermal/clock bias.
        if isodd(event);ta=mg_timed(baseline,a);tb=mg_timed(candidate,b)
        else;tb=mg_timed(candidate,b);ta=mg_timed(baseline,a);end
        push!(warms,(event=event,site=Tuple(I),baseline=ta,candidate=tb,speedup=ta.seconds/tb.seconds))
    end
    au=Array(first(a).u);bu=Array(first(b).u)
    diff=maximum(abs.(au.-bu));diff<=1e-5 || error("warm potential mismatch $diff")
    p.potential.=bu;r=Interior.potential_residual(p)
    r<=tol || error("independent residual $r exceeds $tol")
    all(==(1f0),bu[p.occupied]) || error("conductor potential changed")
    all(iszero,bu[1,:,:])&&all(iszero,bu[end,:,:])&&all(iszero,bu[:,1,:])&&all(iszero,bu[:,end,:])&&all(iszero,bu[:,:,1])&&all(iszero,bu[:,:,end]) || error("grounded boundary changed")
    (n=n,events=events,tolerance=tol,workgroup=workgroup,variant=variant,fused_norm=workspace!==nothing,baseline_group=baselinegroup,rbgroup=rbgroup,densegroup=densegroup,cold_baseline=cold_a,cold_candidate=cold_b,
     cold_speedup=cold_a.seconds/cold_b.seconds,cold_potential_difference=initialdiff,
     final_potential_difference=diff,independent_residual=r,warm_samples=warms,
     warm_baseline_median=median(x.baseline.seconds for x in warms),
     warm_candidate_median=median(x.candidate.seconds for x in warms),
     warm_speedup=median(x.baseline.seconds/x.candidate.seconds for x in warms),
     baseline_seconds_per_cycle=median(x.baseline.seconds/x.baseline.cycles for x in warms),
     candidate_seconds_per_cycle=median(x.candidate.seconds/x.candidate.cycles for x in warms))
end
function mg_growth_parity(upload=mg_upload;n=33,events=64,tol=3e-7,fused_norm=false)
    mask=specimen(n);p=process(n,mask);q=process(n,mask)
    p.eta=q.eta=4.;p.seed=q.seed=UInt64(7331)
    a=device_levels(hierarchy(p.potential,p.occupied),upload)
    b=device_levels(hierarchy(q.potential,q.occupied),upload)
    workspace=fused_norm ? DeviceMGNorm(first(b)) : nothing
    device_solve!(a;tol,maxcycles=500);device_solve!(b;tol,maxcycles=500,workgroup=(32,4,2),norm_workspace=workspace)
    trajectory=Tuple{Int,Int,Int}[];started=time_ns()
    for event in 1:events
        mg_gather!(p,a,upload);mg_gather!(q,b,upload)
        I=select_site(p,(1.,1.,1.));J=select_site(q,(1.,1.,1.))
        I==J || error("layout growth diverged at event $event: $I vs $J")
        push!(trajectory,Tuple(I))
        Interior.electrical_mark!(p,I);p.events+=1;device_mark!(a,I)
        Interior.electrical_mark!(q,J);q.events+=1;device_mark!(b,J)
        device_solve!(a;tol,maxcycles=500);device_solve!(b;tol,maxcycles=500,workgroup=(32,4,2),norm_workspace=workspace)
    end
    au=Array(first(a).u);bu=Array(first(b).u);q.potential.=bu
    (n=n,events=events,tolerance=tol,trajectory=trajectory,
     maximum_potential_difference=maximum(abs.(au.-bu)),
     independent_residual=Interior.potential_residual(q),seconds=(time_ns()-started)*1e-9)
end

function main()
    output=length(ARGS)>=2 ? ARGS[2] : "outputs/mg-performance-20261015/results.json"
    ns=length(ARGS)>=3 ? parse.(Int,ARGS[3:end]) : [33,65,129]
    # Compile both complete paths at the same small problem before timing.
    mask=specimen(33)
    for group in ((8,8,4),(32,4,2))
        levels=device_levels(hierarchy(process(33,mask).potential,mask),mg_upload)
        device_solve!(levels;tol=3e-7,maxcycles=500,workgroup=group)
    end
    rows=Any[]
    sourcepaths=["studies/branch_solver_prototype.jl","studies/branch_metal_multigrid.jl","scripts/benchmark_electrical_multigrid.jl"]
    sources=Dict(p=>bytes2hex(sha256(read(joinpath(@__DIR__,"..",p)))) for p in sourcepaths)
    for n in ns
        row=mg_case(n;events=parse(Int,get(ENV,"MG_EVENTS","3")));push!(rows,row);println(JSON3.write(row));flush(stdout)
        mkpath(dirname(output));open(io->JSON3.pretty(io,(backend=MG_BACKEND,julia=string(VERSION),sources=sources,results=rows)),output,"w")
        GC.gc()
    end
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
