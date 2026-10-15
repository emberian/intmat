using Test, Interior
include(joinpath(@__DIR__,"..","studies","branch_solver_prototype.jl"))
include(joinpath(@__DIR__,"..","studies","branch_metal_multigrid.jl"))
@testset "Device multigrid workgroup layout" begin
    n=33;mask=specimen(n)
    a=device_levels(hierarchy(process(n,mask).potential,mask),identity)
    b=device_levels(hierarchy(process(n,mask).potential,mask),identity)
    device_vcycle!(a);device_vcycle!(b;workgroup=(32,4,2))
    @test first(a).u==first(b).u
    # A complete converged solve, not only matching a smoother launch.
    sa=device_solve!(a;tol=3e-7,maxcycles=500)
    sb=device_solve!(b;workgroup=(32,4,2),tol=3e-7,maxcycles=500)
    @test sa==sb
    @test first(a).u==first(b).u
    @test Interior.potential_residual(process_with_u(n,mask,first(b).u))<=3e-7
    @test all(==(1f0),first(b).u[mask])
    @test all(iszero,first(b).u[1,:,:])
    @test all(iszero,first(b).u[end,:,:])
    # Independent SOR law and exact Dirichlet geometry on a changed conductor.
    p=process(n,mask);p.tolerance=3e-7
    Interior.solve_potential!(p)
    @test maximum(abs.(p.potential.-first(b).u))<1e-5
    c=(n+1)÷2;I=CartesianIndex(c+1,c,c)
    for levels in (a,b);device_mark!(levels,I);end
    sa=device_solve!(a;tol=3e-7,maxcycles=500)
    sb=device_solve!(b;workgroup=(32,4,2),tol=3e-7,maxcycles=500)
    @test sa==sb && first(a).u==first(b).u
    # Coarse and fine hierarchies retain the same solve law.
    shortmask=specimen(17)
    x=device_levels(hierarchy(process(17,shortmask).potential,shortmask),identity)
    y=device_levels(hierarchy(process(17,shortmask).potential,shortmask),identity)
    @test device_solve!(x)==device_solve!(y;workgroup=(32,4,2))
    @test x[1].u==y[1].u
end
@testset "Fused convergence residual norm" begin
    # Non-multiple final workgroups and unrelaxed fields exercise the reduction,
    # rather than only checking an already tiny converged residual.
    for n in (5,17,33,65)
        mask=specimen(n);p=process(n,mask)
        for I in CartesianIndices(p.potential)
            p.potential[I]=Float32(sin(I[1]+.3I[2]-.2I[3]))
        end
        levels=device_levels(hierarchy(p.potential,mask),identity)
        workspace=DeviceMGNorm(first(levels))
        device_residual!(first(levels))
        reference=Float64(maximum(abs,first(levels).residual))/6
        @test device_residual_norm!(workspace,first(levels))==reference
    end
    mask=specimen(33)
    a=device_levels(hierarchy(process(33,mask).potential,mask),identity)
    b=device_levels(hierarchy(process(33,mask).potential,mask),identity)
    workspace=DeviceMGNorm(first(b))
    @test device_solve!(a;tol=3e-7)==device_solve!(b;tol=3e-7,workgroup=(32,4,2),norm_workspace=workspace)
    @test first(a).u==first(b).u
    short=device_levels(hierarchy(process(17,specimen(17)).potential,specimen(17)),identity)
    @test_throws DimensionMismatch device_residual_norm!(workspace,first(short))
end
