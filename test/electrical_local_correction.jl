using Test
include(joinpath(@__DIR__,"..","studies","electrical_local_correction.jl"))
@testset "local electrical correction certificate" begin
 n=17;p=process(n,specimen(n));levels=hierarchy(p.potential,p.occupied)
 prior=mgsolve!(levels;tol=1e-6);I=select_site(p,(1.,1.,1.))
 before=copy(p.potential)
 Interior.electrical_mark!(p,I);mark_hierarchy!(levels,I)
 electrical_patch!(p.potential,p.occupied,I;radius=2,sweeps=8)
 @test all(K->maximum(abs.(Tuple(K-I)))<=2 || p.potential[K]===before[K],CartesianIndices(before))
 localbound=max(prior.residual,Float64(electrical_patch_residual(p.potential,p.occupied,I,2)))
 actual=residual!(first(levels).residual,p.potential,first(levels).rhs,p.occupied)
 @test actual<=localbound
 @test actual>=prior.residual
 solved=electrical_patch_solve!(levels,I;radius=2,sweeps=8,tol=1e-6,previous_bound=localbound)
 @test solved.residual<=1e-6
 @test Interior.potential_residual(p)<1.1e-6
 @test all(p.potential[p.occupied].==1)
 @test all(p.potential[1,:,:].==0)
 # A missing prior certificate always uses the full residual/fallback.
 first(levels).u[3,3,3]=.9f0
 stats=electrical_patch_solve!(levels,I;radius=1,sweeps=1,tol=1e-6)
 @test stats.cycles>0
 @test Interior.potential_residual(p)<1.1e-6
end
