using Test
include(joinpath(@__DIR__,"..","studies","electrical_temporal_solver.jl"))
@testset "tracked harmonic certificate" begin
 p=process(17,specimen(17));levels=device_levels(hierarchy(p.potential,p.occupied),identity)
 solver=ElectricalTemporalSolver(levels)
 initial=electrical_temporal_solve!(solver)
 @test !initial.patch_used && initial.cycles>0
 @test electrical_temporal_solve!(solver).cached
 I=select_site(p,(1.,1.,1.));Interior.electrical_mark!(p,I);electrical_temporal_mark!(solver,I)
 @test !solver.certificate&&solver.patch_pending
 stats=electrical_temporal_solve!(solver)
 @test stats.patch_used
 p.potential.=first(levels).u
 @test Interior.potential_residual(p)<=2e-5
 @test solver.certificate&&solver.certified_generation==solver.generation
 # Unknown field edits are never silently reused as a valid certificate.
 first(levels).u[3,3,3]=.9f0;electrical_temporal_invalidate!(solver,"unknown fine field edit")
 @test !solver.patch_enabled&&!solver.certificate
 @test !electrical_temporal_solve!(solver).patch_used
 @test !solver.patch_enabled
 @test_throws ArgumentError electrical_temporal_recertify!(solver)
 @test !electrical_temporal_recertify!(solver;homogeneous_operator=true).patch_used
 @test solver.patch_enabled
 # Multiple additions invalidate the one-patch support argument.
 for _ in 1:2
  I=select_site(p,(1.,1.,1.));Interior.electrical_mark!(p,I);p.events+=1
  electrical_temporal_mark!(solver,I)
 end
 @test !solver.patch_pending&&solver.pending_count==2
 @test !electrical_temporal_solve!(solver).patch_used
 @test solver.certificate
 # Even rejected recertification must poison the old certificate.
 first(levels).rhs[3,3,3]=1f0
 @test_throws ArgumentError electrical_temporal_recertify!(solver;homogeneous_operator=true)
 @test !solver.patch_enabled&&!solver.certificate
 @test_throws ArgumentError electrical_temporal_mark!(solver,CartesianIndex(1,2,2))
end
