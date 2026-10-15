using Test
include(joinpath(@__DIR__,"..","studies","lattice_ecology_run.jl"))
cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","alife","lattice-volume32.toml"))
@test_throws ErrorException lattice_run_limits(cfg)
limits=lattice_run_limits(cfg;coordinated_volume=true)
@test limits.frames==42&&limits.native_frame_bytes==99090432
for patch in (Dict("dimensions"=>[64,64,64]),Dict("steps"=>3201),Dict("capture_every"=>40),Dict("seeds"=>[11,29]),Dict("arms"=>["full","no_mutation"]))
 @test_throws ErrorException lattice_run_limits(merge(cfg,patch);coordinated_volume=true)
end
small=TOML.parsefile(joinpath(@__DIR__,"..","experiments","alife","lattice-smoke.toml"))
@test !lattice_run_limits(small).coordinated_volume
