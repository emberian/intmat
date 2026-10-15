using Test, Interior, TOML
for file in ("cellular_implicit_mechanics.jl","cellular_pair_proximal.jl","cellular_rkc_mechanics.jl","cellular_rkc_funding.jl","cellular_rkc_ecology.jl")
 Base.include(Interior,joinpath(@__DIR__,"..","studies",file))
end
function rkc_composition_fixture()
 cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","arena","constructive-volume.toml"))
 cfg["grid"]["shape"]=[8,8,8];cfg["grid"]["extent"]=[8.,8.,8.];cfg["dt"]=.01
 cfg["cellular"]["population"]=4;cfg["cellular"]["capacity"]=32;cfg["cellular"]["division_mass"]=.6
 for name in ("resources","feedstock")
  cfg[name]["source_rates"].*=0
  haskey(cfg[name],"sink_rates")&&(cfg[name]["sink_rates"].*=0)
 end
 a=build_arena(cfg);model=only(p for p in a.processes if p isa ConstructiveEcology)
 model.cells.positions.=[2. 2.5 5. 5.5;2. 2. 5. 5.;2. 2. 5. 5.]
 model.cells.biomass[1]=.001;model.cells.metabolites[:,1].=0;model.cells.energy[1]=0
 Interior.rebuild_cell_stencil!(model.cells.stencil,model.cells.positions)
 Interior.cellular_porosity!(a.world,model.cells)
 a
end
@testset "RKC full ecology adapter preserves unchanged mechanics oracle" begin
 a=rkc_composition_fixture();b=deepcopy(a)
 modela=only(p for p in a.processes if p isa ConstructiveEcology)
 modelb=only(p for p in b.processes if p isa ConstructiveEcology)
 b.processes=map(p->p===modelb ? Interior.RKCExperimentalEcology(p;policy=Interior.cellular_mechanics!) : p,b.processes)
 for _ in 1:6
  arena_step!(a);arena_step!(b)
  @test isequal(modela.cells.positions,modelb.cells.positions)
  @test isequal(modela.cells.energy,modelb.cells.energy)
  @test modela.cells.ids==modelb.cells.ids
  @test isequal(modela.cells.hidden,modelb.cells.hidden)
  @test isequal(modela.cells.weights,modelb.cells.weights)
  @test isequal(modela.habitat.enzymes,modelb.habitat.enzymes)
  @test a.world.time===b.world.time
 end
 @test modela.cells.births>0&&modela.cells.deaths>0
 @test modelb.cells.resources===only(p for p in b.processes if p isa ResourcePool&&p===modelb.cells.resources)
end
@testset "RKC biology, construction and closed ledgers compose" begin
 a=rkc_composition_fixture();model=only(p for p in a.processes if p isa ConstructiveEcology)
 policy=Interior.RKCMechanicalPolicy()
 a.processes=map(p->p===model ? Interior.RKCExperimentalEcology(p;policy) : p,a.processes)
 physical=Interior.rkc_physical_processes(a.processes)
 @test only(p for p in physical if p isa ConstructiveEcology)===model
 @test_throws ArgumentError Arena(a.world,a.processes;dt=a.dt)
 mktempdir() do dir
  path=joinpath(dir,"unsupported-study.jls")
  @test_throws ArgumentError checkpoint!(a,path)
  @test !ispath(path)
 end
 mass=Interior.arena_inventory(a.world,physical);energy=Interior.arena_energy(physical)
 for _ in 1:20;arena_step!(a);end
 @test policy.calls==20&&policy.force_evaluations>0
 @test model.cells.births>0&&model.cells.deaths>0
 @test sum(model.habitat.secreted)>0
 @test all(>=(0),model.cells.energy)&&all(>=(0),model.cells.metabolites)
 @test abs(Interior.arena_inventory(a.world,Interior.rkc_physical_processes(a.processes))-mass)<1e-4
 @test abs(Interior.arena_energy(Interior.rkc_physical_processes(a.processes))-energy)<1e-4
 @test a.world.time≈20Float64(Float32(.01))
end
