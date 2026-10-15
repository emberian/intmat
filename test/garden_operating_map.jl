using Test, TOML
include(joinpath(@__DIR__,"..","studies","garden_operating_map.jl"))
using .GardenOperatingMap: operating_map
@testset "Native recipe dimensional operating map" begin
    cfg=TOML.parsefile(joinpath(@__DIR__,"..","experiments","native","constructed-garden-64.toml"))
    b=operating_map(cfg)
    @test b.mean_raw_source_rate≈[1/64,1/64]
    @test b.clocks.bare_exchange_diffusion≈[5.4,5.4]
    @test b.clocks.finite_source_duration≈[1500,1500]
    @test b.clocks.plastic_mark_relaxation_range[1]≈1/.022
    @test b.clocks.plastic_mark_relaxation_range[2]≈50
    @test b.supply_matching_enzyme≈fill((1/64)*2.2/(6*2),2)
    # Pure spatial refinement changes no physical characteristic scale.
    finer=deepcopy(cfg);finer["grid"]["shape"].*=2
    @test operating_map(finer)==b
    # A larger domain with proportionately larger sources has the same local
    # operating point and source duration; world volume itself differs.
    larger=deepcopy(cfg);larger["grid"]["extent"].*=2
    for key in ("source_rates","source_budgets");larger["feedstock"][key].*=8;end
    l=operating_map(larger)
    @test l.clocks==b.clocks
    @test l.supply_matching_enzyme==b.supply_matching_enzyme
    @test l.world_volume==8b.world_volume
    slow=operating_map(cfg;catalytic_rate=.6)
    @test slow.supply_matching_total_enzyme≈10b.supply_matching_total_enzyme
    @test slow.permeability_at_supply_match<b.permeability_at_supply_match
    @test slow.clocks.matched_exchange_diffusion[1]>b.clocks.matched_exchange_diffusion[1]
    # Null transport/supply is explicit, never a JSON Infinity or a false zero.
    closed=deepcopy(cfg);closed["resources"]["diffusion"].=0
    closed["feedstock"]["source_rates"].=0
    c=operating_map(closed)
    @test c.clocks.bare_exchange_diffusion==[nothing,nothing]
    @test c.clocks.finite_source_duration==[nothing,nothing]
    @test c.permeability_at_supply_match==1
end
