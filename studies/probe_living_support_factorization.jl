using Test,Interior
include("living_support_worklists.jl")
include("living_support_factorization.jl")
const B=LivingSupportWorklists
const F=LivingSupportFactorization
@testset "ordered support factorization" begin
    for n in (4,16,64),founders in (24,46)
        a=B.fixture(n,founders);b=deepcopy(a);owner=F.Carrier(b.p.cells)
        Interior.device_cellular_launch!(a.p.cells,Interior.dc_support!,a.p.cells.state,a.p.cells.environment,a.p.cells.parameters,a.p.cells.geometry;ndrange=512)
        Interior.device_cellular_launch!(owner,Interior.dc_support!,b.p.cells.state,b.p.cells.environment,b.p.cells.parameters,b.p.cells.geometry;ndrange=512)
        @test B.exact_compare(a,b)>0
    end
end
