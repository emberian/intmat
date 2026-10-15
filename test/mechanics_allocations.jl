using Test, Random, LinearAlgebra
include("fixtures/bond_network_reference.jl")
function mechanics_bitwise_equal(a,b)
    for name in fieldnames(BondNetwork)
        x=getfield(a,name);y=getfield(b,name)
        if x isa AbstractArray{Float64}
            size(x)==size(y) && reinterpret(UInt64,vec(x))==reinterpret(UInt64,vec(y)) || return false
        elseif x isa Float64
            reinterpret(UInt64,x)==reinterpret(UInt64,y) || return false
        else
            x==y || return false
        end
    end
    true
end
function reference_repair!(p)
    broken=findall(==(1.),p.damage)
    for e in first(broken,min(3,length(broken)))
        i=p.bonds[1,e];j=p.bonds[2,e]
        p.bonds=hcat(p.bonds,[i,j]);push!(p.bond_ids,maximum(p.bond_ids)+1)
        push!(p.rest_length,norm(p.position[:,j]-p.position[:,i]))
        push!(p.stiffness,.7p.stiffness[e]);push!(p.threshold,.02)
        push!(p.damage,0.);push!(p.peak_strain,0.)
    end
end
@testset "Allocation-free bond arithmetic preserves fracture and repair" begin
    w=World(Grid(4,4,4);seed=911)
    candidate=BondNetwork(w;size=(4,4,3),stiffness=5.,failure_strain=.025,stretch_rate=.5,shear_rate=.2)
    candidate.position[1,14]+=.04
    for i in eachindex(candidate.grip)
        candidate.grip[i]==0 && (candidate.load[:,i].=[.002,-.001,.003])
    end
    reference=deepcopy(candidate);initial=length(candidate.bond_ids)
    @test reinterpret(UInt64,vec(bond_forces(candidate)))==reinterpret(UInt64,vec(BondNetworkReference.bond_forces(reference)))
    @test mechanical_energy(candidate)===BondNetworkReference.mechanical_energy(reference)
    @test stability_limit(w,candidate)===BondNetworkReference.stability_limit(w,reference)
    for step in 1:250
        Interior.advance!(w,candidate,.0005);BondNetworkReference.advance!(w,reference,.0005)
        if step in (80,160)
            reference_repair!(candidate);reference_repair!(reference)
        end
        @test mechanics_bitwise_equal(candidate,reference)
    end
    @test count(==(1.),candidate.damage)>0
    @test length(candidate.bond_ids)>initial
    @test candidate.damage_energy>0 && candidate.damping_energy>0
end
