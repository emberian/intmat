using Test, Interior, Random, LinearAlgebra
include(joinpath(@__DIR__,"..","studies","pair_transport_ecology.jl"))
using .PairTransportEcology
@testset "Edge exponential transport study" begin
    r=transport_fixture(size=8,founders=2,capacity=16,diffusion=1.)
    p=r.process.cells;g=p.geometry;q=p.parameters
    rng=MersenneTwister(91);u=rand(rng,Float32,512,3);perm=rand(rng,Float32,512)
    e=merge(p.environment,(fields=copy(u),scratch=similar(u),permeability=perm))
    before=sum(Float64,e.fields);pair_transport!(p,e,q,g,10.)
    @test all(isfinite,e.fields)&&minimum(e.fields)>=0
    @test abs(sum(Float64,e.fields)-before)<1e-4
    @test e.fields==e.scratch
    frozen=merge(e,(fields=copy(u),permeability=zeros(Float32,512)))
    pair_transport!(p,frozen,q,g,10.)
    @test frozen.fields==u
    @test_throws ErrorException pair_transport!(p,e,q,merge(g,(dims=(Int32(7),Int32(8),Int32(8)),)),.1)
    # Fourier mode compares against exact semi-discrete heat flow, not continuum.
    initial=repeat(Float32[1+.4cos(2pi*(i-1)/8) for i in 1:8],64)
    duration=.05;lambda=-4Float64(q.diffusion[1])*sin(pi/8)^2/Float64(g.spacing[1])^2
    exact=1 .+(Float64.(initial).-1).*exp(lambda*duration)
    errors=Float64[]
    for cycles in (1,2,4,8)
        ee=merge(e,(fields=repeat(initial,1,3),scratch=similar(u),permeability=ones(Float32,512)))
        pair_transport!(p,ee,q,g,duration;cycles)
        push!(errors,sqrt(sum(abs2,Float64.(ee.fields[:,1]).-exact)/512))
    end
    @test all(diff(errors).<0)
    @test errors[1]/errors[2]>3
    @test errors[2]/errors[3]>3
    @test p.environment.status[1]==0
    # Complete cells/uptake/junction/metabolism/construction/chemistry/lifecycle.
    for _ in 1:4;transport_step!(r,.002);end
    audit=Interior.device_habitat_audit(Interior.device_habitat_snapshot(r.process))
    @test audit.finite
    @test audit.population==2
end
@testset "Cyclic implicit line transport" begin
    rng=MersenneTwister(92)
    for n in (2,3,7,8),dt in (.001,1.,100.)
        f=rand(rng,Float32,n,3);perm=rand(rng,Float32,n);perm[1]=0
        e=(fields=copy(f),scratch=similar(f),permeability=perm,status=zeros(Int32,1))
        g=(dims=(Int32(n),Int32(1),Int32(1)),spacing=(.25f0,1f0,1f0));q=(diffusion=(.2f0,.3f0,0f0),)
        expected=copy(f)
        for k in 1:3
            A=Matrix{Float64}(I,n,n)
            for i in 1:n,direction in (-1,1)
                j=mod1(i+direction,n);rate=dt*q.diffusion[k]*(Float64(perm[i])+perm[j])*.5/.25^2
                A[i,i]+=rate;A[i,j]-=rate
            end
            expected[:,k]=A\Float64.(f[:,k])
        end
        PairTransportEcology.line_transport!(e,q,g,dt)
        @test maximum(abs.(e.fields.-expected))<2e-6
        @test minimum(e.fields)>=0
        @test abs(sum(Float64,e.fields)-sum(Float64,f))<2e-6
    end
    @test PairTransportEcology.nondiffusive_limit((response_time=100f0,adaptation_time=100f0))==.25
    @test PairTransportEcology.nondiffusive_limit((response_time=.1f0,adaptation_time=1f0))<=.05
    # Exact lattice Fourier decay separates line time error from grid error.
    n=8;f=Float32[1+.4cos(2pi*(i-1)/n) for i in 1:n]
    g=(dims=(Int32(n),Int32(1),Int32(1)),spacing=(.25f0,1f0,1f0));q=(diffusion=(1f0,1f0,1f0),)
    duration=.02;lambda=-4sin(pi/n)^2/.25^2;exact=1 .+(Float64.(f).-1).*exp(lambda*duration)
    errors=Float64[]
    for cycles in (1,2,4,8)
        e=(fields=repeat(f,1,3),scratch=zeros(Float32,n,3),permeability=ones(Float32,n),status=zeros(Int32,1))
        PairTransportEcology.line_transport!(e,q,g,duration;cycles)
        push!(errors,maximum(abs.(Float64.(e.fields[:,1]).-exact)))
    end
    @test all(diff(errors).<0)
    @test errors[2]/errors[3]>1.8
    @test errors[3]/errors[4]>1.8
end
