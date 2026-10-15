using Interior, Test, Random, KernelAbstractions
isdefined(Interior,:DeviceCellular) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_cellular.jl"))
isdefined(Interior,:DeviceCatalyticHabitat) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_habitat.jl"))
isdefined(Interior,:dh_sorption_target) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_habitat_sorption.jl"))

@testset "Device finite-capacity sorption target" begin
    cases=NTuple{7,Float32}[
        (.2,1,2,2,.1,1,.3), (.2,1,1,2,0,1,.3), # ordinary / repeated root
        (.2,1,2,0,.1,1,.3), (.2,1,2,0,0,1,.3), # no adsorption / frozen
        (.2,1,2,2,0,1,.3), (0,1,0,2,.1,1,.3), # no desorption / no capacity
        (.2,1,2,2,.1,1,0), (0,0,2,2,.1,1,.3),
        (prevfloat(1f0),1,1,2,0,1,1), (.99,1,1,1,1f-12,1,100),
        (1f-12,1f-10,1f-8,1f-4,1f-6,.01,10)]
    equilibrium=Float32(Interior.sorption_target(0.,1.,2.,2.,.1,1.,1e6))
    push!(cases,(prevfloat(equilibrium),1,2,2,.1,1,.3),
        (nextfloat(equilibrium),1,2,2,.1,1,.3))
    rng=MersenneTwister(1927)
    for _ in 1:1000
        m,c=Float32.(10. .^ (4rand(rng,2).-2));b=rand(rng,Float32)*min(m,c)
        on,off,v,t=Float32.(10. .^ (4rand(rng,4).-2))
        push!(cases,(b,m,c,on,off,v,t))
    end
    for args in cases
        got=Interior.dh_sorption_target(args...)
        oracle=Float32(Interior.sorption_target(map(Float64,args)...))
        @test got isa Float32
        @test got==oracle
        @test isfinite(got) && 0<=got<=min(args[2],args[3])
    end
    for (b,m,c,on,off,v,t) in cases[1:13]
        whole=Interior.dh_sorption_target(b,m,c,on,off,v,t)
        half=Interior.dh_sorption_target(b,m,c,on,off,v,t/2)
        twice=Interior.dh_sorption_target(half,m,c,on,off,v,t/2)
        @test isapprox(whole,twice;rtol=3eps(Float32),atol=3eps(Float32)*min(m,c))
    end
    for _ in 1:400
        m,c=Float32.(10. .^ (2rand(rng,2).-1));b=rand(rng,Float32)*min(m,c)
        on,off,v,t=Float32.(10. .^ (2rand(rng,4).-1))
        whole=Interior.dh_sorption_target(b,m,c,on,off,v,t)
        half=Interior.dh_sorption_target(b,m,c,on,off,v,t/2)
        twice=Interior.dh_sorption_target(half,m,c,on,off,v,t/2)
        @test isapprox(whole,twice;rtol=3eps(Float32),atol=3eps(Float32)*min(m,c))
    end
end

@testset "Representable sorption transfer and compulsory collapse" begin
    rng=MersenneTwister(518)
    cases=Tuple{Float32,Float32,Float32}[(1,0,.5),(0,1,.5),(1,1,1),
        (1f20,1,0),(1,1f20,0),(0,1f-30,0),(floatmax(Float32),1,0),
        (1f-20,1f20,nextfloat(1f20)),(1f20,1f-20,nextfloat(1f-20))]
    for _ in 1:1000
        f,b,t=Float32.(10. .^ (24rand(rng,3).-12));push!(cases,(f,b,t))
    end
    for (free,bound,target) in cases
        result=Interior.dh_sorption_exchange(free,bound,target)
        @test result.free isa Float32 && result.bound isa Float32 && result.inert isa Float64
        @test isfinite(result.free) && isfinite(result.bound) && isfinite(result.inert)
        @test min(result.free,result.bound,result.inert)>=0
        before=Float64(free)+Float64(bound)
        after=Float64(result.free)+Float64(result.bound)+result.inert
        @test isapprox(before,after;rtol=4eps(Float64),atol=0)
        if target>bound
            @test bound<=result.bound<=target
            @test result.free<=free
        else
            @test target<=result.bound<=bound
            @test result.free>=free
        end
        capacity=min(bound,target)
        collapsed=Interior.dh_sorption_collapse(free,bound,capacity)
        @test 0<=collapsed.bound<=capacity
        @test min(collapsed.free,collapsed.inert)>=0
        @test isapprox(before,Float64(collapsed.free)+Float64(collapsed.bound)+collapsed.inert;rtol=4eps(Float64),atol=0)
        # The same potential for both endpoints makes inert potential the heat.
        potential=3.25
        @test isapprox(potential*before,potential*(Float64(collapsed.free)+Float64(collapsed.bound))+potential*collapsed.inert;rtol=5eps(Float64),atol=0)
    end
    # Preserve sub-ulp inert mass even when no additional compulsory debit is
    # needed after ordinary desorption (large bound + tiny free endpoint).
    free,bound=1f-12,nextfloat(1f0)
    precise=Interior.dh_sorption_collapse(free,bound,1f0)
    @test BigFloat(free)+BigFloat(bound)==BigFloat(precise.free)+BigFloat(precise.bound)+BigFloat(precise.inert)
    @test precise.inert>0
    collapsed=Interior.dh_sorption_collapse(1f20,1f0,0f0)
    @test collapsed.free==1f20 && collapsed.bound==0 && collapsed.inert==1
end

# CPU KA compilation establishes that the pure scalar routines compose inside
# a site kernel; no host loop performs any chemistry in this execution probe.
@kernel function device_sorption_scalar_probe!(out,free,bound,capacity)
    i=@index(Global,Linear)
    collapsed=Interior.dh_sorption_collapse(free[i],bound[i],capacity[i])
    target=Interior.dh_sorption_target(collapsed.bound,collapsed.free+collapsed.bound,capacity[i],2f0,.1f0,1f0,.25f0)
    value=Interior.dh_sorption_exchange(collapsed.free,collapsed.bound,target)
    out[i,1]=value.free;out[i,2]=value.bound;out[i,3]=collapsed.inert+value.inert
end
@testset "Sorption scalar composition in CPU KA kernel" begin
    free=Float32[1,.1,1f20];bound=Float32[.2,.3,1];capacity=Float32[2,.1,0]
    out=zeros(Float64,3,3)
    device_sorption_scalar_probe!(CPU(),32)(out,free,bound,capacity;ndrange=3)
    KernelAbstractions.synchronize(CPU())
    @test all(isfinite,out) && all(>=(0),out)
    @test all(out[:,2].<=capacity)
    for i in eachindex(free)
        @test isapprox(Float64(free[i])+Float64(bound[i]),sum(out[i,:]);rtol=4eps(Float64),atol=0)
    end
end

@testset "Finite Float32 pools can have a Float64 total" begin
    largest=floatmax(Float32)
    target=Interior.dh_sorption_target(largest/2,Float64(largest)+Float64(largest/2),largest,1f0,1f0,1f0,.1f0)
    @test isfinite(target)
    @test 0<=target<=largest
end
