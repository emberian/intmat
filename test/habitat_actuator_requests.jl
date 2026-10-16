using Interior, Test, Random
@testset "shared habitat actuator request conversion" begin
    rng=MersenneTwister(61)
    for exponent in (1f0,2f0,3f0), investment in (0f0,.8f0,1f0), allocation in (0f0,.5f0,1f0), k in 1:32
        mass=rand(rng,Float32);dt=rand(rng,Float32);rate=rand(rng,Float32)
        a=(expression=zeros(Float32,1,3),requests=zeros(Float32,1,2));s=(mass=[mass],)
        Interior.dh_set_requests!(a,s,1,(secretion_rate=rate,allocation_exponent=exponent),dt,investment,allocation)
        @test a.expression[1,:]==Float32[investment,allocation,1f0-allocation]
        @test a.requests[1,1]===dt*rate*mass*investment*allocation^exponent
        @test a.requests[1,2]===dt*rate*mass*investment*(1f0-allocation)^exponent
    end
end
