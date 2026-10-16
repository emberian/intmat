using Test, Interior, Random
@testset "representable proportional metabolic debit" begin
    a=0.006728765367793111;b=0.00251608945166251;stock=a+b
    @test a-stock*a/stock<0 # the original full-oxidation failure
    @test Interior.dc_proportional_debit(a,b,stock)==(0.,0.,a,b)
    rng=MersenneTwister(31017)
    for T in (Float32,Float64)
        # Exercise depletion, near depletion, tiny shares, and wide inventory
        # ratios. Endpoint differences, not nominal requests, must fund output.
        cases=Tuple{T,T,T}[]
        for _ in 1:2000
            a=T(exp2(rand(rng)*80-40));b=T(exp2(rand(rng)*80-40));stock=a+b
            for amount in (zero(T),stock,prevfloat(stock),stock*T(rand(rng)),stock/2)
                push!(cases,(a,b,amount))
            end
        end
        append!(cases,[(T(0),T(0),T(0)),(T(.1),T(.5),T(.6)),
            (T(1),nextfloat(zero(T)),T(.5))])
        for (a,b,amount) in cases
            aa,bb,da,db=Interior.dc_proportional_debit(a,b,amount)
            @test 0<=aa<=a && 0<=bb<=b
            @test da==a-aa && db==b-bb
            @test 0<=da<=a && 0<=db<=b
            @test da+db<=min(a+b,amount)
            @test amount>=a+b || aa+bb>0 || a+b==0
            @test isfinite(aa+bb+da+db)
        end
    end
end
