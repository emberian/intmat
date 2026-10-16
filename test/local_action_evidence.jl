using Test, LinearAlgebra
include("../studies/local_action_evidence.jl")
using .LocalActionEvidence
@testset "source-specific conditional action evidence" begin
    sources=[Evidence(),Evidence()]
    slopes=[.4,-.7];biases=[-2.,3.]
    for i in 1:2, x in (-1.,1.), y in (-1.,1.), noise in (-.3,.3)
        context=[x,y];target=noise*(biases[i]+dot(slopes,context))+10+i+x
        observe!(sources[i],noise,context,target)
    end
    a=fit(sources,1);b=fit(sources,2)
    @test a.slope≈slopes atol=1e-5
    @test b.slope==a.slope
    @test a.intercept≈biases[1] atol=1e-5
    @test b.intercept≈biases[2] atol=1e-5
    @test a.intercept+dot(a.slope,[1.,1.])<0<b.intercept+dot(b.slope,[1.,1.])
    emptychild=Evidence();withchild=fit([sources;[emptychild]],1)
    @test withchild.slope==a.slope
    @test withchild.intercept==a.intercept
    @test withchild.count==a.count
    # Against full block ridge regression: private intercepts, shared slopes.
    A=zeros(4,4);r=zeros(4)
    for i in 1:2
        indices=[i,3,4];A[indices,indices]+=sources[i].gram;r[indices]+=sources[i].rhs
    end
    direct=(A+1e-6I)\r
    @test [a.intercept,b.intercept,a.slope...]≈direct atol=1e-12
    @test fit([Evidence()],1).slope==zeros(2)
    @test_throws ArgumentError observe!(Evidence(),NaN,[0.,0.],0.)
    @test_throws ArgumentError fit(sources,0)
end
