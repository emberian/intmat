using Test,Random,LinearAlgebra,KernelAbstractions
include(joinpath(@__DIR__,"..","studies","shared_recurrent_runtime.jl"))
const R=SharedRecurrentRuntime
function oracle(theta,x,h0,mask)
    H=size(h0,1);B=size(h0,2);Tsteps=size(x,3);O=size(theta.Wy,1)
    h=zeros(eltype(h0),H,B,Tsteps);y=zeros(eltype(h0),O,B,Tsteps)
    for t in 1:Tsteps,b in 1:B
        previous=mask[b,t]==0 ? zeros(eltype(h0),H) : t==1 ? h0[:,b] : h[:,b,t-1]
        h[:,b,t]=tanh.(theta.Wx*x[:,b,t]+theta.Wh*previous+theta.b)
        y[:,b,t]=theta.Wy*h[:,b,t]+theta.by
    end
    h,y
end
@testset "Shared RNN exact analytic TBPTT and reset boundaries" begin
    rng=MersenneTwister(933);theta=R.initialize(rng,2,3,2;T=Float64)
    x=randn(rng,2,2,4);h0=randn(rng,3,2);dy=randn(rng,2,2,4)
    for mask in (ones(Int32,2,4),Int32[1 0 1 1;0 1 1 1])
        c=R.allocate(theta,2,4);R.forward!(c,theta,x,h0,mask);oh,oy=oracle(theta,x,h0,mask)
        @test c.h≈oh atol=2e-15 rtol=2e-15
        @test c.y≈oy atol=3e-15 rtol=3e-15
        g=deepcopy(R.backward!(c,theta,x,h0,mask,dy));dh0=copy(c.dh0)
        loss()=sum(R.forward!(c,theta,x,h0,mask).*dy)
        epsilon=2e-6
        for key in keys(theta),i in eachindex(theta[key])
            before=theta[key][i];theta[key][i]=before+epsilon;plus=loss()
            theta[key][i]=before-epsilon;minus=loss();theta[key][i]=before
            @test isapprox(g[key][i],(plus-minus)/(2epsilon);atol=3e-8,rtol=3e-7)
        end
        for i in eachindex(h0)
            before=h0[i];h0[i]=before+epsilon;plus=loss();h0[i]=before-epsilon;minus=loss();h0[i]=before
            @test isapprox(dh0[i],(plus-minus)/(2epsilon);atol=3e-8,rtol=3e-7)
        end
    end
    # A reset cuts hidden-state credit, while the shared parameters still train.
    mask=ones(Int32,2,4);mask[:,2].=0;dy[:,:,1:3].=0
    c=R.allocate(theta,2,4);R.forward!(c,theta,x,h0,mask);R.backward!(c,theta,x,h0,mask,dy)
    @test all(iszero,c.dh0)
    @test all(iszero,c.da[:,:,1])
    mask.=0;R.forward!(c,theta,x,h0,mask);R.backward!(c,theta,x,h0,mask,dy)
    @test all(iszero,c.g.Wh)
    @test any(!iszero,c.g.Wx)
    # Repeated one-step use has the same causal recurrent trajectory.
    mask.=1;R.forward!(c,theta,x,h0,mask);previous=copy(h0)
    for t in 1:4
        hn=similar(previous);y=zeros(2,2)
        R.step!(hn,y,theta,x[:,:,t],previous,mask[:,t])
        @test hn==c.h[:,:,t]
        @test y==c.y[:,:,t]
        previous=hn
    end
    @test_throws ArgumentError R.step!(previous,zeros(2,2),theta,x[:,:,1],previous,mask[:,1])
end
@testset "Deterministic chunked gradient and Float32 path" begin
    rng=MersenneTwister(171);theta=R.initialize(rng,3,4,2;T=Float64)
    x=randn(rng,3,129,2);h0=randn(rng,4,129);mask=ones(Int32,129,2);c=R.allocate(theta,129,2)
    dy=ones(2,129,2);R.forward!(c,theta,x,h0,mask);g=deepcopy(R.backward!(c,theta,x,h0,mask,dy))
    @test size(c.partials,2)==3
    @test c.g.by==fill(258.,2)
    R.backward!(c,theta,x,h0,mask,dy);@test isequal(g,c.g)
    f32=map(a->Float32.(a),theta);c32=R.allocate(f32,129,2)
    R.forward!(c32,f32,Float32.(x),Float32.(h0),mask)
    R.backward!(c32,f32,Float32.(x),Float32.(h0),mask,Float32.(dy))
    @test c32.h≈c.h atol=1e-6 rtol=1e-5
    for key in keys(g);@test c32.g[key]≈g[key] atol=1e-4 rtol=2e-5;end
end
@testset "Preallocated shared optimizers and clipping" begin
    theta=R.initialize(MersenneTwister(212),2,3,2;T=Float64);g=map(a->fill(.25,size(a)),theta);P=sum(length,values(theta))
    before=deepcopy(theta);opt=R.optimizer(theta;kind=:sgd,rate=.01,clip=1e6)
    R.update!(theta,g,opt);R.check_status!(opt)
    for key in keys(theta);@test theta[key]≈before[key].-.0025 atol=1e-15;end
    theta=deepcopy(before);opt=R.optimizer(theta;kind=:sgd,rate=.01,clip=.5)
    R.update!(theta,g,opt);R.check_status!(opt)
    for key in keys(theta);@test theta[key]≈before[key].-.005/sqrt(P) atol=1e-15;end
    theta=deepcopy(before);opt=R.optimizer(theta;kind=:adam,rate=.001,clip=1e6)
    for step in 1:2;R.update!(theta,g,opt);R.check_status!(opt);end
    for key in keys(theta);@test theta[key]≈before[key].-2*.001*.25/(.25+1e-8) atol=1e-14;end
    @test opt.step==2
    before=deepcopy(theta);g.Wx[1]=NaN;R.update!(theta,g,opt)
    @test_throws ArgumentError R.check_status!(opt)
    @test isequal(theta,before)
    @test_throws ArgumentError R.optimizer(theta;kind=:unknown)
    @test_throws ArgumentError R.optimizer(theta;rate=Inf)
    # Large finite gradients are clipped safely instead of overflowing sum(g²).
    theta=deepcopy(before);huge=map(a->fill(1e200,size(a)),theta)
    opt=R.optimizer(theta;kind=:sgd,rate=.01,clip=.5)
    R.update!(theta,huge,opt);R.check_status!(opt)
    for key in keys(theta);@test theta[key]≈before[key].-.005/sqrt(P) atol=1e-15;end
end
@testset "Recurrent runtime rejects mismatched workspaces and precision loss" begin
    theta=R.initialize(MersenneTwister(17),2,3,2);c=R.allocate(theta,2,3)
    x=zeros(Float32,2,2,3);h0=zeros(Float32,3,2);mask=ones(Int32,2,3)
    @test_throws ArgumentError R.initialize(MersenneTwister(17),2,3,2;recurrent_gain=Inf)
    @test_throws ArgumentError R.initialize(MersenneTwister(17),2,3,2;recurrent_gain=1e100)
    @test_throws ArgumentError R.optimizer(theta;clip=Inf)
    @test_throws ArgumentError R.optimizer(theta;epsilon=1e100)
    @test_throws ArgumentError R.optimizer(theta;epsilon=1e-100)
    @test_throws ArgumentError R.optimizer(theta;beta1=1-1e-12)
    @test_throws DimensionMismatch R.forward!(merge(c,(partials=zeros(Float32,1,1),)),theta,x,h0,mask)
    @test_throws ArgumentError R.forward!(merge(c,(partials=Float64.(c.partials),)),theta,x,h0,mask)
    @test_throws DimensionMismatch R.forward!(merge(c,(g=merge(c.g,(b=zeros(Float32,2),)),)),theta,x,h0,mask)
    @test_throws ArgumentError R.forward!(c,theta,Float64.(x),h0,mask)
end

@testset "Continuation ABI accepts native Float32 zero/one masks" begin
    theta=R.initialize(MersenneTwister(722),2,3,2);c=R.allocate(theta,2,3)
    x=fill(.1f0,2,2,3);h0=fill(.1f0,3,2);mask=Int32[0 1 1;1 0 1];dy=fill(.3f0,2,2,3)
    R.forward!(c,theta,x,h0,mask);h=copy(c.h);R.backward!(c,theta,x,h0,mask,dy);g=deepcopy(c.g)
    R.forward!(c,theta,x,h0,Float32.(mask));R.backward!(c,theta,x,h0,Float32.(mask),dy)
    @test isequal(h,c.h)
    @test isequal(g,c.g)
end
