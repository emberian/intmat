using Interior,Test
@testset "Graded-index ray observer" begin
    f=GradedIndex(fill(1.2,12,12,12))
    r=trace_grin(f,(.2,-.1,-.8),(.1,0.,1.);target_z=.8,ds=.013)
    @test r.status==:reached
    @test r.position[1]≈.36 atol=1e-10
    @test r.position[2]≈-.1 atol=1e-10
    @test maximum(abs.(r.direction.-Interior._unit((.1,0.,1.))))<1e-12
    @test r.optical_length≈1.2r.arc_length atol=1e-10
    n=32;a=[1+.1*((i-.5)*2/n-1) for i in 1:n,j in 1:n,k in 1:n]
    g=GradedIndex(a)
    value,grad=index_gradient(g,(.2,.1,0.))
    @test value≈1.02 atol=1e-12
    @test grad[1]≈.1 atol=1e-12
    @test abs(grad[2])+abs(grad[3])<1e-12
    bend=trace_grin(g,(0.,0.,-.7),(0.,0.,1.);target_z=.7,ds=.01)
    @test bend.direction[1]>0
    @test bend.position[1]>0
    # Stratified n(x): longitudinal optical momentum n*v_z is invariant.
    momentum_error(h)=begin
        r=trace_grin(g,(0.,0.,-.7),(0.,0.,1.);target_z=.7,ds=h)
        ne,_=index_gradient(g,r.position)
        abs(ne*r.direction[3]-1.)
    end
    @test momentum_error(.02)<momentum_error(.04)
    @test momentum_error(.01)<2e-6
    ref=trace_grin(g,(0.,0.,-.7),(0.,0.,1.);target_z=.7,ds=.0005).position
    coarse=trace_grin(g,(0.,0.,-.7),(0.,0.,1.);target_z=.7,ds=.04).position
    fine=trace_grin(g,(0.,0.,-.7),(0.,0.,1.);target_z=.7,ds=.02).position
    err(p)=sqrt(sum((p[i]-ref[i])^2 for i in 1:3))
    @test err(fine)<.6err(coarse)
    @test_throws ArgumentError GradedIndex(fill(0.,2,2,2))
    @test_throws ArgumentError GradedIndex(zeros(0,2,2))
    invalid=GradedIndex(ones(2,2,2));invalid.index[1]=NaN
    @test_throws ArgumentError trace_grin(invalid,(0.,0.,0.),(0.,0.,1.))
    @test_throws ArgumentError trace_grin(g,(0.,0.,0.),(0.,0.,0.))
    @test_throws ArgumentError trace_grin(g,(0.,0.,0.),(0.,0.,1.);ds=0)
    limited=trace_grin(g,(0.,0.,-.7),(0.,0.,-1.);max_length=.1)
    @test limited.status==:unresolved
    @test isempty(trace_grin(g,(0.,0.,0.),(0.,0.,1.);save_path=false).path)
    image=render_grin(GradedIndex(ones(8,8,8));width=8,height=8,ds=.03)
    @test image.unresolved==0
    @test all(isfinite,image.rgb)
end
