using Interior, Test

@testset "Catalytic deposition rejects invalid requests and terminates at tiny scales" begin
    stencil=(offsets=[1,4],sites=[1,2,3],weights=[.2,.3,.5])
    for invalid in (-4.6543e-319,-nextfloat(0.),NaN,Inf,-Inf)
        field=zeros(3);before=copy(field)
        @test_throws ArgumentError Interior._catalytic_deposit!(field,stencil,1,invalid)
        @test field==before
    end
    for amount in (0.,nextfloat(0.),3nextfloat(0.),1e-319,floatmin(Float64),1e-100,.1),
            initial in (0.,1e-320,1.,1e100)
        field=fill(initial,3);before=copy(field)
        actual=Interior._catalytic_deposit!(field,stencil,1,amount)
        @test isfinite(actual) && 0<=actual<=amount
        @test all(>=(initial),field)
        # Sum endpoint changes at higher precision to expose subnormal losses;
        # only aggregate Float64 summation roundoff is allowed.
        received=sum(BigFloat(field[i])-BigFloat(before[i]) for i in eachindex(field))
        @test abs(received-BigFloat(actual))<=4BigFloat(eps(max(actual,nextfloat(0.))))
    end
end

@testset "Complete respiration cannot create a negative complementary donor" begin
    w=World(Grid(1,1;extent=(1.,1.,1.));seed=44);fill!(w.fields.nutrient,0)
    r=ResourcePool(w;potentials=(4,4,0),initial=(0,0,0),diffusion=(0,0,0))
    c=CellularEcology(w,reshape([.5,.5,0.],3,1),r;capacity=2,hidden=2,
        energy=0.,metabolites=0.,uptake_rate=0.,basal_rate=0.,neural_cost=0.,
        catabolic_rate=1e6,growth_rate=0.,speed=0.,repulsion=0.,stiffness=0.,
        division_mass=10.,mutation=0.,developmental_partition=0.)
    # The sum rounds to the larger donor, which is exactly Float32. The old
    # multiply-then-divide withdrawal overdraws the smaller donor by one ulp.
    tiny=4.808073959212107e-36;large=3.92216432942849e-11
    @test large==Float64(Float32(large)) && tiny+large==large
    @test (tiny+large)*tiny/(tiny+large)>tiny
    c.metabolites[:,1].=(tiny,large)
    f=ResourcePool(w;names=(:raw_a,:raw_b),potentials=(6,6),initial=(0,0),diffusion=(0,0))
    h=CatalyticHabitat(w,c,f;regulation=false,secretion_rate=1.,catalytic_rate=0.,decay_rate=0.,mutation=0.)
    p=ConstructiveEcology(w,c,h)
    mass=Interior.arena_inventory(w,(p,));potential=Interior.arena_energy((p,))
    step!(w,p,.01)
    @test all(x->isfinite(x)&&x>=0,c.metabolites) && all(>=(0),c.energy)
    @test w.time==Float64(.01f0)
    @test Interior.arena_inventory(w,(p,))≈mass atol=1e-15 rtol=0
    @test Interior.arena_energy((p,))≈potential atol=128eps(potential) rtol=0
    @test all(x->isfinite(x)&&x>=0,h.secreted)
    for (a,b,amount) in ((nextfloat(0.),nextfloat(0.),2nextfloat(0.)),
            (3nextfloat(0.),5nextfloat(0.),5nextfloat(0.)),
            (floatmin(Float64),floatmin(Float64),floatmin(Float64)),
            (tiny,large,tiny+large),(.1,.2,.15))
        c.metabolites[:,1].=(a,b);c.inert_residual=0.;fill!(r.fields[3],0)
        before=copy(c.metabolites[:,1])
        degraded=Interior._cellular_degrade_precursors!(c,1,amount)
        removed=sum(BigFloat(before[k])-BigFloat(c.metabolites[k,1]) for k in 1:2)
        received=sum(BigFloat,r.fields[3])+BigFloat(c.inert_residual)
        @test all(>=(0),c.metabolites) && degraded.mass>=0 && c.inert_residual>=0
        if max(a,b)<=floatmin(Float64)
            @test removed==received
        else
            @test abs(removed-received)<=4BigFloat(eps(Float64(removed)))
        end
        endpoint_potential=sum((BigFloat(before[k])-BigFloat(c.metabolites[k,1]))*BigFloat(r.potentials[k]) for k in 1:2)
        @test abs(endpoint_potential-BigFloat(degraded.potential))<=4BigFloat(eps(max(degraded.potential,nextfloat(0.))))
    end
end
