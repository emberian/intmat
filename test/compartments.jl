@testset "active compartments" begin
    function medium_fixture(;barrier=0f0)
        w=World(Grid(16,8;extent=(16,8,1)))
        p=ImplicitPhaseSeparation(w;mean=0,noise=0,mineral_pinning=0,
            nutrient_diffusion=0,interface_activity=0,activity_relaxation=0)
        m=ActiveMedium(w,p;interface_barrier=0,membrane_barrier=barrier,
            signal_diffusion=.1,nutrient_diffusion=.2)
        w,p,m
    end
    w,p,m=medium_fixture()
    fill!(w.fields.nutrient,0); w.fields.nutrient[1:8,:,:].=1
    before=sum(w.fields.nutrient)
    advance!(w,m,20)
    @test isapprox(sum(w.fields.nutrient),before;rtol=1e-6)
    @test minimum(w.fields.nutrient)>=0
    @test maximum(w.fields.nutrient)<=1
    @test minimum(m.permeability)>0
    open,_,mo=medium_fixture()
    closed,_,mc=medium_fixture(barrier=100)
    for world in (open,closed)
        fill!(world.fields.nutrient,0); world.fields.nutrient[1:8,:,:].=1
    end
    closed.fields.deposit[[1,8,9,16],:,:].=1
    advance!(open,mo,1); advance!(closed,mc,1)
    @test sum(closed.fields.nutrient[9:16,:,:])<sum(open.fields.nutrient[9:16,:,:])/10
    @test isapprox(sum(closed.fields.nutrient),64;atol=1e-5)
    m.signal_source[4,4,1]=2
    advance!(w,m,.5)
    @test isapprox(sum(m.signal),2exp(-.5m.signal_decay);rtol=1e-6)
    @test all(iszero,m.signal_source)
    previous=sum(m.signal)
    advance!(w,m,.5)
    @test isapprox(sum(m.signal),previous*exp(-.5m.signal_decay);rtol=1e-6)
    m.remodeling[4,4,1]=1
    phase_mass=sum(p.phi)
    advance!(w,m,.05)
    @test all(iszero,m.remodeling)
    @test p.chemical_forcing[4,4,1]!=0
    oldphi=copy(p.phi)
    advance!(w,p,.05)
    @test maximum(abs.(p.phi.-oldphi))>0
    @test isapprox(sum(p.phi),phase_mass;atol=1e-5)
    single=World(Grid(1,1,1))
    singlephase=ImplicitPhaseSeparation(single;noise=0,mineral_pinning=0,nutrient_diffusion=0,interface_activity=0,activity_relaxation=0)
    singlemedium=ActiveMedium(single,singlephase)
    advance!(single,singlemedium,100)
    @test single.fields.nutrient[1]==1
    # Extreme finite barrier coefficients can overflow the denominator and
    # close faces entirely; zero/zero harmonic faces must remain well defined.
    closedzero=ActiveMedium(single,singlephase;membrane_barrier=floatmax(Float32))
    single.fields.deposit[1]=2
    advance!(single,closedzero,100)
    @test all(isfinite,single.fields.nutrient)
    @test closedzero.permeability[1]==0
    tiny=World(Grid(2,1,1))
    tinyphase=ImplicitPhaseSeparation(tiny;noise=0,mineral_pinning=0,nutrient_diffusion=0,interface_activity=0,activity_relaxation=0)
    tinymedium=ActiveMedium(tiny,tinyphase;membrane_barrier=floatmax(Float32))
    fill!(tiny.fields.deposit,2)
    tiny.fields.nutrient[1]=0
    advance!(tiny,tinymedium,100)
    @test vec(tiny.fields.nutrient)==Float32[0,1]
    mktempdir() do dir
        manifest=record_state(w,(m,p),joinpath(dir,"medium");step=0)
        state=manifest["processes"][1]["state"]
        @test state["membrane_alias"]=="world.fields.deposit"
        @test !haskey(state,"phase")
        for name in (:membrane,:permeability,:signal,:signal_source,:remodeling,:forcing,:scratch)
            desc=state[string(name)]
            @test Tuple(desc["shape"])==w.grid.size
            @test reinterpret(Float32,read(joinpath(dir,"medium",desc["file"])))==vec(getproperty(m,name))
        end
    end
    m.signal[1]=-1
    @test_throws ArgumentError advance!(w,m,.1)
    m.signal[1]=0
    w.fields.nutrient[1]=-1
    @test_throws ArgumentError advance!(w,m,.1)
end
