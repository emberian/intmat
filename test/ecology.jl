@testset "heritable finite ecology" begin
    w=World(Grid(8,8;extent=(8,8,1)))
    p=Ecology(Float32[1 3; 1 3; 0 0];storage=2,energy=5,mutation=0.1,capacity=8)
    p.weights[end,:,:].=0; p.weights[end,end,:].=2
    m=ecology_mass(w,p)
    step!(w,p,0.1)
    @test p.births==2
    @test p.parents[3:4]==[1,2]
    @test p.lineages[3:4]==[1,2]
    @test p.generations[3:4]==[1,1]
    @test p.weights[:,:,3]!=p.weights[:,:,1]
    @test ecology_mass(w,p)≈m rtol=1e-6
    @test sum(Float64,p.energy)+p.heat≈p.energy_input rtol=1e-6
    q=deepcopy(p); v=deepcopy(w)
    reorder_ecology!(q,reverse(eachindex(q.ids)))
    step!(w,p,0.1); step!(v,q,0.1)
    @test p.ids==q.ids && p.positions==q.positions && p.weights==q.weights
    p.capacity=length(p.ids); p.storage.=2; p.energy.=5
    before=ecology_mass(w,p)
    step!(w,p,0.1)
    @test p.blocked_births>0
    @test ecology_mass(w,p)≈before rtol=1e-6
    p.age.=p.lifespan
    before=ecology_mass(w,p)
    step!(w,p,0.1)
    @test isempty(p.ids)
    @test ecology_mass(w,p)≈before rtol=1e-6
end
@testset "ecology inheritance and invalid controls" begin
    w=World(Grid(4,4;extent=(4,4,1)))
    p=Ecology(zeros(Float32,3,1);storage=2,energy=5,mutation=0,capacity=2)
    p.weights[end,:,:].=0; p.weights[end,end,:].=2
    founder=copy(p.weights[:,:,1])
    step!(w,p,0.1)
    @test p.weights[:,:,2]==founder
    @test p.hidden[:,2]==p.hidden[:,1]
    p.uptake_rate=-1
    before=copy(w.fields.nutrient)
    @test_throws ArgumentError step!(w,p,0.1)
    @test before==w.fields.nutrient
end
@testset "ecology controller uses physical time" begin
    base=Ecology(zeros(Float32,3,1);uptake_rate=0,basal_cost=0,work_rate=0,secretion_rate=0,speed=0,response_time=1)
    fill!(base.weights,0)
    for j in 1:size(base.hidden,1)
        base.weights[j,j,1]=0.9f0
        base.weights[j,end,1]=0.2f0
    end
    function after(dt)
        p=deepcopy(base); w=World(Grid(4,4))
        for k in 1:round(Int,0.1/dt); step!(w,p,dt); end
        p.hidden[1,1]
    end
    coarse=after(0.1); medium=after(0.01); fine=after(0.001)
    @test abs(medium-fine)<abs(coarse-fine)
    @test abs(coarse-fine)<0.002
    @test fine<0.025
end
