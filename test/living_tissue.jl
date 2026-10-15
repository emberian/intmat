using Test, Interior
include(joinpath(@__DIR__,"..","studies","living_tissue.jl"))
using .LivingTissueStudy

struct TissuePreflightGuard end
Interior.stability_limit(w,::TissuePreflightGuard)=.001
Interior.advance!(w,::TissuePreflightGuard,dt)=error("preflight should reject first")

@testset "continuum tissue inventory, funding and recording" begin
    w=tissue_world(Grid(12,10,1;extent=(12,10,1)));p=LivingTissue(w)
    initial=tissue_statistics(w,p)
    experiment=Experiment(w,(tissue=p,);dt=.04)
    report=run!(experiment;steps=500,capture_every=500)
    final=tissue_statistics(w,p)
    @test final.finite && final.physical_min>=0 && final.neural_maxabs<=1
    @test final.material≈initial.material atol=2e-4
    @test final.potential+final.spent≈initial.potential atol=2e-4
    @test final.spent>0 && final.matrix>0
    @test report.samples[end].diagnostics.material_mass≈final.material
    @test isempty(Interior.process_ports(p).writes)
    @test record_process_fields(p).state===p.state
    before=copy(p.state);t=w.time
    @test_throws ArgumentError step!(w,(p,TissuePreflightGuard()),.04)
    @test before==p.state && w.time==t
    @test_throws ArgumentError step!(w,p,1.)
    @test_throws ArgumentError tissue_parameters(metabolic_work=0)
    @test_throws ArgumentError tissue_parameters(activity_control=2)
    @test_throws ArgumentError LivingTissue(w;initial=fill(-1f0,12,10,1,7))
    mktempdir() do dir
        manifest=record_state(w,p,joinpath(dir,"recording");step=500)
        @test manifest["arrays"]["process-1-state"]["shape"]==[12,10,1,7]
        path=joinpath(dir,"checkpoint.jls");tissue_checkpoint(w,p,path)
        v,q=restore_tissue!(path)
        @test q.state==p.state && v.time==w.time
        for _ in 1:17;step!(w,p,.04);step!(v,q,.04);end
        @test q.state==p.state && v.time==w.time
        @test_throws ArgumentError tissue_checkpoint(w,p,path)
    end
end

@testset "symmetric transport and finite fuel" begin
    w=tissue_world(Grid(5,4,1;extent=(5,4,1)))
    state=zeros(Float32,5,4,1,7)
    state[:,:,1,1].=reshape(Float32.(1:20),5,4)./20
    state[2,2,1,3]=.8f0
    p=LivingTissue(w;initial=state)
    n0=copy(state[:,:,:,1]);m0=copy(state[:,:,:,3]);dt=.04f0;q=p.parameters
    expected=copy(n0)
    for j in 1:4,i in 1:5
        for (a,b) in ((mod1(i-1,5),j),(mod1(i+1,5),j),(i,mod1(j-1,4)),(i,mod1(j+1,4)))
            face=1f0/(1f0+4f0*q.barrier*(m0[i,j,1]+m0[a,b,1]))
            expected[i,j,1]+=dt*q.precursor_diffusion*face*(n0[a,b,1]-n0[i,j,1])
        end
        # Matrix erosion returns precursor in the same local update.
        expected[i,j,1]+=m0[i,j,1]*(dt*q.erosion)/(1f0+dt*q.erosion)
    end
    step!(w,p,dt)
    @test maximum(abs.(p.state[:,:,:,1].-expected))<2e-7
    @test all(iszero,p.state[:,:,:,4]) # No biomass means no metabolic charge.
    # With no precursor and no recycling, neither excitation nor construction
    # can obtain work despite an initially positive activator.
    s=zeros(Float32,1,1,1,7);s[1,1,1,2]=1;s[1,1,1,5]=.6
    q=tissue_parameters(turnover=0,starvation_turnover=0,growth=0)
    z=tissue_world(Grid(1,1,1));r=LivingTissue(z;initial=s,parameters=q)
    for _ in 1:100;step!(z,r,.04);end
    @test r.state[1,1,1,2]==1 && r.state[1,1,1,3]==0 && r.state[1,1,1,4]==0
    @test 0<r.state[1,1,1,5]<s[1,1,1,5]
    @test r.state[1,1,1,7]>0 # Passive local memory can outlast a drive.
    # A two-site periodic domain has two shared faces. Donor-directed transport
    # moves biomass toward the precursor-rich site without creating biomass.
    s=zeros(Float32,2,1,1,7);s[1,1,1,1]=.2;s[2,1,1,1]=1;s[1,1,1,2]=1
    q=tissue_parameters(precursor_diffusion=0,biomass_diffusion=0,excitation_diffusion=0,
        taxis=1,turnover=0,starvation_turnover=0,growth=0,construction=0,erosion=0)
    z=tissue_world(Grid(2,1,1;extent=(2,1,1)));r=LivingTissue(z;initial=s,parameters=q)
    step!(z,r,.04)
    @test sum(r.state[:,:,:,2])≈1 atol=1e-7
    @test r.state[2,1,1,2]>.079 && minimum(r.state[:,:,:,1:4])>=0
end

@testset "delayed inhibition cycles at a fixed site" begin
    q=tissue_parameters(growth=0,turnover=0,starvation_turnover=0,
        construction=0,erosion=0,matrix_feedback=0,neural_gain=2.2)
    s=zeros(Float32,1,1,1,7);s[1,1,1,1]=100;s[1,1,1,2]=1;s[1,1,1,5]=.1
    w=tissue_world(Grid(1,1,1));p=LivingTissue(w;initial=s,parameters=q)
    crossing=0;previous=1f0;lo=1f0;hi=-1f0
    for k in 1:3000
        step!(w,p,.04)
        if k>1000
            x=p.state[1,1,1,5];crossing+=x*previous<0;previous=x
            lo=min(lo,x);hi=max(hi,x)
        end
    end
    @test crossing>=4 && lo<-.2 && hi>.2
    @test p.state[1,1,1,4]>0
    # Refinement checks the complete coupled law over a short fixed duration.
    a=tissue_world(Grid(1,1,1));b=tissue_world(a.grid)
    p=LivingTissue(a;initial=s,parameters=q);r=LivingTissue(b;initial=s,parameters=q)
    for _ in 1:100;step!(a,p,.04);end
    for _ in 1:200;step!(b,r,.02);end
    @test maximum(abs.(p.state.-r.state))<.015
end
