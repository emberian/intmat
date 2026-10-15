module PackingTissueChecks
using Test,Interior
include(joinpath(@__DIR__,"..","studies","living_tissue_packing.jl"))
using .LivingTissuePacking
const LP=LivingTissuePacking

function packing_stage(s,r,q,dt;h=(1f0,1f0,1f0))
    out=similar(s);rout=similar(r)
    LP.packing_kernel!(CPU(),256)(out,s,rout,r,h,Float32(dt),q,0f0,.15f0;ndrange=size(r))
    out,rout
end
@testset "packing funding and free-space limited growth" begin
    q=packing_parameters(growth=1000,renewal=0,motility=0,inhibition_diffusion=0)
    s=zeros(Float32,1,1,1,7);s[1,1,1,:].=(2.,.7,.1,0,0,0,0);r=zeros(Float32,1,1,1)
    out,ro=packing_stage(s,r,q,.08)
    @test out[1,1,1,2]>s[1,1,1,2] # Paid basal growth does not need positive x.
    @test out[1,1,1,2]+out[1,1,1,3]<=q.capacity+eps(Float32)
    @test out[1,1,1,4]>0 && out[1,1,1,1]<s[1,1,1,1]
    @test sum(Float64,out[:,:,:,1:4])≈sum(Float64,s[:,:,:,1:4]) atol=3e-7
    @test ro==r
    for x in (-.9f0,.9f0)
        s[1,1,1,:].=(0,.4,.3,0,x,0,0)
        out,ro=packing_stage(s,r,q,.08)
        @test out==s && all(iszero,ro)
        s[1,1,1,1]=1e-7
        out,ro=packing_stage(s,r,q,.08)
        @test minimum(out[:,:,:,1:4])>=0 && out[1,1,1,4]≈1f-7
        @test sum(Float64,out[:,:,:,1:4])≈sum(Float64,s[:,:,:,1:4]) atol=2e-7
    end
    @test_throws ArgumentError packing_parameters(capacity=0)
    @test_throws ArgumentError packing_parameters(basal_activity=0)
    w=LP.tissue_world(Grid(1,1,1));s[1,1,1,2]=1
    @test_throws ArgumentError PackingTissue(w;initial=s)
    @test_throws ArgumentError PackingTissue(w;base_parameters=LP.turnover_base_parameters())
end

@testset "six competing inflows respect remaining solid capacity" begin
    q=packing_parameters(capacity=1,growth=0,renewal=0,motility=1,
        construction=0,unbuilding=0,withdrawal=0,inhibition_diffusion=0,phase_work=.02)
    s=zeros(Float32,5,5,5,7);s[:,:,:,1].=1;s[:,:,:,2].=1;s[:,:,:,5].=-1
    s[3,3,3,2]=0;s[3,3,3,3]=.5;s[3,3,3,1]=8
    r=zeros(Float32,5,5,5);initial=sum(Float64,s[:,:,:,1:4]);biomass=sum(Float64,s[:,:,:,2])
    out,ro=packing_stage(s,r,q,.14) # 2*D*dt*sum(h^-2)=.84: all six faces compete.
    @test out[3,3,3,2]≈.14f0*6*.96f0*.5f0 atol=2e-7
    @test maximum(out[:,:,:,2].+out[:,:,:,3])<=1+2eps(Float32)
    @test minimum(out[:,:,:,1:4])>=0
    @test sum(Float64,out[:,:,:,2])≈biomass atol=3e-6
    @test sum(Float64,out[:,:,:,1:4])≈initial atol=1e-5
    for _ in 1:25
        s,r=packing_stage(out,ro,q,.14);out,ro=s,r
    end
    @test maximum(out[:,:,:,2].+out[:,:,:,3])<=1+4eps(Float32)
    @test minimum(out[:,:,:,1:4])>=0
    @test sum(Float64,out[:,:,:,2])≈biomass atol=2e-5
    @test sum(Float64,out[:,:,:,1:4])≈initial atol=3e-5
    # A completely occupied recipient cannot gain biomass despite a strong
    # precursor bias. Operating cost may still be paid when actions are blocked.
    s=zeros(Float32,3,3,3,7);s[:,:,:,1].=2;s[:,:,:,2].=.3;s[:,:,:,3].=.7;s[:,:,:,5].=-1
    s[2,2,2,1]=8;r=zeros(Float32,3,3,3)
    out,ro=packing_stage(s,r,q,.14)
    @test out[:,:,:,2:3]==s[:,:,:,2:3]
    @test sum(out[:,:,:,4])>0
end

@testset "hopping uses old donors and anisotropic spacing" begin
    q=packing_parameters(capacity=1,growth=0,renewal=0,motility=1,
        construction=0,unbuilding=0,withdrawal=0,inhibition_diffusion=0)
    s=zeros(Float32,5,1,1,7);s[:,:,:,1].=2;s[1,1,1,2]=.8;s[:,:,:,5].=-1;r=zeros(Float32,5,1,1)
    out,ro=packing_stage(s,r,q,.08)
    @test out[2,1,1,2]>0 && out[5,1,1,2]>0
    @test out[3,1,1,2]==out[4,1,1,2]==0 # No arriving material is spent twice.
    @test sum(out[:,:,:,2])≈.8f0 atol=1e-7
    finer,_=packing_stage(s,r,q,.08;h=(2f0,1f0,1f0))
    @test finer[2,1,1,2]≈out[2,1,1,2]/4 atol=1e-8
end

@testset "packing composition recording and exact continuation" begin
    w=LP.tissue_world(Grid(8,7,1;extent=(8,7,1)));p=PackingTissue(w;reserve=2)
    initial=packing_statistics(w,p)
    for _ in 1:400;step!(w,p,.04);end
    final=packing_statistics(w,p)
    @test final.finite && final.physical_min>=0 && final.neural_maxabs<=1
    @test final.capacity_excess<=16eps(Float32)*final.capacity
    @test final.material≈initial.material atol=1e-3
    @test final.potential+final.spent≈initial.potential+initial.spent atol=1e-3
    @test final.reserve<initial.reserve && final.spent>0
    @test Interior.experiment_diagnostics(w,(p,)).material_mass≈final.material
    before=copy(p.tissue.state);reserve=copy(p.reserve);t=w.time
    @test_throws ArgumentError step!(w,(p,p.flow),.04)
    @test_throws ArgumentError step!(w,(p,p.tissue),.04)
    @test p.tissue.state==before && p.reserve==reserve && w.time==t
    @test_throws ArgumentError step!(w,p,1.)
    mktempdir() do dir
        path=joinpath(dir,"checkpoint.jls");packing_checkpoint(w,p,path)
        v,q=restore_packing!(path)
        @test q.tissue.state==p.tissue.state && q.reserve==p.reserve && v.time==w.time
        for _ in 1:11;step!(w,p,.04);step!(v,q,.04);end
        @test q.tissue.state==p.tissue.state && q.reserve==p.reserve && v.time==w.time
        manifest=record_state(w,p,joinpath(dir,"state");step=411)
        @test manifest["arrays"]["process-1-reserve"]["shape"]==[8,7,1]
    end
end
@testset "optional paid basal transport and matrix remodeling" begin
    q=packing_parameters(capacity=1,growth=0,renewal=0,motility=1,
        construction=0,unbuilding=1,withdrawal=0,inhibition_diffusion=0,
        basal_transport=.2,basal_remodeling=.1)
    s=zeros(Float32,3,1,1,7);s[:,:,:,1].=2;s[1,1,1,2]=.6
    r=zeros(Float32,3,1,1);out,ro=packing_stage(s,r,q,.08)
    @test out[2,1,1,2]>0 && out[3,1,1,2]>0 # Quiescent paid tissue can relocate.
    @test sum(out[:,:,:,2])≈.6f0 atol=1e-7
    s=zeros(Float32,1,1,1,7);s[1,1,1,:].=(1,0,.5,0,0,0,0);r=zeros(Float32,1,1,1)
    out,ro=packing_stage(s,r,q,.08)
    @test out[1,1,1,3]<s[1,1,1,3] && out[1,1,1,4]>0
    @test sum(Float64,out[:,:,:,1:4])≈sum(Float64,s[:,:,:,1:4]) atol=1e-7
    s[1,1,1,1]=0;out,ro=packing_stage(s,r,q,.08)
    @test out==s # Basal actions also stop without paid precursor/reserve.
    q=packing_parameters(capacity=1,growth=0,renewal=0,motility=1,
        construction=0,unbuilding=0,withdrawal=0,inhibition_diffusion=0,basal_transport=100)
    s=zeros(Float32,5,5,5,7);s[:,:,:,1].=2;s[:,:,:,2].=1;s[:,:,:,5].=-1
    s[3,3,3,2]=0;s[3,3,3,1]=8;r=zeros(Float32,5,5,5)
    out,ro=packing_stage(s,r,q,.14)
    @test maximum(out[:,:,:,2])<=1+4eps(Float32)
    @test minimum(out[:,:,:,1:4])>=0
end
@testset "joint biomass and matrix inflows share free occupancy" begin
    q=packing_parameters(capacity=1,growth=0,renewal=0,motility=1,matrix_transport=1,
        construction=0,unbuilding=0,withdrawal=0,inhibition_diffusion=0,phase_work=.02)
    s=zeros(Float32,5,5,5,7);s[:,:,:,1].=2;s[:,:,:,2].=.4;s[:,:,:,3].=.6;s[:,:,:,5].=-1
    s[3,3,3,2:3].=0;s[3,3,3,1]=8;r=zeros(Float32,5,5,5)
    initial_b=sum(Float64,s[:,:,:,2]);initial_m=sum(Float64,s[:,:,:,3]);initial=sum(Float64,s[:,:,:,1:4])
    out,ro=packing_stage(s,r,q,.14)
    @test out[3,3,3,2]≈.14f0*6*.96f0*.4f0 atol=2e-7
    @test out[3,3,3,3]≈.14f0*6*.6f0 atol=2e-7
    @test maximum(out[:,:,:,2].+out[:,:,:,3])<=1+4eps(Float32)
    @test minimum(out[:,:,:,1:4])>=0
    @test sum(Float64,out[:,:,:,2])≈initial_b atol=3e-6
    @test sum(Float64,out[:,:,:,3])≈initial_m atol=3e-6
    for _ in 1:25;out,ro=packing_stage(out,ro,q,.14);end
    @test maximum(out[:,:,:,2].+out[:,:,:,3])<=1+4eps(Float32)
    @test minimum(out[:,:,:,1:4])>=0
    @test sum(Float64,out[:,:,:,2])≈initial_b atol=2e-5
    @test sum(Float64,out[:,:,:,3])≈initial_m atol=2e-5
    @test sum(Float64,out[:,:,:,1:4])≈initial atol=3e-5
    # Every recipient is full: neither species can use a second capacity budget.
    s=zeros(Float32,3,3,3,7);s[:,:,:,1].=2;s[:,:,:,2].=.4;s[:,:,:,3].=.6;s[:,:,:,5].=-1
    r=zeros(Float32,3,3,3);out,ro=packing_stage(s,r,q,.14)
    @test out[:,:,:,2:3]==s[:,:,:,2:3]
    w=LP.tissue_world(Grid(5,5,5;extent=(5,5,5)))
    q=packing_parameters(motility=0,matrix_transport=10,inhibition_diffusion=0)
    p=PackingTissue(w;parameters=q)
    @test stability_limit(w,p)<.02
end
@testset "paid matrix donors do not spend incoming gel twice" begin
    q=packing_parameters(capacity=1,growth=0,renewal=0,motility=0,matrix_transport=1,
        construction=0,unbuilding=0,withdrawal=0,inhibition_diffusion=0)
    s=zeros(Float32,5,1,1,7);s[:,:,:,1].=2;s[1,1,1,3]=.6;r=zeros(Float32,5,1,1)
    out,ro=packing_stage(s,r,q,.08)
    @test out[2,1,1,3]>0 && out[5,1,1,3]>0
    @test out[3,1,1,3]==out[4,1,1,3]==0
    @test sum(out[:,:,:,3])≈.6f0 atol=1e-7
    s[1,1,1,1]=0;out,ro=packing_stage(s,r,q,.08)
    @test out[:,:,:,3]==s[:,:,:,3] # Unfunded matrix cannot launch an exchange.
end
end # module PackingTissueChecks
