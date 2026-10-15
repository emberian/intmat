using Test,Interior
include(joinpath(@__DIR__,"..","studies","living_tissue_turnover.jl"))
using .LivingTissueTurnover
const LT=LivingTissueTurnover
include(joinpath(@__DIR__,"..","studies","living_tissue_observations.jl"))
using .TissueObservations

@testset "reformation requires an observed opening" begin
    p=TissueTurnoverObservation(fill(.239f0,1,1))
    @test observe_turnover!(p,fill(.3f0,1,1)).reformation_events==0
    @test observe_turnover!(p,fill(.09f0,1,1)).reopening_events==1
    @test observe_turnover!(p,fill(.25f0,1,1)).reformation_events==1
    @test observe_turnover!(p,fill(.3f0,1,1)).reformation_events==1
    @test observe_turnover!(p,fill(.09f0,1,1)).reopening_events==2
    @test observe_turnover!(p,fill(.2f0,1,1)).reformation_events==1
    @test observe_turnover!(p,fill(.25f0,1,1)).reformation_events==2
end

@testset "paid phase transfers and conservative relocation" begin
    for excitation in (-.6f0,.6f0)
        s=zeros(Float32,1,1,1,7);s[1,1,1,:].=(.8,.5,.4,.1,excitation,0,0)
        reserve=ones(Float32,1,1,1);out=similar(s);rout=similar(reserve)
        q=turnover_parameters(renewal=0,motility=0,inhibition_diffusion=0)
        dt=.04f0
        LT.active_turnover_kernel!(CPU(),256)(out,s,rout,reserve,(1f0,1f0,1f0),dt,q,8f0;ndrange=(1,1,1))
        phase=abs(excitation)-q.phase_threshold
        work=dt*q.phase_work*(s[1,1,1,2]+s[1,1,1,3])*phase
        @test out[1,1,1,4]≈s[1,1,1,4]+work atol=1e-7
        @test sum(out[:,:,:,1:4])≈sum(s[:,:,:,1:4]) atol=2e-7
        @test rout==reserve
        if excitation>0
            @test out[1,1,1,3]>s[1,1,1,3] && out[1,1,1,2]<s[1,1,1,2]
        else
            @test out[1,1,1,1]>s[1,1,1,1] && out[1,1,1,3]<s[1,1,1,3] && out[1,1,1,2]<s[1,1,1,2]
        end
        s[1,1,1,1]=0;reserve.=0
        LT.active_turnover_kernel!(CPU(),256)(out,s,rout,reserve,(1f0,1f0,1f0),dt,q,8f0;ndrange=(1,1,1))
        @test out==s && all(iszero,rout)
    end
    s=zeros(Float32,2,1,1,7);s[:,1,1,1].=(.2,1.);s[1,1,1,2]=.8;s[:,1,1,5].=-.6
    reserve=zeros(Float32,2,1,1);out=similar(s);rout=similar(reserve)
    q=turnover_parameters(renewal=0,motility=1,construction=0,unbuilding=0,withdrawal=0,inhibition_diffusion=0,phase_threshold=0)
    LT.active_turnover_kernel!(CPU(),256)(out,s,rout,reserve,(1f0,1f0,1f0),.04f0,q,0f0;ndrange=(2,1,1))
    @test sum(out[:,:,:,2])≈.8f0 atol=1e-7
    @test out[2,1,1,2]>.037 && minimum(out[:,:,:,1:4])>=0
    s[1,1,1,1]=0
    LT.active_turnover_kernel!(CPU(),256)(out,s,rout,reserve,(1f0,1f0,1f0),.04f0,q,0f0;ndrange=(2,1,1))
    @test out[:,:,:,2]==s[:,:,:,2] # Unfunded donors cannot relocate.
end

@testset "physical initialization is consistent under refinement" begin
    coarse=Grid(8,8,1;extent=(24,24,1));fine=Grid(16,16,1;extent=(24,24,1))
    a=turnover_initial_state(coarse,921;noise_model="bandlimited")
    b=turnover_initial_state(fine,921;noise_model="bandlimited")
    for c in (1,2,5)
        restricted=[sum(b[2i-1:2i,2j-1:2j,1,c])/4 for i in 1:8,j in 1:8]
        @test maximum(abs.(restricted.-a[:,:,1,c]))<2e-7
    end
    w=tissue_world(coarse);p=TurnoverTissue(w;initial=a,reserve=4)
    v=tissue_world(fine);q=TurnoverTissue(v;initial=b,reserve=4)
    @test turnover_statistics(w,p).physical_material≈turnover_statistics(v,q).physical_material atol=1e-4
    @test turnover_statistics(w,p).physical_reserve≈turnover_statistics(v,q).physical_reserve
    coarse=Grid(4,4,4;extent=(24,24,24));fine=Grid(8,8,8;extent=(24,24,24))
    a=turnover_initial_state(coarse,921;noise_model="bandlimited")
    b=turnover_initial_state(fine,921;noise_model="bandlimited")
    for c in (1,2,5)
        restricted=[sum(b[2i-1:2i,2j-1:2j,2k-1:2k,c])/8 for i in 1:4,j in 1:4,k in 1:4]
        @test maximum(abs.(restricted.-a[:,:,:,c]))<2e-7
    end
end

@testset "finite renewal, composition and continuation" begin
    w=tissue_world(Grid(8,7,1;extent=(8,7,1)));p=TurnoverTissue(w)
    initial=turnover_statistics(w,p)
    for _ in 1:300;step!(w,p,.04);end
    final=turnover_statistics(w,p)
    @test final.finite && final.physical_min>=0 && final.neural_maxabs<=1
    @test final.material≈initial.material atol=1e-3
    @test final.potential+final.spent≈initial.potential+initial.spent atol=1e-3
    @test final.reserve<initial.reserve && final.spent>initial.spent
    @test final.matrix>0
    @test Interior.experiment_diagnostics(w,(p,)).material_mass≈final.material
    before=copy(p.tissue.state);reserve=copy(p.reserve);t=w.time
    @test_throws ArgumentError step!(w,(p,p.tissue),.04)
    @test p.tissue.state==before && p.reserve==reserve && w.time==t
    @test_throws ArgumentError TurnoverTissue(w;reserve=-1)
    @test_throws ArgumentError turnover_parameters(phase_work=0)
    mktempdir() do dir
        path=joinpath(dir,"checkpoint.jls");turnover_checkpoint(w,p,path)
        v,q=restore_turnover!(path)
        @test v.time==w.time && q.reserve==p.reserve && q.tissue.state==p.tissue.state
        for _ in 1:11;step!(w,p,.04);step!(v,q,.04);end
        @test v.time==w.time && q.reserve==p.reserve && q.tissue.state==p.tissue.state
        manifest=record_state(w,p,joinpath(dir,"state");step=311)
        @test manifest["arrays"]["process-1-reserve"]["shape"]==[8,7,1]
    end
end
