using Test
include(joinpath(@__DIR__,"../studies/native_temporal_living.jl"))
const V=NativeTemporalLiving
@testset "continuing temporal living policy" begin
    for balance in (-1f0,-.2f0,0f0,.2f0,1f0),trend in (-.1f0,0f0,.1f0)
        a=V.feature_bank(balance,trend,false);b=V.feature_bank(balance,trend,true)
        @test sum(abs2,a)≈2f0 atol=2eps(Float32)
        @test sum(abs2,a)==sum(abs2,b)
        @test a[1:2]==b[1:2]
        @test V.feature_bank(balance,trend,true,true)==b
    end
    profile=joinpath(@__DIR__,"../experiments/native/living-temporal-allocation-16.toml")
    cfg=V.TOML.parsefile(profile);cfg["cellular"]["biomass"]=1.79
    birthprofile=tempname()*".toml";open(io->V.TOML.print(io,cfg),birthprofile,"w")
    f=V.N.fixture(birthprofile);g=deepcopy(f)
    h=V.sidecar(f.p;rate=0.,learning=false,critic=false,sigma=.3,act_time=4.,inheritance_fraction=.5)
    j=V.sidecar(g.p;rate=0.,learning=false,critic=false,sigma=.3,act_time=4.,reactive=true,inheritance_fraction=.5)
    @test size(h.state.learned,2)==3
    @test h.options.critic_scale==4f0
    V.L.advance!(f.w,f.p,h,1.;dt=.02);V.L.advance!(g.w,g.p,j,1.;dt=.02)
    a=V.Interior.device_habitat_snapshot(f.p);b=V.Interior.device_habitat_snapshot(g.p)
    @test isequal(a.cells.state,b.cells.state)
    @test isequal(a.arrays,b.arrays)
    @test isequal(a.feedstock,b.feedstock)
    @test count(==(1),a.cells.state.active)>24
    @test h.temporal.ids==f.p.cells.state.ids
    @test h.temporal.fast==j.temporal.fast
    @test all(h.state.drive[:,1][a.cells.state.active.==1].==h.options.investment_logit)
    @test all(h.state.drive[:,3][a.cells.state.active.==1].==h.options.reserve_logit)
    @test all(h.state.queue_ids.==0)
    @test all(h.state.paid.==0)
    # An independently copied sidecar owns the temporal history, including its clock.
    copied=deepcopy((w=f.w,p=f.p,h=h))
    V.L.advance!(f.w,f.p,h,.2;dt=.02);V.L.advance!(copied.w,copied.p,copied.h,.2;dt=.02)
    @test h.temporal.fast==copied.h.temporal.fast
    @test h.temporal.last_step==copied.h.temporal.last_step
    @test isequal(V.Interior.device_habitat_snapshot(f.p),V.Interior.device_habitat_snapshot(copied.p))
end
@testset "independent observation and held-action clocks" begin
    profile=joinpath(@__DIR__,"../experiments/native/living-temporal-construction-16.toml")
    cfg=V.TOML.parsefile(profile);cfg["cellular"]["biomass"]=1.79
    path=tempname()*".toml";open(io->V.TOML.print(io,cfg),path,"w")
    f=V.N.fixture(path)
    a=V.C.sidecar(f.p);b=V.C.sidecar(f.p;decision_time=.2)
    @test isequal(a.options,b.options)
    @test isequal(a.state,b.state)
    h=V.sidecar(f.p;decision_time=1.,observation_time=.2,horizon=4.,act_time=12.,reserve=.02,rate=.1)
    @test h.options.decision_stride==50
    @test h.options.observation_stride==10
    @test h.options.queue_length==4
    V.L.advance!(f.w,f.p,h,.02;dt=.02)
    ids=copy(f.p.cells.state.ids);drive=copy(h.state.drive);fast=copy(h.temporal.fast)
    V.L.advance!(f.w,f.p,h,.4;dt=.02)
    same=(ids.==f.p.cells.state.ids).&(ids.!=0)
    @test h.state.drive[same,:]==drive[same,:]
    @test h.temporal.fast!=fast
    @test all(h.temporal.last_step[same].==20)
    copied=deepcopy((w=f.w,p=f.p,h=h))
    V.L.advance!(f.w,f.p,h,8.58;dt=.02)
    V.L.advance!(copied.w,copied.p,copied.h,8.58;dt=.02)
    @test isequal(V.Interior.device_habitat_snapshot(f.p),V.Interior.device_habitat_snapshot(copied.p))
    @test h.temporal==copied.h.temporal
    @test h.state==copied.h.state
    @test sum(h.state.matured)>0
    @test sum(h.state.actor_updates)>0
    @test sum(h.state.paid)>0
    @test count(==(1),f.p.cells.state.active)>24
    @test all(h.temporal.ids[f.p.cells.state.active.==1].==f.p.cells.state.ids[f.p.cells.state.active.==1])
end
@testset "whole initial cohort survival accounting" begin
    initial=(cells=[(id=1,usable=2.),(id=2,usable=3.)],)
    final=(cells=[(id=2,usable=4.)],)
    events=[(parent_id=1,child_usable_before_copy=.5,time=1.)]
    r=V.returns(initial,final,events,0.)
    @test r[1].alive==false && r[1].value===nothing
    @test r[1].all_initial_value==-1.5
    @test r[2].all_initial_value==r[2].value==1.
end

@test_throws ErrorException V.frozen_probe(nothing,nothing,(backend=nothing,),nothing,1.,:learned)
