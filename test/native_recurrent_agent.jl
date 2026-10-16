using Test
include("../studies/native_recurrent_agent.jl")
const A=NativeRecurrentAgent
@testset "native recurrent adapter" begin
    e=A.make_env(joinpath(@__DIR__,"../experiments/native/living-temporal-construction-16.toml"))
    @test size(e.x,1)==25
    @test_throws ErrorException A.make_env("unused";lower=(-.1f0,0f0,.02f0))
    @test_throws ErrorException A.make_env("unused";phase_time=Inf)
    @test_throws ErrorException A.act!(e,zeros(Float64,3,size(e.x,2)))
    @test all(iszero,e.incoming)
    old=copy(e.incoming);A.observe!(e);@test e.incoming==old
    @test A.decode(A.INITIAL_BIAS[3],A.LOWER[3],A.UPPER[3])≈.035f0
    @test A.Interior.device_controller_override!(nothing,e.w,e.p.cells,.02)===nothing
    # The observer/wrapper without installed commands preserves every native cell array.
    a=deepcopy(e);b=deepcopy(e)
    A.L.advance!(a.w,a.p,a.h,.2;dt=.02)
    A.L.advance!(b.w,b.p,b.h.habitat,.2;dt=.02)
    @test all(getproperty(a.p.cells.state,k)==getproperty(b.p.cells.state,k) for k in keys(a.p.cells.state))
    raw=zeros(Float32,3,size(e.x,2));raw[3,:].=A.INITIAL_BIAS[3]
    A.act!(e,raw);before=copy(e.start_value);A.advance!(e)
    @test sum(e.funded)>0
    @test sum(e.outgoing)==24
    @test e.target==e.x[1:12,:]
    @test all(e.reward[i]≈Float32(A.value(e.p.cells.state,e.p.cells.parameters,i)-before[i]+e.gifts[i]) for i in 1:24)
    # All13 mappings use the shared native controller/regulator boundary.
    h=e.h;p=e.p.cells;g=p.geometry
    for k in 1:13;h.agent.commands[1,k]=k==11 ? .1f0 : .01f0*k;end
    A.override!(h.backend,64)(p.state,p.environment.status,g,h.agent,Tuple(1:13);ndrange=Int(g.capacity))
    @test p.state.outputs[1,Int(g.hidden)+1:Int(g.hidden)+10]≈Float32.(.01:.01:.10)
    @test p.state.outputs[1,Int(g.hidden)+11]≈-.8f0
    A.Interior.device_anabolic_prepare!(h,e.w,p,.02)
    @test h.habitat.arrays.expression[1,1]≈.12f0
    @test h.habitat.arrays.expression[1,2]≈.13f0
    # Exact identity transfer is once per child; no slot-reuse bootstrap.
    c=deepcopy(e);s=c.p.cells.state;s.new_parent[2]=1;s.ids[2]=3000
    c.h.agent.gift_recorded[3000]=0
    beforegift=c.gifts[1];gift=A.value(s,c.p.cells.parameters,2)
    A.birth_gift!(c.h.backend,64)(s,c.p.cells.environment.status,c.p.cells.parameters,c.h.agent;ndrange=Int(g.capacity))
    @test c.gifts[1]≈beforegift+gift
    A.birth_gift!(c.h.backend,64)(s,c.p.cells.environment.status,c.p.cells.parameters,c.h.agent;ndrange=Int(g.capacity))
    @test c.gifts[1]≈beforegift+gift
    A.outcome!(c.h.backend,64)(s,c.p.cells.environment.status,c.p.cells.parameters,c.h.agent;ndrange=Int(g.capacity))
    @test c.outgoing[2]==0
    @test c.reward[2]==Float32(-c.start_value[2]+c.gifts[2])
    c.h.agent.commands[2,12]=.75f0;c.h.agent.funded[2]=9.;c.h.agent.funded_ids[2]=2
    c.h.habitat.arrays.last_secretion[2,:].=0
    A.track_funded!(c.h.backend,64)(s,c.p.cells.environment.status,c.h.habitat.arrays,c.h.agent;ndrange=Int(g.capacity))
    c.h.steps+=1;A.observe!(c)
    @test c.x[24,2]==0 && c.x[25,2]==0
    @test c.funded[2]==0 && c.x[12,2]==0
    @test c.incoming[2]==0
    frame=A.summary(e)
    @test abs(frame.ledger.material_residual)<1e-8
    @test abs(frame.ledger.potential_residual)<1e-8
end
