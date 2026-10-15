using Interior, Test, KernelAbstractions
isdefined(Interior,:DeviceCellular) || Base.include(Interior,joinpath(@__DIR__,"..","src","device_cellular.jl"))

function energy_audit_motor(T;energy=1.)
    (position=zeros(Float32,1,3),position_next=zeros(Float32,1,3),active=Int32[1],
     energy=T[energy],actions=fill(.2f0,1,11),force=zeros(Float32,1,3),rows=zeros(Float32,1),
     ledger=zeros(Float32,1,9),ledger_low=zeros(Float32,1,9))
end
energy_audit_work(s)=Float64(s.ledger[1,1])+Float64(s.ledger_low[1,1])
function energy_audit_fixture(;energy_type=Float32,kwargs...)
    w=World(Grid(4,4;extent=(4,4,1)));fill!(w.fields.nutrient,0)
    opts=merge((;capacity=4,hidden=2,energy_type,biomass=1.6,energy=1.,metabolites=0.,
        initial=(0.,0.,0.),diffusion=(0.,0.,0.),speed=0.,uptake_rate=0.,growth_rate=0.,
        catabolic_rate=0.,basal_rate=0.,neural_cost=0.,repulsion=0.,stiffness=0.,
        mutation=0.,developmental_partition=0.,mechanical_substeps=1), (;kwargs...))
    w,Interior.device_cellular_initialize(w,reshape([2.,2.,0.],3,1);opts...)
end

@testset "Actual motor microcost survives selective energy precision" begin
    e=(status=Int32[0],);q=(speed=.06f0,motor_cost=.5f0);g=(nd=Int32(3),extent=(4f0,4f0,4f0))
    totals=Float64[]
    for substeps in (1,16,64)
        s=energy_audit_motor(Float64);old=energy_audit_motor(Float32)
        dt=.001f0/Float32(substeps)
        for _ in 1:1100*substeps
            Interior.dc_move_at!(s,e,q,g,dt,1);s.position.=s.position_next
            Interior.dc_move_at!(old,e,q,g,dt,1);old.position.=old.position_next
        end
        work=energy_audit_work(s);push!(totals,work)
        @test s.energy[1]<1
        @test s.energy[1]+work≈1 atol=2e-12 rtol=0
        @test s.position==old.position # Same motion while funding is not binding.
        @test work>0
        substeps>=16 && @test old.energy[1]==1f0 # The original demonstrated failure.
    end
    @test maximum(totals)-minimum(totals)<2e-10
    # A sub-Float32-ULP battery must still fund affordable fractional movement.
    s=energy_audit_motor(Float64;energy=1e-12);before=s.energy[1]
    Interior.dc_move_at!(s,e,q,g,.001f0/16,1)
    @test 0<=s.energy[1]<before
    @test maximum(abs,s.position_next)>0
    @test s.energy[1]+energy_audit_work(s)≈before atol=1e-24 rtol=0
    @test e.status==[0]
end

@testset "Actual maintenance endpoints and legacy default" begin
    w,p=energy_audit_fixture(;energy_type=Float64,basal_rate=.008,neural_cost=.0002)
    for _ in 1:1100
        Interior.dc_metabolism!(CPU(),1)(p.state,p.environment,p.parameters,p.geometry,nothing,.001f0;ndrange=4)
    end
    @test p.state.energy[1]+energy_audit_work(p.state)≈1 atol=2e-12 rtol=0
    @test eltype(p.state.metabolites)==eltype(p.state.mass)==Float32
    w,a=energy_audit_fixture();v,b=energy_audit_fixture(;energy_type=Float32)
    @test a.geometry.energy_bits==32
    for _ in 1:3
        Interior.device_cellular_step!(w,a,.001f0)
        Interior.device_cellular_step!(v,b,.001f0)
    end
    @test Interior.device_cellular_snapshot(a)==Interior.device_cellular_snapshot(b)
end

@testset "Respiration credit and growth preserve intended material types" begin
    for energy in (0.,1.,prevfloat(Float64(floatmax(Float32)))),released in (1f-8,.01f0)
        after,heat=Interior.dc_energy_credit(energy,released,.65f0)
        @test after-energy+heat==Float64(released)
        @test heat>=0 && after>=energy
    end
    for energy in (1.,1e-12)
        a,b,mass,remaining,work,heat,built=Interior.device_anabolic_build!(nothing,1,.2f0,.3f0,1f0,energy,.01f0,.4f0,4f0,4f0)
        @test a isa Float32 && b isa Float32 && mass isa Float32 && built isa Float32
        @test remaining>=0 && isapprox(remaining+work,energy;atol=1e-16)
        @test built<=.01f0
    end
end

@testset "Wider battery division death and strict continuation" begin
    w,p=energy_audit_fixture(;energy_type=Float64,biomass=2.)
    p.state.energy[1]=1+2.0^-40;p.state.actions[1,10]=1
    initial=p.state.energy[1]
    Interior.dc_lifecycle!(CPU(),1)(p.state,p.environment,p.parameters,p.geometry,UInt32(1);ndrange=1)
    Interior.dc_divide!(CPU(),1)(p.state,p.environment,p.parameters,p.geometry,UInt32(1);ndrange=4)
    @test count(==(1),p.state.active)==2
    @test sum(p.state.energy)==initial
    @test p.state.energy[1]==p.state.energy[2]
    # A normally completed step rebuilds geometric scratch before checkpointing.
    Interior.device_cellular_step!(w,p,.001f0)
    cp=Interior.device_cellular_checkpoint(w,p)
    v=World(w.grid);fill!(v.fields.nutrient,0)
    restored=Interior.device_cellular_restore(v,cp)
    @test eltype(restored.state.energy)==Float64
    @test restored.geometry.energy_bits==64
    @test Interior.device_cellular_snapshot(restored)==Interior.device_cellular_snapshot(p)
    for _ in 1:3
        Interior.device_cellular_step!(w,p,.001f0)
        Interior.device_cellular_step!(v,restored,.001f0)
    end
    @test Interior.device_cellular_snapshot(restored)==Interior.device_cellular_snapshot(p)
    w,dead=energy_audit_fixture(;energy_type=Float64,biomass=.26,energy=1e-12,
        basal_rate=1.,starvation_rate=100.)
    initial=dead.state.energy[1]
    Interior.dc_metabolism!(CPU(),1)(dead.state,dead.environment,dead.parameters,dead.geometry,nothing,.2f0;ndrange=4)
    @test dead.state.active[1]==0 && dead.state.energy[1]==0
    @test energy_audit_work(dead.state)≈initial atol=1e-25 rtol=0
end
