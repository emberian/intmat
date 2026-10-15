function device_protocol_fixture(;compact=false)
    grid=Grid(4,4;extent=(4,4,1))
    w=compact ? World(grid;fields=(),scratch_fields=()) : World(grid)
    hasproperty(w.fields,:nutrient) && fill!(w.fields.nutrient,0)
    r=ResourcePool(w;potentials=(4,4,0),initial=(.1,.2,.3),diffusion=(.02,.02,.01))
    c=CellularEcology(w,reshape([1.,1.,0.],3,1),r;capacity=2,hidden=2,
        uptake_rate=0,basal_rate=0,neural_cost=0,growth_rate=0,catabolic_rate=0,
        speed=0,repulsion=0,stiffness=0,mutation=0)
    p=Interior.DeviceCellular(w,c;mechanical_substeps=1)
    w,p
end

@testset "device generic process protocol" begin
    w,p=device_protocol_fixture()
    ports=process_ports(p)
    @test isempty(ports.reads) && isempty(ports.writes) && isempty(ports.publishes)
    @test ports.owned==(:state,:environment)
    @test stability_limit(w,p)==Interior.device_cellular_timestep_limit(p)
    initial=Interior.owned_material_mass(p)
    @test initial==Interior.device_cellular_audit(Interior.device_cellular_snapshot(p)).material
    e=Experiment(w,(cells=p,);dt=.01)
    report=run!(e;steps=2,capture_every=2)
    @test w.time==2Float64(.01f0)
    @test p.environment.metadata[6]==2
    @test all(s->isapprox(s.diagnostics.owned_mass,initial;atol=1e-5),report.samples)
    @test all(s->s.diagnostics.finite,report.samples)
    @test report.metadata.ports[1].owned==(:state,:environment)
end

@testset "device generic preflight precedes sibling mutation" begin
    w,p=device_protocol_fixture()
    fill!(w.fields.nutrient,1);fill!(w.fields.activity,1)
    before=deepcopy(w.fields);state=Interior.device_cellular_snapshot(p)
    @test_throws ArgumentError step!(w,(Deposition(threshold=0),p),2stability_limit(w,p))
    @test w.fields==before && w.time==0
    @test Interior.device_cellular_snapshot(p)==state
    p.environment.status[1]=8
    @test_throws ArgumentError step!(w,(Deposition(threshold=0),p),.01)
    @test_throws ArgumentError Experiment(w,(first=Deposition(),cells=p);dt=.01)
    @test_throws ArgumentError Interior.owned_material_mass(p)
    @test_throws ArgumentError record_process_fields(p)
    @test w.fields==before && w.time==0
end

@testset "device recording flattens detached array roles" begin
    w,p=device_protocol_fixture();before=Interior.device_cellular_snapshot(p)
    saved=record_process_fields(p)
    @test saved.schema=="interior-device-cellular-observation"
    @test saved.parameter_names==keys(p.parameters)
    @test saved.geometry_names==keys(p.geometry)
    @test all(name->hasproperty(saved,Symbol("state_",name)),keys(p.state))
    @test all(name->hasproperty(saved,Symbol("environment_",name)),keys(p.environment))
    @test saved.state_mass!==p.state.mass
    @test Interior.device_cellular_snapshot(p)==before
    saved.state_mass[1]+=1
    @test p.state.mass[1]==before.state.mass[1]
    mktempdir() do dir
        manifest=record_state(w,p,joinpath(dir,"observation");step=0)
        state=only(manifest["processes"])["state"]
        @test haskey(state,"state_weights") && haskey(state,"environment_fields")
        @test state["state_ids"]["dtype"]=="Int32"
        @test state["state_mass"]["shape"]==[2]
        @test isfile(joinpath(dir,"observation",state["state_mass"]["file"]))
    end
    @test Interior.device_cellular_snapshot(p)==before
end

struct DeviceProtocolOwner{P}
    child::P
end
Interior.process_children(p::DeviceProtocolOwner)=(p.child,)

@testset "device Arena adapter remains explicitly unsupported" begin
    w,p=device_protocol_fixture()
    @test_throws ArgumentError Arena(w,(p,);dt=.01)
    @test_throws ArgumentError Arena(w,(DeviceProtocolOwner(p),);dt=.01)
    # An unchecked positional constructor must not bypass the checkpoint guard.
    raw=Arena(w,(p,),.01,0,false,false,0.,"",Dict{String,Any}())
    mktempdir() do dir
        path=joinpath(dir,"raw-device.jls")
        @test_throws ArgumentError checkpoint!(raw,path)
        @test !ispath(path)
    end
end


@testset "device compact World diagnostics inspect hidden state" begin
    w,p=device_protocol_fixture(;compact=true)
    @test isempty(w.fields) && isempty(w.scratch)
    report=run!(Experiment(w,(cells=p,);dt=.01);steps=1)
    d=last(report.samples).diagnostics
    @test d.grid_mass==0 && d.owned_mass>0
    @test d.activity_max===nothing && d.nutrient_min===nothing
    @test d.finite && w.time==Float64(.01f0)
    p.state.weights[1,1,1]=NaN32
    @test !Interior.process_state_finite(p)
    @test !experiment_diagnostics(w,(p,)).finite
end
