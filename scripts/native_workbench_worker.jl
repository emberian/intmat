#!/usr/bin/env julia
# A finite command adapter around the native Julia composition. No network/eval.
include("native_workbench_startup.jl")
const NATIVE_WORKBENCH_STARTUP_TOKEN = if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("PROFILE FRESH_OUTPUT required")
    NativeWorkbenchStartup.begin_startup(ARGS[2],@__FILE__)
else
    nothing
end
module NativeWorkbench
using ..NativeWorkbenchStartup
using Interior, Random, TOML, Serialization, SHA
using Interior: JSON3
include("native_constructed_garden.jl")
include("device_habitat_capture.jl")
# Pin the code loaded by this resident process, never whatever a concurrent editor
# happens to leave on disk at the next observation/acknowledgement.
const LOADED_CODE_PINS=(
    driver=bytes2hex(sha256(read(@__FILE__))),
    factory=bytes2hex(sha256(read(joinpath(@__DIR__,"native_constructed_garden.jl")))),
    capture=bytes2hex(sha256(read(joinpath(@__DIR__,"device_habitat_capture.jl")))),
    startup=bytes2hex(sha256(read(joinpath(@__DIR__,"native_workbench_startup.jl")))))
const FIELD_NAMES=(:free_a,:free_b,:raw_a,:raw_b,:enzyme_a,:enzyme_b,:bound_a,:bound_b,:permeability)
const ACTIONS=("describe","observe","step","pause","checkpoint","receptor_scales","fork","stop")
Base.@kwdef struct Limits
    max_steps::Int=4096
    max_steps_per_command::Int=256
    max_model_time::Float64=64.
    wall_seconds::Float64=900.
    idle_seconds::Float64=120.
    artifact_bytes::Int=512_000_000
    max_commands::Int=128
    max_branches::Int=4
end
mutable struct Branch{W,P,B,A}
    world::W
    process::P
    buffer::B
    initial::A
    failed::Bool
    observation::Any
    checkpoint::Any
    history::Vector{Any}
    interventions::Vector{Any}
    parent::Any
end
mutable struct Session
    root::String
    branches::Dict{String,Any}
    config::Dict{String,Any}
    limits::Limits
    dt::Float32
    side::Int
    started::Float64
    last_command::Float64
    sequence::Int
    cache::Dict{String,Any}
    stop_reason::Union{Nothing,String}
    busy::Bool
end
hashfile(path)=bytes2hex(sha256(read(path)))
function atomic_json(path,value)
    temp=path*".pending"
    open(io->JSON3.write(io,value),temp,"w");mv(temp,path;force=true)
end
diskbytes(s)=sum(filesize(joinpath(root,name)) for (root,_,files) in walkdir(s.root) for name in files;init=0)
function reserve!(s,bytes)
    diskbytes(s)+bytes<=s.limits.artifact_bytes || throw(ArgumentError("artifact budget reserve exhausted"))
end
function branches_state(s)
    Dict(name=>(step=b.process.cells.attempted_steps,time=b.world.time,paused=true,failed=b.failed,
        latest_observation=b.observation,latest_checkpoint=b.checkpoint,
        material_parameters=b.process.habitat.parameters,interventions=b.interventions,parent=b.parent) for (name,b) in s.branches)
end
function manifest!(s)
    atomic_json(joinpath(s.root,"worker.json"),(schema="interior-native-workbench-1",status=s.stop_reason===nothing ? (s.busy ? "running" : "ready") : "stopped",
        stop_reason=s.stop_reason,actions=ACTIONS,limits=NamedTuple{fieldnames(Limits)}(Tuple(getfield(s.limits,k) for k in fieldnames(Limits))),dt=s.dt,side=s.side,
        source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
        driver_sha256=LOADED_CODE_PINS.driver,factory_sha256=LOADED_CODE_PINS.factory,capture_helper_sha256=LOADED_CODE_PINS.capture,startup_helper_sha256=LOADED_CODE_PINS.startup,config=s.config,branches=branches_state(s),
        commands_completed=s.sequence,started_unix=s.started,wall_deadline_unix=s.started+s.limits.wall_seconds,wall_elapsed=time()-s.started,
        note="Paused between finite commands; observations are block means, never restart state."))
end
function accounting(b)
    a=b.initial;v=Interior.device_habitat_audit(Interior.device_habitat_snapshot(b.process))
    material=v.material-a.material-v.supplied+a.supplied+v.removed-a.removed
    potential=v.potential-a.potential+v.work-a.work+v.heat-a.heat-v.supplied_potential+a.supplied_potential+v.removed_potential-a.removed_potential
    tolerance=max(1e-8,1e-6*max(1.,a.material+v.supplied-a.supplied))
    if !(v.finite && abs(material)<=tolerance && abs(potential)<=4tolerance)
        b.failed=true;error("completed-state accounting guard")
    end
    (audit=v,material_error=material,potential_error=potential)
end
function observe!(s,name;event=nothing)
    b=s.branches[name];b.failed && error("branch has a partial failed state")
    p=b.process;c=p.cells;h=p.habitat;Interior.device_cellular_status!(c)
    side=s.side;reserve!(s,2side^3*9*4+Int(c.geometry.capacity)*1024+65536)
    habitat_coarse!(c.backend,128)(b.buffer,c.environment,h.arrays,h.feedstock.arrays,c.geometry,Int32(side);ndrange=side^3)
    fields=reshape(Array(b.buffer),side,side,side,9)
    active=findall(==(1),Array(c.state.active));select(a)=Array(a)[active,:]
    values=(fields=fields,ids=Array(c.state.ids)[active],parents=Array(c.state.parents)[active],lineages=Array(c.state.lineages)[active],
        generations=Array(c.state.generations)[active],position=select(c.state.position),radii=Array(c.state.radii)[active],
        mass=Array(c.state.mass)[active],energy=Array(c.state.energy)[active],metabolites=select(c.state.metabolites),
        bond_ids=select(c.state.bond_ids),hidden=select(c.state.hidden),adaptation=select(c.state.adaptation),actions=select(c.state.actions),
        material_inputs=select(c.state.inputs)[:,2Int(c.geometry.hidden)+16:2Int(c.geometry.hidden)+18],
        catalytic_expression=select(h.arrays.expression),last_step_secretion=select(h.arrays.last_secretion))
    all(a->all(isfinite,a),values) || error("nonfinite observation")
    folder=joinpath(s.root,"branches",name,"frames",lpad(s.sequence,5,'0')*(event===nothing ? "" : "-"*event.phase));mkdir(folder)
    arrays=Dict{String,Any}();offset=0;binary=joinpath(folder,"state.bin")
    open(binary*".pending","w") do io
        for (key,value) in pairs(values)
            data=key in (:ids,:parents,:lineages,:generations,:bond_ids) ? Int32.(value) : key===:energy ? Float64.(value) : Float32.(value)
            arrays[string(key)]=(offset=offset,bytes=sizeof(data),shape=size(data),type=eltype(data)===Int32 ? "i32" : eltype(data)===Float64 ? "f64" : "f32")
            write(io,data);offset+=sizeof(data)
        end
    end
    mv(binary*".pending",binary)
    audit=accounting(b)
    meta=(schema="native-garden-observatory-frame-1",branch=name,time=b.world.time,step=c.attempted_steps,
        shape=(side,side,side),native_shape=c.geometry.dims,extent=c.geometry.extent,population=length(active),field_names=FIELD_NAMES,
        field_units="concentration except dimensionless permeability",field_sampling=side==c.geometry.dims[1] ? "native voxel observations" : "explicit concentration block means",
        source_sha256=hashfile(binary),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
        binary="state.bin",binary_sha256=hashfile(binary),order="column-major",endian=ENDIAN_BOM==0x04030201 ? "little" : "big",arrays=arrays,
        secretion_note="Actual final completed step secretion, not an interval total.",accounting=audit,
        material_parameters=p.habitat.parameters,intervention=event,
        material_input_note="Receptor inputs are from the last completed step; new scales first act on the next physical step.")
    atomic_json(joinpath(folder,"frame.json"),meta)
    b.observation=(file=relpath(joinpath(folder,"frame.json"),s.root),sha256=hashfile(joinpath(folder,"frame.json")),
        step=c.attempted_steps,time=b.world.time,population=length(active),binary_bytes=filesize(binary),
        source_sha256=hashfile(joinpath(folder,"frame.json")),side=side,native_shape=c.geometry.dims,
        observation_event=event===nothing ? "" : event.id*"-"*event.phase)
    push!(b.history,merge(b.observation,(file=relpath(joinpath(folder,"frame.json"),joinpath(s.root,"branches",name)),)))
    observation_index!(s,name)
end
function observation_index!(s,name)
    b=s.branches[name]
    atomic_json(joinpath(s.root,"branches",name,"observatory.json"),
        (schema="native-garden-observatory-1",run_complete=s.stop_reason!==nothing,
         native_shape=b.process.cells.geometry.dims,physical_dt=s.dt,frames=b.history,interventions=b.interventions,
         observer=(stored_bytes=sum(x.binary_bytes for x in b.history),status="Paused between finite commands",stop_reason=s.stop_reason)))
end
function checkpoint!(s,name;event=nothing)
    b=s.branches[name];b.failed && error("branch has a partial failed state")
    reserve!(s,2resident_bytes(b.process)+2_000_000)
    ledger=accounting(b);payload=merge(Interior.device_habitat_checkpoint(b.world,b.process),
        (workbench=(branch=name,intervention=event,interventions=copy(b.interventions),parent=b.parent,accounting_baseline=b.initial),))
    path=joinpath(s.root,"checkpoints",name*"-"*lpad(s.sequence,5,'0')*(event===nothing ? "" : "-"*event.phase)*".jls")
    serialize(path*".pending",payload);mv(path*".pending",path)
    b.checkpoint=(file=relpath(path,s.root),sha256=hashfile(path),bytes=filesize(path),step=b.process.cells.attempted_steps,time=b.world.time,accounting=ledger)
end
function create_session(profile,root;backend=CPU(),side=8,limits=Limits(),startup_token=nothing)
    started=ispath(root) ? NativeWorkbenchStartup.validate(root,startup_token) : time()
    all(x->isfinite(x)&&x>0,(limits.max_steps,limits.max_steps_per_command,limits.max_model_time,limits.wall_seconds,limits.idle_seconds,limits.artifact_bytes,limits.max_commands)) || error("finite positive limits required")
    2<=limits.max_branches<=4 || error("branch cap must be between two and four")
    limits.max_steps_per_command<=limits.max_steps || error("command cap exceeds session cap")
    cfg=profile isa AbstractString ? TOML.parsefile(profile) : deepcopy(profile)
    haskey(cfg,"material_sensing") || error("profile must declare the sensing branch receptors")
    shape=cfg["grid"]["shape"];0<side<=minimum(shape)&&all(n->n%side==0,shape) || error("observation side must divide native shape")
    dt=Float32(cfg["native"]["dt"]);bs=Dict{String,Any}()
    for name in ("off","sensing")
        trial=deepcopy(cfg);name=="off"&&delete!(trial,"material_sensing")
        w,p,funding=native_constructed_garden(backend,trial)
        buffer=Interior.KernelAbstractions.zeros(backend,Float32,side^3,9)
        bs[name]=Branch(w,p,buffer,funding.initial,false,nothing,nothing,Any[],Any[],nothing)
    end
    # Same physiological state; receptor semantics are the sole branch-law edit.
    a=Interior.device_habitat_snapshot(bs["off"].process);b=Interior.device_habitat_snapshot(bs["sensing"].process)
    a.cells==b.cells && a.arrays==b.arrays && a.feedstock==b.feedstock || error("fresh paired state mismatch")
    root=abspath(root);mkpath(root);mkdir(joinpath(root,"checkpoints"))
    for name in ("off","sensing");mkpath(joinpath(root,"branches",name,"frames"));end
    s=Session(root,bs,cfg,limits,dt,side,started,time(),0,Dict{String,Any}(),nothing,false)
    reserve!(s,sum(2resident_bytes(b.process)+2_000_000 for b in values(bs))+2_000_000)
    for name in ("off","sensing");checkpoint!(s,name);observe!(s,name);end
    s.last_command=time();manifest!(s);s
end
function validate_command(c,s)
    c isa AbstractDict || throw(ArgumentError("command must be an object"))
    Set(keys(c))<=Set(("id","action","branch","steps","scales","requested_at")) || throw(ArgumentError("unknown command fields"))
    id=get(c,"id",nothing);id isa String && occursin(r"^[A-Za-z0-9_-]{1,80}$",id) || throw(ArgumentError("invalid command id"))
    action=get(c,"action",nothing);action in ACTIONS || throw(ArgumentError("unknown action"))
    branch=get(c,"branch","both");(branch=="both" || haskey(s.branches,branch)) || throw(ArgumentError("unknown branch"))
    if action=="step"
        n=get(c,"steps",nothing);n isa Integer && !(n isa Bool) && 0<n<=s.limits.max_steps_per_command || throw(ArgumentError("step count exceeds finite command cap"))
    else
        haskey(c,"steps") && throw(ArgumentError("steps only applies to step"))
    end
    if action=="receptor_scales"
        branch!="both" && hasproperty(s.branches[branch].process.habitat.parameters,:material_scales) || throw(ArgumentError("selected branch has no material receptors"))
        receptor_scales(get(c,"scales",nothing))
    else
        haskey(c,"scales") && throw(ArgumentError("scales only applies to receptor_scales"))
    end
    action=="fork" && branch=="both" && throw(ArgumentError("fork requires one complete branch"))
    id,action,branch=="both" ? ["off","sensing"] : [branch]
end
function receptor_scales(values)
    values isa AbstractVector && length(values)==3 || throw(ArgumentError("exactly three receptor scales required"))
    all(x->x isa Real && !(x isa Bool) && isfinite(x) && x>0,values) || throw(ArgumentError("positive finite numeric receptor scales required"))
    scales=Tuple(Float32.(values))
    all(x->isfinite(x)&&x>0,scales) || throw(ArgumentError("receptor scales must remain positive and finite in Float32"))
    scales
end
function receptor_process(b,scales)
    p=b.process;h=p.habitat;c=p.cells
    hasproperty(h.parameters,:material_scales) || error("branch lacks material receptors")
    q=merge(h.parameters,(material_scales=scales,))
    Interior.dh_validate_parameters(q,h.feedstock.parameters)
    argument=Interior.device_habitat_argument(c,h.arrays,q)
    ownership=merge(h.ownership,(parameters=q,argument=argument))
    habitat=typeof(h)(c,h.arrays,h.feedstock,q,argument,ownership)
    candidate=typeof(p)(c,habitat)
    Interior.validate(b.world,candidate)
    candidate
end
function intervene_receptors!(s,name,id,scales)
    b=s.branches[name];b.failed&&error("branch has a partial failed state")
    candidate=receptor_process(b,scales)
    candidate.habitat.parameters.material_channels==b.process.habitat.parameters.material_channels || error("fixed receptor channels changed")
    reserve!(s,4resident_bytes(b.process)+4_000_000+4s.side^3*9*4)
    event=(id=id,action="receptor_scales",sequence=s.sequence,time=b.world.time,step=b.process.cells.attempted_steps,
        channels=candidate.habitat.parameters.material_channels,before=b.process.habitat.parameters.material_scales,after=scales)
    before=merge(event,(phase="before",));after=merge(event,(phase="after",))
    checkpoint!(s,name;event=before);observe!(s,name;event=before)
    before_checkpoint=b.checkpoint;before_observation=b.observation
    # Single host-wrapper swap; physical arrays, cells, and controller state are shared unchanged.
    b.process=candidate
    push!(b.interventions,merge(event,(before_checkpoint=before_checkpoint,before_observation=before_observation)))
    manifest!(s)
    checkpoint!(s,name;event=after);observe!(s,name;event=after)
end
function fork_current!(s,name)
    length(s.branches)<s.limits.max_branches || error("branch budget exhausted")
    b=s.branches[name];b.failed&&error("cannot fork a partial failed branch")
    child=first(filter(n->!haskey(s.branches,n),("fork1","fork2")))
    reserve!(s,4resident_bytes(b.process)+4_000_000+2s.side^3*9*4)
    checkpoint!(s,name)
    parent_checkpoint=b.checkpoint;path=joinpath(s.root,parent_checkpoint.file)
    hashfile(path)==parent_checkpoint.sha256 || error("parent checkpoint hash mismatch")
    payload=deserialize(path)
    w=World(b.world.grid,b.world.backend,NamedTuple(),NamedTuple(),b.world.time,b.world.seed)
    process=Interior.device_habitat_restore(w,payload)
    Interior.device_habitat_snapshot(process)==Interior.device_habitat_snapshot(b.process) || error("strict fork differs from parent state")
    provenance=(branch=name,checkpoint=parent_checkpoint,step=b.process.cells.attempted_steps,time=b.world.time,
        source=payload.source,cellular_source=payload.cells.source,parameters=payload.parameters,sequence=s.sequence)
    buffer=Interior.KernelAbstractions.zeros(w.backend,Float32,s.side^3,9)
    fresh=Branch(w,process,buffer,b.initial,false,nothing,nothing,Any[],deepcopy(b.interventions),provenance)
    mkpath(joinpath(s.root,"branches",child,"frames"));s.branches[child]=fresh
    checkpoint!(s,child);observe!(s,child)
    child
end
function advance_one!(b,dt)
    try
        Interior.device_cellular_step!(b.world,b.process.cells,dt;anabolism=b.process.habitat)
        b.world.time=Float64(b.process.cells.attempted_steps)*Float64(dt)
    catch
        b.failed=true;rethrow()
    end
end
function process_command!(s,c)
    # Bad commands get an acknowledgement but never a physics step. IDs are
    # idempotent within this bounded worker; altered reuse is explicitly refused.
    id=string(get(c,"id","invalid"));digest=bytes2hex(sha256(JSON3.write(c)))
    if haskey(s.cache,id)
        saved=s.cache[id]
        ack=saved.digest==digest ? saved.ack : (id=id,ok=false,accepted=false,status="rejected",message="command id reused with changed content")
        atomic_json(joinpath(s.root,"control.json"),ack);return ack
    end
    ok=true;accepted=false;message="completed";s.sequence+=1
    try
        s.stop_reason===nothing || error("worker is stopped")
        s.sequence<=s.limits.max_commands || error("command budget exhausted")
        time()-s.started<s.limits.wall_seconds || error("wall deadline reached")
        id,action,names=validate_command(c,s)
        reserve!(s,65536+sum(Int(b.process.cells.geometry.capacity)*1024+2s.side^3*9*4 for b in values(s.branches)))
        if action=="step"
            n=Int(c["steps"])
            for name in names
                b=s.branches[name];b.failed&&error("branch has a partial failed state")
                b.process.cells.attempted_steps+n<=s.limits.max_steps || error("branch step budget exceeded")
                b.world.time+n*Float64(s.dt)<=s.limits.max_model_time || error("branch model-time budget exceeded")
            end
            accepted=true;s.busy=true;manifest!(s)
            for _ in 1:n, name in names
                time()-s.started<s.limits.wall_seconds || error("wall deadline reached at complete step boundary")
                advance_one!(s.branches[name],s.dt)
            end
            for name in names;observe!(s,name);end
        elseif action=="fork"
            length(s.branches)<s.limits.max_branches || error("branch budget exhausted")
            s.branches[only(names)].failed && error("cannot fork a partial failed branch")
            accepted=true;s.busy=true;manifest!(s)
            fork_current!(s,only(names))
        elseif action=="receptor_scales"
            accepted=true;s.busy=true;manifest!(s)
            intervene_receptors!(s,only(names),id,receptor_scales(c["scales"]))
        elseif action=="observe"
            accepted=true;
            for name in names;observe!(s,name);end
        elseif action=="checkpoint"
            accepted=true;
            for name in names;checkpoint!(s,name);end
        elseif action=="stop"
            accepted=true;
            s.stop_reason="stop command"
        else
            accepted=true
        end
    catch e
        ok=false;message=sprint(showerror,e)
    end
    s.busy=false;s.last_command=time()
    if s.sequence>=s.limits.max_commands;s.stop_reason="command budget";end
    if time()-s.started>=s.limits.wall_seconds;s.stop_reason="wall deadline";end
    ack=(schema="interior-native-workbench-ack-1",id=id,ok=ok,accepted=accepted,status=ok ? "completed" : accepted ? "failed" : "rejected",message=message,sequence=s.sequence,branches=branches_state(s))
    s.cache[id]=(digest=digest,ack=ack)
    open(joinpath(s.root,"commands.jsonl"),"a") do io;JSON3.write(io,(command=c,ack=ack));write(io,'\n');end
    atomic_json(joinpath(s.root,"control.json"),ack);manifest!(s);ack
end
function serve!(s;poll_seconds=.1)
    last_hash=""
    while s.stop_reason===nothing
        if time()-s.started>=s.limits.wall_seconds;s.stop_reason="wall deadline";break;end
        if time()-s.last_command>=s.limits.idle_seconds;s.stop_reason="idle deadline";break;end
        path=joinpath(s.root,"command.json")
        if islink(path);s.stop_reason="command path must be an owned regular file";break;end
        if isfile(path)
            filesize(path)<=4096 || (s.stop_reason="command file exceeds byte cap";break)
            raw=read(path,String);digest=bytes2hex(sha256(raw))
            if digest!=last_hash
                last_hash=digest
                try
                    c=JSON3.read(raw,Dict{String,Any});process_command!(s,c)
                catch e
                    atomic_json(joinpath(s.root,"control.json"),(id="invalid",ok=false,message=sprint(showerror,e)))
                end
            end
        end
        sleep(poll_seconds)
    end
    for name in keys(s.branches);observation_index!(s,name);end
    manifest!(s);s
end
end # module
if abspath(PROGRAM_FILE)==@__FILE__
    try
        include("device_cellular_setup.jl")
        session=NativeWorkbench.create_session(ARGS[1],ARGS[2];backend,
            side=parse(Int,get(ENV,"WORKBENCH_SIDE","8")),startup_token=NATIVE_WORKBENCH_STARTUP_TOKEN,
            limits=NativeWorkbench.Limits(wall_seconds=parse(Float64,get(ENV,"WORKBENCH_WALL_SECONDS","900")),
                idle_seconds=parse(Float64,get(ENV,"WORKBENCH_IDLE_SECONDS","120"))))
        NativeWorkbench.serve!(session)
    catch e
        NativeWorkbenchStartup.failed(ARGS[2],NATIVE_WORKBENCH_STARTUP_TOKEN,e)
        rethrow()
    end
end
