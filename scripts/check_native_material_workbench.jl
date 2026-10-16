#!/usr/bin/env julia
# Bounded host adapter gate; all physical work uses the resident native law.
include("native_workbench_worker.jl")
include("device_cellular_setup.jl")
using Serialization
const M=NativeWorkbench
function gate(profile,root)
    s=M.create_session(profile,root;backend,living=true,side=4,limits=M.Limits(max_steps=64,max_steps_per_command=16,max_model_time=2.,wall_seconds=900.,idle_seconds=120.))
    acks=Any[]
    function command(id,action,branch;kwargs...)
        request=merge(Dict{String,Any}("id"=>id,"action"=>action,"branch"=>branch),Dict(string(k)=>v for (k,v) in kwargs))
        ack=M.process_command!(s,request);push!(acks,(id=id,action=action,branch=branch,ok=ack.ok,accepted=ack.accepted,message=ack.message))
        ack
    end
    try
    rates=Dict("decay_rate"=>2.,"desorption_rate"=>1.,"raw_diffusion"=>.03)
    @assert !command("no-original-edit","material_rates","sensing";rates).accepted
    @assert command("seed","step","sensing";steps=4).ok
    @assert command("fork","fork","sensing").ok
    p=s.branches["sensing"];c=s.branches["fork1"];original=Interior.device_habitat_snapshot(c.process)
    private=M.living_snapshot(c.policy);parent=Interior.device_habitat_snapshot(p.process)
    parent_checkpoint=p.checkpoint;parent_hash=M.hashfile(joinpath(s.root,parent_checkpoint.file))
    @assert !command("unsafe-diffusion","material_rates","fork1";rates=merge(rates,Dict("raw_diffusion"=>1.))).accepted
    @assert command("rates","material_rates","fork1";rates).ok
    changed=Interior.device_habitat_snapshot(c.process)
    @assert changed.cells==original.cells && changed.arrays==original.arrays && changed.feedstock==original.feedstock
    @assert M.living_snapshot(c.policy)==private && Interior.device_habitat_snapshot(p.process)==parent
    frames=c.history[end-1:end];metadata=[JSON3.read(read(joinpath(s.root,"branches/fork1",f.file),String)) for f in frames]
    @assert metadata[1].binary_sha256==metadata[2].binary_sha256
    cp=deserialize(joinpath(s.root,c.checkpoint.file));wr=World(c.world.grid,backend,NamedTuple(),NamedTuple(),c.world.time,c.world.seed)
    restored=Interior.device_habitat_restore(wr,cp);@assert Interior.device_habitat_snapshot(restored)==changed
    @assert command("child","step","fork1";steps=8).ok
    @assert Interior.device_habitat_snapshot(p.process)==parent
    @assert command("parent","step","sensing";steps=8).ok
    @assert c.world.time==p.world.time && c.observation.material_summary.catalyst!=p.observation.material_summary.catalyst
    @assert M.hashfile(joinpath(s.root,parent_checkpoint.file))==parent_hash
    @assert command("stop","stop","both").ok
    result=(schema="native-material-workbench-gate-1",runtime=device_cellular_runtime(),source=Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=Interior.DEVICE_HABITAT_SOURCE_PIN,
        driver_sha256=M.LOADED_CODE_PINS.driver,limits=NamedTuple{fieldnames(M.Limits)}(Tuple(getfield(s.limits,k) for k in fieldnames(M.Limits))),
        commands=acks,original_state_preserved=true,private_state_preserved=true,strict_restore=true,
        same_time_binary_sha256=metadata[2].binary_sha256,parent_checkpoint_sha256=parent_hash,
        intervention=last(c.interventions),left_summary=p.observation.material_summary,right_summary=c.observation.material_summary,
        clocks=Dict(n=>(step=b.process.cells.attempted_steps,time=b.world.time) for (n,b) in s.branches),stopped=s.stop_reason)
    M.atomic_json(joinpath(s.root,"gate.json"),result);println(JSON3.write((ok=true,clocks=result.clocks,left=result.left_summary,right=result.right_summary)))
    catch err
        s.stop_reason="gate failure: "*sprint(showerror,err);M.manifest!(s);rethrow()
    end
end
gate(ARGS[1],ARGS[2])
