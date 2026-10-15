#!/usr/bin/env julia
# Long-running CPU experiment driver. All observations/checkpoints occur only
# between complete physics steps. Source migration requires an explicit fresh fork.
using Interior, TOML, JSON3, Serialization, SHA

function soak_main(args)
    length(args)>=2 || error("usage: cellular_soak.jl PROFILE OUTPUT [--seconds N] [--resume CHECKPOINT] [--budget-bytes N] [--metrics-seconds N] [--archive-seconds N] [--throughput] [--allow-source-change]")
    profile,output=abspath(args[1]),abspath(args[2]); options=Dict{String,String}();flags=Set{String}()
    i=3
    while i<=length(args)
        if args[i] in ("--throughput","--allow-source-change")
            push!(flags,args[i]);i+=1;continue
        end
        args[i] in ("--seconds","--resume","--budget-bytes","--metrics-seconds","--archive-seconds") || error("unknown argument $(args[i])")
        i<length(args) || error("missing argument value")
        options[args[i]]=args[i+1];i+=2
    end
    seconds=parse(Float64,get(options,"--seconds","43200"))
    budget=parse(Int,get(options,"--budget-bytes","4000000000"))
    metric_seconds=parse(Float64,get(options,"--metrics-seconds","20"))
    archive_seconds=parse(Float64,get(options,"--archive-seconds","600"))
    all(x->isfinite(x)&&x>0,(seconds,metric_seconds,archive_seconds)) && budget>=1_000_000 || error("invalid budget/cadence")
    config=TOML.parsefile(profile)
    get(config,"model","")=="cellular" || error("cellular recipe required")
    resume=get(options,"--resume",nothing)
    throughput="--throughput" in flags
    allow_source_change="--allow-source-change" in flags
    allow_source_change && resume===nothing && error("--allow-source-change requires --resume")
    if allow_source_change && ispath(output)
        isdir(output) && isempty(readdir(output)) || error("source migration requires a fresh empty output directory")
    end
    migration=nothing;parent_bytes=nothing
    arena=if resume===nothing
        build_arena(config)
    elseif allow_source_change
        # Capture one immutable parent image: rotating source slots can change
        # between reads. Both restoration and provenance use these exact bytes.
        parent_bytes=read(abspath(resume))
        temporary,io=mktemp()
        try
            write(io,parent_bytes);close(io)
            restored=restore_arena(temporary;allow_source_change=true)
            envelope=open(Serialization.deserialize,temporary)
            migration=Dict("parent_checkpoint"=>abspath(resume),"parent_checkpoint_sha256"=>bytes2hex(sha256(parent_bytes)),
                "parent_sources"=>envelope.sources,"loaded_sources"=>Interior.arena_sources(),"initial_step"=>restored.steps,
                "preserved_parent_file"=>"parent-checkpoint.jls")
            restored
        finally
            isopen(io) && close(io)
            isfile(temporary) && rm(temporary)
        end
    else
        restore_arena(abspath(resume))
    end
    if resume!==nothing
        saved=Dict(k=>v for (k,v) in arena.config if !startswith(k,"_"))
        saved==config || error("resume profile differs from the checkpoint configuration")
    end
    arena.stopped=false
    arena.paused && error("paused checkpoint requires an explicit reviewed resume before this soak")
    lock=Interior.arena_lock(output;resume=resume!==nothing && !allow_source_change)
    session=string(round(Int,time()*1e6),"-",getpid())
    reserve=min(128_000_000,budget÷8)
    bytes_used()=sum(filesize(joinpath(dir,file)) for (dir,_,files) in walkdir(output) for file in files;init=0)
    function room(path,n;final=false)
        old=isfile(path) ? filesize(path) : 0
        bytes_used()-old+n<=budget-(final ? 0 : reserve) || throw(OverflowError("artifact budget reached; existing evidence retained"))
    end
    function jsonwrite(name,value;final=false)
        data=Vector{UInt8}(codeunits(JSON3.write(value)))
        path=joinpath(output,name);room(path,length(data);final)
        Interior.arena_atomic(io->write(io,data),path)
    end
    function checkpointwrite(name;final=false)
        # Only this newly created temporary is discarded if budget admission fails.
        temporary,io=mktemp(dirname(output));close(io)
        try
            checkpoint!(arena,temporary)
            target=joinpath(output,name);room(target,filesize(temporary);final)
            mv(temporary,target;force=true)
        finally
            isfile(temporary) && rm(temporary)
        end
    end
    slot=0
    function restartsave(;final=false)
        slot=mod(slot,3)+1;name="checkpoint-$slot.jls"
        checkpointwrite(name;final)
        jsonwrite("checkpoint.json",Dict("file"=>name,"step"=>arena.steps,"session"=>session);final)
    end
    initial_step=arena.steps
    started=time();lastmetrics=-Inf;lastcheckpoint=started;lastarchive=-Inf
    rows=0;chunk=1;reason="wall_time_complete";consistent=true;last_completed_step=arena.steps
    function observe(;final=false)
        snapshot=arena_snapshot(arena)
        jsonwrite("latest.json",snapshot;final)
        row=Dict("session"=>session,"observed_at_unix"=>time(),"step"=>arena.steps,"time"=>arena.world.time,"metrics"=>snapshot["metrics"])
        line=Vector{UInt8}(codeunits(JSON3.write(row)*"\n"))
        path=joinpath(output,"metrics-$session-$(lpad(chunk,4,'0')).jsonl")
        # Append only to this session's numbered chunk; total budget includes it.
        bytes_used()+length(line)<=budget-(final ? 0 : reserve) || throw(OverflowError("artifact budget reached"))
        open(io->write(io,line),path,"a")
        rows+=1
        if rows>=256;chunk+=1;rows=0;end
        jsonwrite("status.json",Dict("session"=>session,"pid"=>getpid(),"state"=>final ? "stopped" : "running","reason"=>final ? reason : nothing,
            "step"=>arena.steps,"simulation_time"=>arena.world.time,"started_unix"=>started,"updated_unix"=>time(),"deadline_unix"=>started+seconds,
            "artifact_bytes"=>bytes_used(),"artifact_budget_bytes"=>budget,"julia_threads"=>Threads.nthreads(),"throughput"=>throughput,"initial_step"=>initial_step,"elapsed_wall_seconds"=>time()-started);final)
    end
    try
        if parent_bytes!==nothing
            parent=joinpath(output,"parent-checkpoint.jls");room(parent,length(parent_bytes))
            Interior.arena_atomic(io->write(io,parent_bytes),parent)
        end
        jsonwrite("session-$session.json",Dict("schema"=>"interior-cellular-soak","profile"=>profile,"profile_sha256"=>bytes2hex(sha256(read(profile))),
            "julia_version"=>string(VERSION),"loaded_sources"=>Interior.arena_sources(),"wrapper_sha256"=>bytes2hex(sha256(read(@__FILE__))),
            "started_unix"=>started,"wall_seconds"=>seconds,"artifact_budget_bytes"=>budget,"resumed_from"=>resume,"initial_step"=>arena.steps,
            "backend"=>"Julia CPU","arguments"=>args,"throughput"=>throughput,
            "batch_sleep_seconds"=>throughput ? 0.0 : get(config,"sleep_seconds",0.01),"source_migration"=>migration))
        restartsave()
        while time()-started<seconds && !arena.stopped
            if isfile(joinpath(output,"STOP"));reason="operator_stop_file";break;end
            now=time()
            if now-lastmetrics>=metric_seconds;observe();lastmetrics=now;end
            if now-lastarchive>=archive_seconds
                checkpointwrite("state-$(lpad(arena.steps,12,'0'))-$session.jls");lastarchive=now
            end
            if now-lastcheckpoint>=get(config,"checkpoint_seconds",90.0);restartsave();lastcheckpoint=now;end
            for _ in 1:get(config,"batch",2)
                consistent=false
                arena_step!(arena)
                last_completed_step=arena.steps
                consistent=true
            end
            throughput || sleep(get(config,"sleep_seconds",0.01))
        end
    catch error
        if error isa InterruptException;reason="interrupted"
        elseif error isa OverflowError;reason="artifact_budget_reached"
        else;reason="error: "*sprint(showerror,error);rethrow();end
    finally
        try
            if consistent
                restartsave(;final=true)
                observe(;final=true)
                println("Soak saved step ",arena.steps," reason=",reason)
            else
                # A stage may have mutated fields before throwing. Publish no
                # executable checkpoint or field observation of partial state.
                jsonwrite("status.json",Dict("session"=>session,"pid"=>getpid(),"state"=>"failed",
                    "reason"=>reason,"last_completed_step"=>last_completed_step,"unsafe_in_memory_state"=>true,
                    "checkpoint_preserved"=>true,"started_unix"=>started,"updated_unix"=>time(),
                    "deadline_unix"=>started+seconds,"artifact_bytes"=>bytes_used(),"artifact_budget_bytes"=>budget);final=true)
                println("Mid-step failure; prior checkpoint preserved at or before step ",last_completed_step)
            end
        finally
            Interior.arena_unlock(lock)
        end
    end
    arena
end
abspath(PROGRAM_FILE)==(@__FILE__) && soak_main(ARGS)
