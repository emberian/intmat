using Interior, TOML, JSON3
function main(args)
    profile=joinpath(@__DIR__,"..","..","experiments","arena","active-2d.toml"); dir=joinpath(@__DIR__,"..","..","outputs","arena","active-2d"); resume=nothing; steps=100; continuous=false; allow_source_change=false
    i=1
    while i<=length(args)
        arg=args[i]
        if arg=="--continuous"; continuous=true
        elseif arg=="--allow-source-change"; allow_source_change=true
        elseif arg in ("--profile","--output","--resume","--steps")
            i+=1; i<=length(args) || error("missing value for $arg")
            if arg=="--profile"; profile=args[i]
            elseif arg=="--output"; dir=args[i]
            elseif arg=="--resume"; resume=args[i]
            else; steps=parse(Int,args[i]); end
        else; error("unknown argument $arg"); end
        i+=1
    end
    allow_source_change && resume===nothing && error("--allow-source-change requires --resume")
    config=TOML.parsefile(profile)
    arena=resume===nothing ? build_arena(config) : restore_arena(resume;allow_source_change)
    resume===nothing || (arena.stopped=false)
    println("Arena output: ",abspath(dir),"; stop with command.json or Ctrl-C")
    run_arena!(arena,dir;resume=resume!==nothing,steps=continuous ? nothing : steps,
        snapshot_seconds=get(config,"snapshot_seconds",1.0),checkpoint_seconds=get(config,"checkpoint_seconds",60.0),
        checkpoint_slots=get(config,"checkpoint_slots",3),batch=get(config,"batch",5),sleep_seconds=get(config,"sleep_seconds",0.002))
    println("Saved step ",arena.steps," at simulation time ",arena.world.time)
end
main(ARGS)
