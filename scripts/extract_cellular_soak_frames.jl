#!/usr/bin/env julia
# Read trusted archived checkpoints with the original source; never step or save.
using Interior, JSON3, SHA
root=ARGS[1];result=Any[]
for name in ("broad","mosaic")
    dir=joinpath(root,name)
    files=sort(filter(f->startswith(f,"state-")&&endswith(f,".jls"),readdir(dir)))
    selected=[("initial",first(files)),("midpoint",files[cld(length(files),2)]),
              ("final",String(JSON3.read(read(joinpath(dir,"checkpoint.json"),String)).file))]
    for (label,file) in selected
        path=joinpath(dir,file);a=restore_arena(path)
        snapshot=arena_snapshot(a;max_side=32,max_agents=3000)
        push!(result,Dict("run"=>name,"label"=>label,"file"=>file,"sha256"=>bytes2hex(sha256(read(path))),
            "step"=>snapshot["step"],"time"=>snapshot["time"],"grid"=>snapshot["grid"],
            "metrics"=>snapshot["metrics"],"agents"=>snapshot["agents"]))
    end
end
JSON3.write(stdout,result);println()
