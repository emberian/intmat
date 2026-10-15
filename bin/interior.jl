#!/usr/bin/env julia
# Supported application entry point. Core physics remains in Interior.
const INTERIOR_ROOT = normpath(joinpath(@__DIR__, ".."))
const INTERIOR_STUDIES = ("acoustic", "writing", "electrical", "phase", "ecology", "mechanics", "grin", "compositions", "scout", "verify", "collection_phase", "collection_writing", "collection_waves", "collection_branches", "collection_mechanics", "collection_ecology", "alife_adaptation")
function interior_help(command=nothing)
    if command===nothing
        println("""
Interior — material simulations and living arenas
Usage: julia --project bin/interior.jl COMMAND [ARGS...]

  arena   [--profile TOML] [--output DIR] [--steps N | --continuous]
          [--resume CHECKPOINT] [--allow-source-change]
  run     [PROFILE [OUTPUT [cpu|metal]]] [--steps N] [--no-render]
  bigrun  PROFILE OUTPUT [--backend cpu|metal] [--steps N] [--dry-run]
  study   NAME [study-specific arguments]

Recipes: experiments/arena, experiments/materials, experiments/bigruns
Studies: $(join(INTERIOR_STUDIES, ", "))
Use COMMAND --help for command details. Outputs default under the repository's outputs/.
""")
    elseif command=="arena"
        println("arena [--profile TOML] [--output DIR] [--steps N | --continuous] [--resume CHECKPOINT] [--allow-source-change]\nDefault recipe: experiments/arena/active-2d.toml. Fresh output must be empty. Resuming requires an explicit checkpoint. Source migration requires --allow-source-change.")
    elseif command=="run"
        println("run [PROFILE [OUTPUT [cpu|metal]]] [--steps N] [--no-render]\nDefault recipe: experiments/materials/choir.toml. Requires a fresh or empty output directory. --steps overrides recipe duration; --no-render skips optical export and rendering.")
    elseif command=="bigrun"
        println("bigrun PROFILE OUTPUT [--backend cpu|metal] [--steps N] [--dry-run]\nRequires a fresh output directory. --dry-run estimates stability and memory with a small grid at the requested spacing.")
    elseif command=="study"
        println("study NAME [study-specific arguments]\nAvailable: $(join(INTERIOR_STUDIES, ", ")). Programs live in studies/; inspect their parameters before launching an expensive study.")
    else
        error("unknown command $command")
    end
end
function interior_main()
    isempty(ARGS) && return interior_help()
    command=popfirst!(ARGS)
    command in ("--help","-h","help") && return interior_help(isempty(ARGS) ? nothing : first(ARGS))
    command in ("arena","run","bigrun","study") || error("unknown command $command; use --help")
    any(arg->arg in ("--help","-h"),ARGS) && return interior_help(command)
    if command=="study"
        isempty(ARGS) && return interior_help(command)
        name=popfirst!(ARGS)
        name in INTERIOR_STUDIES || error("unknown study $name; use study --help")
        global INTERIOR_RUNNING_STUDY=name
        include(joinpath(INTERIOR_ROOT,"studies",name*".jl"))
    else
        include(joinpath(@__DIR__,"cli",command*".jl"))
    end
end
interior_main()
