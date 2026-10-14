using Interior, Test, Random
include("core.jl")
include("particles.jl")
include("render.jl")
isfile(joinpath(@__DIR__,"scenes.jl")) && include("scenes.jl")

include("neural.jl")
include("audit.jl")
include("clock.jl")

include("recording.jl")
include("acoustics.jl")
include("writing.jl")
include("electrical.jl")
include("phase.jl")
include("mechanics.jl")
include("ecology.jl")
include("experiments.jl")
include("process_composition.jl")

include("grin.jl")

include("implicit_phase.jl")
include("compartments.jl")
include("neural_ecology.jl")
include("arena.jl")
include("arena_composition.jl")

include("neural_batch.jl")
include("neural_workspace.jl")
include("resources.jl")
include("cooperative_boundaries.jl")
include("metabolic_ecology.jl")
include("metabolic_composition.jl")
