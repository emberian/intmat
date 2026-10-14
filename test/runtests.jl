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
