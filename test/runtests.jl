using Interior, Test, Random
include("core.jl")
include("particles.jl")
include("render.jl")
isfile(joinpath(@__DIR__,"scenes.jl")) && include("scenes.jl")

include("neural.jl")
include("audit.jl")
