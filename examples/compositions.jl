if abspath(PROGRAM_FILE)==@__FILE__; global INTERIOR_RUNNING_STUDY="compositions"; end
# Compatibility entry point; canonical interface: bin/interior.jl
include(joinpath(@__DIR__, "..", "studies", "compositions.jl"))
