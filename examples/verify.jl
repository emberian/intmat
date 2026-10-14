if abspath(PROGRAM_FILE)==@__FILE__; global INTERIOR_RUNNING_STUDY="verify"; end
# Compatibility entry point; canonical interface: bin/interior.jl
include(joinpath(@__DIR__, "..", "studies", "verify.jl"))
