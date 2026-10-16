#!/usr/bin/env julia
length(ARGS)==4||error("COHORT_PROFILE COHORT_OUTPUT LIVING_PROFILE LIVING_OUTPUT")
original_arguments=copy(ARGS)
empty!(ARGS);append!(ARGS,original_arguments[1:2])
module CohortGate
include("probe_cohort_construction.jl")
end
empty!(ARGS);append!(ARGS,original_arguments[3:4])
module LivingGate
include("probe_living_construction.jl")
end
