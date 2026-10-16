#!/usr/bin/env julia
include("living_support_worklists.jl")
LivingSupportWorklists.profile(ARGS[1];resolution=length(ARGS)>1 ? parse(Int,ARGS[2]) : 64)
