#!/usr/bin/env julia
include("living_resource_fusion.jl")
LivingResourceFusion.benchmark(ARGS[1];resolution=length(ARGS)>1 ? parse(Int,ARGS[2]) : 64)
