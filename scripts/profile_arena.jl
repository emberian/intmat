#!/usr/bin/env julia
using Interior, JSON3, Statistics
include(joinpath(@__DIR__, "arena_experiments.jl"))
isdefined(Interior, :neural_controller_outputs) || Base.include(Interior, joinpath(@__DIR__, "../src/neural_batch.jl"))
isdefined(Interior, :NeuralWorkspace) || Base.include(Interior, joinpath(@__DIR__, "../src/neural_workspace.jl"))
function measured(f; repeats=8)
    f()
    samples=[@timed f() for _ in 1:repeats]
    (;seconds=median(s.time for s in samples),bytes=median(s.bytes for s in samples))
end
function profile_population(n)
    a,phase,m,c,flux=experiment_arena(;dims=(64,64,1),capacity=n,founders=n,lifespan=10000f0)
    c.storage .= 0; c.energy .= 100
    dims=a.world.grid.size
    snapshot=copy(c.positions); cells=[Interior._neural_cell(view(snapshot,:,i),dims) for i in 1:n]
    bins=Dict{CartesianIndex{3},Vector{Int}}()
    for i in 1:n; push!(get!(bins,cells[i],Int[]),i); end
    inputs=[Interior._neural_input(a.world,c,i,cells[i],bins,snapshot) for i in 1:n]
    sensors=measured(()->[Interior._neural_input(a.world,c,i,cells[i],bins,snapshot) for i in 1:n])
    controller=measured(()->[tanh.((c.weights[:,:,i].+c.plastic_marks[:,:,i])*inputs[i]) for i in 1:n])
    fused_controller=measured(()->Interior.neural_controller_outputs(c.weights,c.plastic_marks,inputs))
    workspace=Interior.NeuralWorkspace(a.world,c)
    frame=measured(()->Interior.workspace_prepare!(a.world,c,workspace))
    process=Dict(string(typeof(p))=>measured(()->advance!(a.world,p,a.dt)) for p in a.processes)
    complete=measured(()->arena_step!(a))
    snapshot_cost=measured(()->arena_snapshot(a))
    (;population=n,sensors,controller,fused_controller,frame,process,complete,snapshot=snapshot_cost)
end
populations=isempty(get(ENV,"ARENA_PROFILE_POPULATIONS","")) ? (256,1000,2400) : parse.(Int,split(ENV["ARENA_PROFILE_POPULATIONS"],","))
results=[profile_population(n) for n in populations]
println(JSON3.write(results))
if !isempty(ARGS); open(io->JSON3.write(io,results),ARGS[1],"w"); end
