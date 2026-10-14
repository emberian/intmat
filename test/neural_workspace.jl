using Test, Random
isdefined(Interior,:NeuralWorkspace) || Base.include(Interior,joinpath(@__DIR__,"../src/neural_workspace.jl"))
@testset "contiguous observation frame reference" begin
    for dims in ((64,64,1),(3,2,1),(5,4,3)), n in (0,1,53)
        w=World(Grid(dims...)); rng=MersenneTwister(66)
        pos=rand(rng,Float32,3,n).*reshape(Float32.(collect(dims)),3,1)
        dims[3]==1 && (pos[3,:].=0)
        p=NeuralEcology(pos;seed=66,capacity=max(n,1)); w.fields.nutrient.=rand(rng,Float32,dims); w.fields.deposit.=rand(rng,Float32,dims)
        snapshot=copy(p.positions); cells=[Interior._neural_cell(view(snapshot,:,i),dims) for i in 1:n]; bins=Dict{CartesianIndex{3},Vector{Int}}()
        for i in 1:n; push!(get!(bins,cells[i],Int[]),i); end
        inputs=[Interior._neural_input(w,p,i,cells[i],bins,snapshot) for i in 1:n]
        out=[tanh.((p.weights[:,:,i]+p.plastic_marks[:,:,i])*inputs[i]) for i in 1:n]
        ws=Interior.NeuralWorkspace(w,p); Interior.workspace_prepare!(w,p,ws)
        @test ws.cells==cells
        @test all(isapprox(ws.inputs[:,i],inputs[i];atol=1f-7,rtol=1f-6) for i in 1:n)
        @test all(isapprox(ws.outputs[:,i],out[i];atol=2f-6,rtol=2f-5) for i in 1:n)
        @test all(isfinite,ws.inputs) && all(isfinite,ws.outputs)
        original_inputs=ws.inputs; Interior.workspace_prepare!(w,p,ws); @test original_inputs===ws.inputs
    end
end
@testset "workspace CSR fair resource competition" begin
    w=World(Grid(3,2,1)); p=NeuralEcology(Float32[0 0 1 2;0 0 1 1;0 0 0 0];capacity=8)
    ws=Interior.NeuralWorkspace(w,p); frame=Interior.neural_frame!(ws,w,p)
    bins=Dict{CartesianIndex{3},Vector{Int}}()
    for i in eachindex(p.ids); push!(get!(bins,frame.cells[i],Int[]),i); end
    field=fill(.3f0,w.grid.size); other=copy(field); requests=Float32[.5,.7,.2,.9]
    reference=Interior._neural_allocate!(field,bins,requests)
    actual=Interior._neural_allocate!(other,ws,requests)
    @test actual==reference
    @test other==field
    # Existing capacity survives population shrink/growth through a replacement core.
    p2=NeuralEcology(zeros(Float32,3,8);capacity=8); Interior.workspace_prepare!(w,p2,ws)
    backing=ws.inputs
    p3=NeuralEcology(zeros(Float32,3,5);capacity=8); frame=Interior.neural_frame!(ws,w,p3)
    @test ws.inputs===backing
    @test size(frame.inputs,2)==5 && length(frame.ids)==5
end
