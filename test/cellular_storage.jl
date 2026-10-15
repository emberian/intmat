using Test, Random, Serialization

function storage_physical_equal(a,b)
    a isa Interior.CellularGenomeView && (a=copy(a))
    b isa Interior.CellularGenomeView && (b=copy(b))
    typeof(a)==typeof(b) || return false
    if a isa Union{Float32,Float64}
        T=sizeof(a)==4 ? UInt32 : UInt64
        return reinterpret(T,a)==reinterpret(T,b)
    elseif a isa AbstractArray
        axes(a)==axes(b) || return false
        return all(storage_physical_equal(a[i],b[i]) for i in eachindex(a,b))
    elseif a isa AbstractDict
        Set(keys(a))==Set(keys(b)) || return false
        return all(storage_physical_equal(a[k],b[k]) for k in keys(a))
    elseif a isa Union{Number,AbstractString,Symbol,Nothing,Type,Function}
        return isequal(a,b)
    end
    all(storage_physical_equal(getfield(a,n),getfield(b,n)) for n in fieldnames(typeof(a)))
end

@testset "Indexed genome lifecycle matches dense physical state" begin
    rng=MersenneTwister(731)
    w=World(Grid(16,16;extent=(8,8,1));seed=731)
    resources=ResourcePool(w;potentials=(4,3,0),initial=(.5,.5,0),diffusion=(.02,.02,.01))
    positions=rand(rng,3,8).*reshape([8.,8.,0.],3,1)
    dense=CellularEcology(w,positions,resources;capacity=128,biomass=1.79,energy=2.,metabolites=.3)
    wi,indexed=deepcopy((w,dense));Interior._cellular_enable_indexed_storage!(indexed)
    @test indexed.weights.bank===indexed.plastic_marks.bank
    for step in 1:80
        if step==7
            permutation=reverse(collect(eachindex(dense.ids)))
            reorder_ecology!(dense,permutation);reorder_ecology!(indexed,permutation)
        elseif step==19
            for cells in (dense,indexed)
                cells.biomass[3]=.1;cells.metabolites[:,3].=0.;cells.energy[3]=0.
            end
        elseif step==23
            for cells in (dense,indexed)
                cells.biomass[1]=2.;cells.energy[1]=2.;cells.last_actions[10,1]=1f0
            end
        end
        step!(w,(dense.resources,dense),.1)
        step!(wi,(indexed.resources,indexed),.1)
        @test storage_physical_equal(dense,indexed)
        @test storage_physical_equal(w,wi)
    end
    @test indexed.births>8 && indexed.deaths>=1
    @test length(indexed.weights.bank.slots)==length(indexed.ids)
    @test indexed.weights.bank.next_slot-1==length(indexed.ids)+length(indexed.weights.bank.free_slots)
    # The bank is an explicit run representation; snapshots materialize values.
    snapshot=record_process_fields(indexed)
    @test snapshot.weights isa Array{Float32,3}
    @test snapshot.weights==copy(indexed.weights)
    arena=Arena(wi,(indexed.resources,indexed);dt=.1)
    mktempdir() do directory
        path=joinpath(directory,"indexed.jls");checkpoint!(arena,path)
        restored=restore_arena(path);cells=restored.processes[2]
        @test cells.weights.bank===cells.plastic_marks.bank
        @test cells.resources===restored.processes[1]
        @test storage_physical_equal(indexed,cells)
        arena_step!(arena);arena_step!(restored)
        @test storage_physical_equal(arena,restored)
    end
end

@testset "Indexed genome array and inactive-slot ownership" begin
    rng=MersenneTwister(92)
    weights=randn(rng,Float32,5,7,4);marks=randn(rng,Float32,5,7,4)
    bank=Interior.CellularGenomeBank(weights,marks)
    w=Interior.CellularGenomeView{:weights}(bank);m=Interior.CellularGenomeView{:marks}(bank)
    @test Base.mightalias(w,bank.weights)
    @test Base.mightalias(m,bank.marks)
    @test !Base.mightalias(w,m)
    @test !Base.mightalias(w,copy(w))
    @test w[:,:,2]==weights[:,:,2]
    @test copy(view(w,2:4,3:6,2))==weights[2:4,3:6,2]
    Interior._cellular_select_genomes!(bank,[4,1,3])
    @test copy(w)==weights[:,:,[4,1,3]]
    @test copy(m)==marks[:,:,[4,1,3]]
    slot=only(bank.free_slots);block=35
    fill!(view(bank.weights,(slot-1)*block+1:slot*block),NaN32)
    fill!(view(bank.marks,(slot-1)*block+1:slot*block),NaN32)
    @test all(isfinite,w) && all(isfinite,m)
    capacity=bank.capacity
    Interior._cellular_append_genomes!(bank,1)
    Interior._cellular_copy_slot!(w,4,1);Interior._cellular_copy_slot!(m,4,1)
    @test bank.capacity==capacity && isempty(bank.free_slots)
    @test w[:,:,4]==weights[:,:,4] && m[:,:,4]==marks[:,:,4]
    @test all(isfinite,w) && all(isfinite,m)
    inputs=randn(rng,Float32,7,4);a=zeros(Float32,5,4);b=similar(a)
    Interior.neural_controller!(a,w,m,inputs)
    Interior.neural_controller!(b,copy(w),copy(m),inputs)
    @test reinterpret(UInt32,vec(a))==reinterpret(UInt32,vec(b))
end

@testset "Birth RNG reuse preserves daughter streams" begin
    reused=MersenneTwister(0)
    for id in 1:20
        seed=string(42,":division:",id,":",id÷2)
        fresh=MersenneTwister(seed)
        Random.seed!(reused,seed)
        # Uniform mutation gates, Gaussian mutation/partition caches and the
        # Float64 placement stream all reset exactly between daughters.
        @test rand(fresh,100)==rand(reused,100)
        @test reinterpret(UInt32,randn(fresh,Float32,100))==reinterpret(UInt32,randn(reused,Float32,100))
        @test randn(fresh,3)==randn(reused,3)
    end
end
