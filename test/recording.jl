using Interior: JSON3
struct RecordingTestProcess
    memory::Vector{Float32}
    rate::Float32
end
Interior.record_process_fields(p::RecordingTestProcess) = (memory=p.memory,rate=p.rate)
@testset "view-independent state recording" begin
    function recording_fixture()
        w=World(Grid(4,3,1;extent=(4,6,2));seed=28)
        seed_activity!(w)
        neural=NeuralField(w;seed=93,mode=:random)
        particles=ParticleCarrier(Float32[0.4 1.2; 0.7 2.1; 0 0];ids=[2^25+1,7])
        (w,(neural,Deposition(),particles))
    end
    dense,p=recording_fixture(); sparse,q=recording_fixture()
    mktempdir() do dir
        for i in 1:5
            step!(dense,p,0.01f0);step!(sparse,q,0.01f0)
            record_state(dense,p,joinpath(dir,"dense-$i");step=i)
            i==5 && record_state(sparse,q,joinpath(dir,"sparse");step=i)
        end
        @test dense.time == sparse.time
        @test all(n->getproperty(dense.fields,n)==getproperty(sparse.fields,n),propertynames(dense.fields))
        @test p[1].latent1 == q[1].latent1
        @test p[1].latent2 == q[1].latent2
        @test p[3].positions == q[3].positions
        @test p[3].hidden == q[3].hidden
        @test p[3].storage == q[3].storage
        @test p[3].ids == q[3].ids
        snapshot=joinpath(dir,"dense-5")
        m=JSON3.read(read(joinpath(snapshot,"state.json"),String))
        @test length(m.fields)==8
        @test m.step==5 && Float32(m.time)==dense.time
        @test collect(m.grid.coordinates.x)==[-1.5,-0.5,0.5,1.5]
        @test collect(m.grid.coordinates.z)==[0.0]
        @test m.processes[3].state.ids.values[1]==2^25+1
        @test length(m.source_sha256["recording.jl"])==64
        for (key,descriptor) in pairs(m.arrays)
            data=reinterpret(Float32,read(joinpath(snapshot,String(descriptor.file))))
            @test length(data)==prod(descriptor.shape)
            @test all(isfinite,data)
        end
        for n in propertynames(dense.fields)
            descriptor=getproperty(m.fields,n)
            data=reinterpret(Float32,read(joinpath(snapshot,String(descriptor.file))))
            @test data == vec(getproperty(dense.fields,n))
        end
        for (index,names) in ((1,(:latent1,:latent2)),(3,(:positions,:hidden,:storage)))
            for n in names
                descriptor=getproperty(m.processes[index].state,n)
                data=reinterpret(Float32,read(joinpath(snapshot,String(descriptor.file))))
                @test data == vec(getproperty(p[index],n))
            end
        end
        @test_throws ArgumentError record_state(dense,p,snapshot;step=5)
        @test_throws ArgumentError record_state(dense,(Ref(1),),joinpath(dir,"unsupported");step=5)
        @test !ispath(joinpath(dir,"unsupported"))
        custom=RecordingTestProcess(Float32[1,2],0.3f0)
        cm=record_state(dense,(custom,),joinpath(dir,"custom");step=5)
        @test cm["processes"][1]["state"]["rate"]==custom.rate
        @test reinterpret(Float32,read(joinpath(dir,"custom","process-1-memory.f32")))==custom.memory
        @test_throws ArgumentError record_state(dense,p,joinpath(dir,"negative");step=-1)
    end
end
