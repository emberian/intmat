using Test
include(joinpath(@__DIR__,"..","scripts/cell_program_observatory.jl"))
using .CellProgramObservatory
@testset "Bounded observation-only program recording" begin
    metadata=(label="Synthetic recorder contract fixture",scope="No simulation or learning claim")
    cell=(id=1,sensed_step=0,raw=[.2,.3],raw_next=[.21,.31],prediction_from_step=0,prediction_for_step=10,prediction_error=nothing)
    frame=(time=.02,step=1,cells=[cell])
    r=Recorder(metadata;ids=[1],max_frames=2,max_cells=2,max_bytes=20_000)
    @test capture!(r,frame)
    frame.cells[1].raw[1]=.9
    @test document(r).frames[1]["cells"][1]["raw"][1]==.2
    @test capture!(r,(time=.22,step=11,cells=[merge(cell,(id=2,))]))
    @test isempty(document(r).frames[2]["cells"])
    @test !capture!(r,(time=.42,step=21,cells=[cell]))
    @test document(r).recorder.stop_reason=="frame cap reached"
    @test !capture!(r,frame)
    @test length(document(r).frames)==2
    @test_throws ArgumentError Recorder(metadata;ids=[1,1])
    @test_throws ArgumentError Recorder(metadata;max_frames=true)
    @test_throws ArgumentError Recorder(metadata;max_bytes=20)
    @test_throws ArgumentError capture!(Recorder(metadata),(time=.1,step=1,cells=[cell,cell]))
    @test_throws ArgumentError capture!(Recorder(metadata),(time=.1,step=1,cells=[merge(cell,(raw=[.1],))]))
    backward=Recorder(metadata);@test capture!(backward,(time=1.,step=50,cells=[]))
    @test_throws ArgumentError capture!(backward,(time=.9,step=49,cells=[]))
    limited=Recorder(metadata;max_cells=1);@test !capture!(limited,(time=.1,step=1,cells=[cell,merge(cell,(id=2,))]))
    @test isempty(document(limited).frames)
    bytes=Recorder(metadata;max_bytes=2300);@test !capture!(bytes,(time=.1,step=1,cells=[merge(cell,(annotation=repeat("x",2000),))]))
    @test document(bytes).recorder.stop_reason=="byte cap reached"
    mktempdir() do dir
        path=joinpath(dir,"recording.json");@test write_observatory(r,path)==path
        @test_throws ArgumentError write_observatory(r,path)
    end
end
