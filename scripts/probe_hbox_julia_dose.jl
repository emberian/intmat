# Julia host launcher only: all thermal/material algorithms come from Interior.
# julia --project=AMD_ENV probe_hbox_julia_dose.jl FIXTURE OUTPUT [N=384] [STEPS=25]
using Interior, AMDGPU
const JSON3 = Interior.JSON3
length(ARGS) >= 2 || error("fixture JSON and output JSON required")
fixture_path=abspath(ARGS[1]); output_path=abspath(ARGS[2])
n=length(ARGS)>=3 ? parse(Int,ARGS[3]) : 384
steps=length(ARGS)>=4 ? parse(Int,ARGS[4]) : 25
16<=n<=768 && 1<=steps<=100 || error("bounded grid16..768 and steps1..100 required")
AMDGPU.functional() || error("AMDGPU runtime unavailable")
backend=AMDGPU.ROCBackend()
fixture=JSON3.read(read(fixture_path,String))
haskey(fixture,:writer_parameters) || error("fixture must include writer_parameters; regenerate with scripts/dose_writer_fixture.jl")
parameters=Dict{Symbol,Any}(Symbol(k)=>v for (k,v) in pairs(fixture.writer_parameters))
parameters[:path]=Symbol(parameters[:path]); parameters[:center]=Tuple(parameters[:center])
kwargs=(;parameters...)
shape=Tuple(Int.(fixture.shape)); extent=Tuple(Float32.(fixture.extent))
dt=Float32(fixture.sources[1].dt)
all(row->Float32(row.dt)==dt,fixture.sources) || error("fixed timestep fixture required")
function build_writer(dims)
    world=World(Grid(dims...;extent),backend,NamedTuple(),NamedTuple(),0.0,91)
    writer=Interior.DoseThermalWriter(world;
        initial_nutrient=fixture.initial_nutrient,initial_deposit=fixture.initial_deposit,kwargs...)
    world,writer
end
function expected_field(filename)
    data=Array{Float32}(undef,shape)
    open(io->read!(io,data),joinpath(dirname(fixture_path),filename))
    data
end
function check_field(actual,expected)
    values=Array(actual); difference=abs.(values.-expected)
    all(difference .<= 2f-6 .+ 2f-5.*abs.(expected)) || error("fixture field mismatch $(maximum(difference))")
    maximum(difference)
end
world,writer=build_writer(shape)
frame_map=Dict(Int(frame.step)=>frame for frame in fixture.frames)
errors=Any[]
started=time()
for step in 1:length(fixture.sources)
    step!(world,writer,dt)
    if haskey(frame_map,step)
        frame=frame_map[step]; q=writer.writer
        row=(;step,temperature_max_abs_error=check_field(q.temperature,expected_field(frame.temperature)),
            dose_max_abs_error=check_field(q.dose,expected_field(frame.dose)))
        push!(errors,row); println(JSON3.write(row));flush(stdout)
    end
end
nutrient=Interior.KernelAbstractions.zeros(backend,Float32,shape...)
deposit=similar(nutrient)
Interior.writer_material!(nutrient,deposit,writer)
material_errors=(;nutrient_max_abs_error=check_field(nutrient,expected_field(fixture.material_readout.nutrient)),
    deposit_max_abs_error=check_field(deposit,expected_field(fixture.material_readout.deposit)))
fixture_seconds=time()-started
world,writer=build_writer((n,n,n))
warm=@elapsed step!(world,writer,dt)
# Warm state is explicitly excluded; benchmark continues the same trajectory.
elapsed=@elapsed for _ in 1:steps
    step!(world,writer,dt)
end
q=writer.writer
result=(;schema="interior-julia-amdgpu-dose-probe",julia=string(VERSION),
    backend="AMDGPU",fixture=fixture_path,fixture_source_hashes=fixture.sources_sha,
    frame_errors=errors,material_errors,fixture_compile_and_check_seconds=fixture_seconds,
    warm_step_seconds=warm,benchmark=(;grid=[n,n,n],timed_steps=steps,dt,seconds=elapsed,
        seconds_per_step=elapsed/steps,clock=world.time,temperature_range=extrema(q.temperature),dose_range=extrema(q.dose)),
    scope="Existing Julia DoseThermalWriter and writer_material!; no reimplemented thermal algorithms")
open(io->JSON3.pretty(io,result),output_path,"w")
println(JSON3.write(result));flush(stdout)
