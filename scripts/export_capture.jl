using Interior, JSON3, TOML

# Export immutable captures, never integrate or modify the simulation.
length(ARGS) >= 2 || error("usage: julia --project scripts/export_capture.jl RUN OUTPUT [TIME=5] [CARRIER=.08] [DEPOSIT_GAIN=.7]")
run_directory=abspath(ARGS[1]); output_directory=abspath(ARGS[2])
requested_time=length(ARGS)>=3 ? parse(Float64,ARGS[3]) : 5.0
carrier=length(ARGS)>=4 ? parse(Float32,ARGS[4]) : .08f0
deposit_gain=length(ARGS)>=5 ? parse(Float32,ARGS[5]) : .7f0
profile=TOML.parsefile(joinpath(run_directory,"profile.toml"))
n=Int(profile["simulation"]["resolution"])
records=[]
for line in eachline(joinpath(run_directory,"frames.jsonl"))
    isempty(strip(line)) && continue
    # The live runner can be appending a line; only complete JSON records qualify.
    endswith(strip(line),"}") || continue
    record=JSON3.read(line)
    paths=[joinpath(run_directory,"$(record.prefix)-$name.f32") for name in ("deposit","activity","nutrient")]
    all(path->isfile(path)&&filesize(path)==4n^3,paths) && push!(records,record)
end
isempty(records) && error("No complete scalar captures available")
selected=records[argmin([abs(Float64(record.world_time)-requested_time) for record in records])]
read_field(name)=reshape(copy(reinterpret(Float32,read(joinpath(run_directory,"$(selected.prefix)-$name.f32")))),n,n,n)
deposit=read_field("deposit");activity=read_field("activity");nutrient=read_field("nutrient")
mask=similar(deposit)
for I in CartesianIndices(mask)
    x,y,z=ntuple(d->2f0*(I[d]-.5f0)/n-1f0,3)
    r=sqrt((x/.87f0)^2+(y/.94f0)^2+(z/.80f0)^2)
    mask[I]=clamp((1f0-r)*30f0,0f0,1f0)
end
export_volume((density=carrier .+ deposit_gain.*deposit,activity=activity,structure=deposit,mask=mask,deposit=deposit,nutrient=nutrient),output_directory)
open(joinpath(output_directory,"capture-source.json"),"w") do io
    JSON3.pretty(io,(source_run=run_directory,capture_prefix=String(selected.prefix),requested_time=requested_time,actual_time=Float64(selected.world_time),shape=(n,n,n),simulation_extent=get(profile["simulation"],"extent",[2,2,2]),density_mapping=(carrier=carrier,deposit_gain=deposit_gain),viewing_crop=(axes=(.87,.94,.80),edge_rate=30),interpretation="Artistic optical extraction from a recorded periodic volume; no dynamics or physical boundary added"))
end
println("Exported ",selected.prefix," at actual t=",selected.world_time," to ",output_directory)
