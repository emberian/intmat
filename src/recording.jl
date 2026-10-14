# Recording is an observation boundary: it never calls advance!, seeds an RNG,
# or modifies a world/process. Captures must occur between complete step! calls.
function _record_value(value, key, dir, arrays)
    if value isa AbstractArray{<:AbstractFloat}
        T = eltype(value)
        T in (Float16,Float32,Float64) || throw(ArgumentError("unsupported floating element type $T at $key"))
        host = Array(value)
        bits = 8sizeof(T)
        filename = key * ".f" * string(bits)
        open(joinpath(dir, filename), "w") do io
            write(io, host)
        end
        entry = Dict("file"=>filename, "dtype"=>"float"*string(bits), "shape"=>collect(size(host)),
                     "source_dtype"=>string(eltype(value)))
        arrays[key] = entry
        return entry
    elseif value isa AbstractArray{<:Integer}
        # Persistent IDs stay exact, including values not representable in Float32.
        return Dict("dtype"=>string(eltype(value)), "shape"=>collect(size(value)),
                    "values"=>collect(value))
    elseif value isa Union{Number,AbstractString,Symbol,Bool} || value === nothing
        return value isa Symbol ? string(value) : value
    elseif value isa Tuple || value isa NamedTuple
        return [_record_value(v, key*"-"*string(i), dir, arrays) for (i,v) in enumerate(value)]
    elseif value isa ParticleRule
        return Dict(string(n)=>_record_value(getfield(value,n),key*"-"*string(n),dir,arrays)
                    for n in fieldnames(typeof(value)))
    end
    throw(ArgumentError("record_state cannot represent $(typeof(value)) at $key"))
end

"""
    record_process_fields(process)

Return a NamedTuple containing every persistent state array and effective rule
parameter to record. Custom processes must explicitly implement this hook;
unsupported processes fail rather than silently losing hidden state. Values may
be numeric/string/symbol scalars, tuples, or numeric arrays. Scratch buffers may
also be included. The hook must not mutate the process.
"""
record_process_fields(p) = throw(ArgumentError("record_state does not support process $(typeof(p)); implement record_process_fields"))
function record_process_fields(p::Union{NeuralField,ParticleCarrier,Excitable,Deposition,Diffusion,WaveForcing})
    names = fieldnames(typeof(p))
    NamedTuple{names}(Tuple(getfield(p,n) for n in names))
end

function _record_process(p, index, dir, arrays)
    recorded = record_process_fields(p)
    recorded isa NamedTuple || throw(ArgumentError("record_process_fields must return a NamedTuple"))
    state = Dict(string(n)=>_record_value(getproperty(recorded,n),"process-$index-"*string(n),dir,arrays)
                 for n in propertynames(recorded))
    entry = Dict{String,Any}("index"=>index,"type"=>string(nameof(typeof(p))),"state"=>state)
    if p isa NeuralField
        entry["provenance"] = Dict("rule"=>"fixed recurrent tanh; no training",
            "weights"=>"effective recurrence and bias are stored in state",
            "initialization"=>"constructor seed and mode are not retained by NeuralField")
    elseif p isa ParticleCarrier
        entry["position_coordinates"] = "periodic cell coordinates; origin at negative domain boundary; physical position = position * spacing - extent / 2"
        entry["particle_axis"] = "columns; IDs, positions, hidden and storage share the same order"
    end
    entry
end

"""
    record_state(w, processes, dir; step)

Write a view-independent state observation to a new directory, returning its
manifest. Includes all world fields, process state and effective rules, physical
cell-center coordinates, simulation time, and implementation SHA-256 hashes.
Floating arrays preserve native Float16/32/64 precision in host binaries in
Julia column-major order; integer
IDs remain exact in JSON. Call between complete simulation steps. This is an
inspection/export format, not a restart checkpoint or a concurrent snapshot.
The destination must not exist; failure removes only the temporary export.
"""
function record_state(w::World, processes::Tuple, dir::AbstractString; step::Integer)
    step >= 0 || throw(ArgumentError("step must be nonnegative"))
    target = abspath(dir)
    ispath(target) && throw(ArgumentError("recording destination already exists: $target"))
    mkpath(dirname(target))
    temporary = mktempdir(dirname(target); prefix=".recording-")
    try
        KernelAbstractions.synchronize(w.backend)
        arrays = Dict{String,Any}()
        fields = Dict(string(n)=>_record_value(getproperty(w.fields,n),"field-"*string(n),temporary,arrays)
                      for n in propertynames(w.fields))
        process_states = [_record_process(p,i,temporary,arrays) for (i,p) in enumerate(processes)]
        coordinates = Dict(string((:x,:y,:z)[d])=>
            [(i-0.5)*Float64(w.grid.spacing[d])-Float64(w.grid.extent[d])/2 for i in 1:w.grid.size[d]]
            for d in 1:3)
        sources = Dict(basename(path)=>bytes2hex(SHA.sha256(read(path)))
                       for path in sort(filter(p->endswith(p,".jl"),readdir(@__DIR__;join=true))))
        manifest = Dict("schema"=>"interior-state", "schema_version"=>2,
            "step"=>step,"time"=>w.time,"seed"=>w.seed,
            "julia_version"=>string(VERSION),"backend"=>string(typeof(w.backend)),
            "byte_order"=>(ENDIAN_BOM==0x04030201 ? "little" : "big"),
            "array_order"=>"column-major; first index varies fastest",
            "grid"=>Dict("shape"=>collect(w.grid.size),"extent"=>collect(w.grid.extent),
                "spacing"=>collect(w.grid.spacing),"boundary"=>"periodic",
                "coordinates"=>coordinates,"centering"=>"cell"),
            "fields"=>fields,"processes"=>process_states,"arrays"=>arrays,
            "source_sha256"=>sources,"purpose"=>"state observation; no restart loader")
        open(joinpath(temporary,"state.json"),"w") do io
            JSON3.pretty(io,manifest)
        end
        mv(temporary,target)
        manifest
    finally
        isdir(temporary) && rm(temporary;recursive=true)
    end
end
record_state(w::World, process, dir::AbstractString; step::Integer) =
    record_state(w,(process,),dir;step)
