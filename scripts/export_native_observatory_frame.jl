#!/usr/bin/env julia
# Detached CPU-only observation export. Never loads Interior or a device backend.
using Serialization, JSON3, SHA
length(ARGS)==3 || error("FRAME EXPECTED_SHA256 FRESH_DIRECTORY required")
source,expected,out=ARGS
bytes2hex(sha256(read(source)))==expected || error("source SHA256 mismatch")
ispath(out)&&error("fresh export directory required")
x=deserialize(source)
hasproperty(x,:schema)&&x.schema=="interior-native-constructive-frame" || error("unsupported observation schema")
length(x.ids)==size(x.position,1)==length(x.radii) || error("inconsistent cellular observation")
all(isfinite,x.fields)&&all(isfinite,x.position)&&all(isfinite,x.radii) || error("nonfinite observation")
mkdir(out);arrays=Dict{String,Any}();offset=Ref(0)
open(joinpath(out,"state.bin"),"w") do io
    for name in (:fields,:ids,:parents,:lineages,:generations,:position,:radii,:mass,:energy,:metabolites,:bond_ids,:hidden,:adaptation,:actions,:catalytic_expression,:last_step_secretion)
        hasproperty(x,name)||continue
        value=getproperty(x,name)
        data=name in (:ids,:parents,:lineages,:generations,:bond_ids) ? Int32.(value) : name==:energy ? Float64.(value) : Float32.(value)
        arrays[string(name)]=Dict("offset"=>offset[],"bytes"=>sizeof(data),"shape"=>collect(size(data)),"type"=>eltype(data)==Int32 ? "i32" : eltype(data)==Float64 ? "f64" : "f32")
        write(io,data);offset[]+=sizeof(data)
    end
end
bytes2hex(sha256(read(source)))==expected || error("source changed during export")
meta=Dict("schema"=>"native-garden-observatory-frame-1","time"=>x.time,"step"=>x.step,
    "shape"=>collect(x.shape),"native_shape"=>collect(x.native_shape),"extent"=>collect(x.extent),
    "population"=>length(x.ids),"field_names"=>string.(x.field_names),"field_units"=>x.field_units,
    "field_sampling"=>x.shape==x.native_shape ? "native voxel concentrations" : "explicit concentration block means",
    "source_sha256"=>expected,"binary"=>"state.bin","binary_sha256"=>bytes2hex(sha256(read(joinpath(out,"state.bin")))),
    "order"=>"column-major","endian"=>ENDIAN_BOM==0x04030201 ? "little" : "big","arrays"=>arrays,
    "secretion_note"=>"Actual realized secretion in the final completed physical step only; not an interval total.")
open(io->JSON3.pretty(io,meta),joinpath(out,"frame.json"),"w")
