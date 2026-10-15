# Allocation-only correction; trajectory/timing evidence remains unchanged.
include(joinpath(@__DIR__,"../scripts/device_cellular_setup.jl"))
include("device_memory_dense.jl")
include(joinpath(@__DIR__,"../scripts/device_habitat_pilot_fixture.jl"))
const DMD=DeviceMemoryDense
host,hp=device_habitat_pilot_fixture(CPU();size=16,founders=8,capacity=64,variant="closed")
cp=Interior.device_habitat_checkpoint(host,hp)
newworld()=World(host.grid,backend,NamedTuple(),NamedTuple(),0.,host.seed)
baseline=Interior.device_habitat_restore(newworld(),cp)
specialized=DMD.restore(newworld(),cp)
Interior.KernelAbstractions.synchronize(backend)
isstorage=x->x isa AMDGPU.ROCArray
key=x->x.buf[]
bytes=x->sizeof(x.buf[])
measure(x)=(logical_unique_array_bytes=DMD.retained_bytes(x,isstorage),
    backing_allocation_bytes=DMD.retained_bytes(x,isstorage;storage_key=key,storage_bytes=bytes))
a=measure(baseline);b=measure(specialized)
report=(;schema="dense-memory-allocation-correction",runtime=device_cellular_runtime(),source=DMD.SOURCE_PIN,
    driver_sha256=bytes2hex(Interior.SHA.sha256(read(@__FILE__))),size=16,capacity=64,baseline=a,specialized=b,
    saved_bytes=a.backing_allocation_bytes-b.backing_allocation_bytes,
    saved_bytes_per_voxel=(a.backing_allocation_bytes-b.backing_allocation_bytes)/16^3,
    vram_used=parse(Int,strip(read("/sys/class/drm/card2/device/mem_info_vram_used",String))),
    scope="initial canonical checkpoint; no evolution/timing; logical array identities versus shared ROCArray managed allocation identities")
open(io->JSON3.pretty(io,report),ARGS[1],"w")
