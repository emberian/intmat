# Working hbox Torch route (15 October 2026)

The earlier successful Torch environment still works, including fused
`torch.compile` kernels and GPU real FFTs. It is separate from the retired
`/home/hbox/h1-ghost/venv` and from the new Julia AMDGPU experiments:

```sh
HSA_OVERRIDE_GFX_VERSION=10.3.0 \
  /tank/chreatures/envs/rocm-dev/bin/python -u PROGRAM.py
```

The recovered environment is Torch **2.9.1+rocm6.3**, HIP
**6.3.42134-a9a80e791**, using the RX 6750 XT's 12,868,124,672 bytes of VRAM.
Torch uses its `cuda` API for this ROCm GPU and reports `AMD Radeon Graphics`
under the required per-process gfx1030 override. Its wheel bundles HIP,
rocBLAS and rocFFT libraries; the missing system-wide ROCm installation was not
evidence that this environment was unavailable. No drivers or existing
environments were modified.

## Recovery evidence

`/Users/ember/dev/chreatures/docs/COMPUTE.md` explicitly records the retirement of
the H1 environment and names this replacement. The actual surrounding earlier
conversation is recoverable with:

```sh
cv show codex:01a0732d-671a-76f0-beb5-7c506b86ef03 --around 300 --context 10
```

That September exchange records a completed 256-circuit GPU acquisition, plus
NumPy/Torch parity of 3.58e-7. It led to the environment; the new measurements
below establish current operation independently. The initial current arithmetic
probe produced the exact expected sum 357389824 for squared Float32 integers
0–1023. The recurring `/opt/amdgpu/share/libdrm/amdgpu.ids` warning did not prevent
any of these kernels from completing.

## Current bounded measurements

Five warmed synchronized repetitions per size; Float32 periodic seven-point
diffusion stencil `u + .05*(six_neighbors - 6u)`. It is a representative thermal
operator, without moving source, cooling, dose, material conversion or arena
ecology. CPU/GPU eager output agreed exactly on 16³; compiled output differed by
at most 5.96e-8. Compilation took 51.49 seconds and is excluded from these times.

| Grid | Eager stencil | Compiled stencil | Resident rFFT + inverse | Upload one field | Download one field |
| --- | ---: | ---: | ---: | ---: | ---: |
| 64³ | .0921 ms | .0558 ms | .0941 ms | .0561 ms | .0491 ms |
| 128³ | .5798 ms | .1749 ms | .5338 ms | .7519 ms | .7393 ms |
| 256³ | 7.1008 ms | 1.2946 ms | 4.7430 ms | 6.2429 ms | 5.8129 ms |

The compiled 256³ stencil advances 12.96 billion voxel outputs/second and is
5.49× faster than the eager expression in this bounded comparison. Actual memory
traffic was not measured, so this is not a claim of peak hardware bandwidth.
Maximum FFT round-trip absolute error was 7.15e-7 for inputs in [0,1]. First FFT
plan/round-trip cost approximately 3.53 seconds per new shape and is excluded
from warmed timing. Peak allocated VRAM at 256³ was 538,984,448 bytes, including
the FFT state and probe buffers. These are component timings, not a complete
writer/arena speedup or a measured 1152³ run. Full values are in
`gpu/hbox-torch-verification.json`.

The probe source is `scripts/probe_hbox_torch.py`; a private working copy and
caches are under `/tank/chreatures/probes/intmat-20261015`. Reproduce from hbox:

```sh
SWARM_MEM_MAX=8G /home/hbox/bin/swarm-build timeout 360s env \
  HSA_OVERRIDE_GFX_VERSION=10.3.0 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 \
  TORCHINDUCTOR_COMPILE_THREADS=1 MAX_JOBS=2 \
  TORCHINDUCTOR_CACHE_DIR=/tank/chreatures/probes/intmat-20261015/inductor \
  TRITON_CACHE_DIR=/tank/chreatures/probes/intmat-20261015/triton \
  PYTORCH_KERNEL_CACHE_PATH=/tank/chreatures/probes/intmat-20261015/pytorch \
  /tank/chreatures/envs/rocm-dev/bin/python -u \
  /tank/chreatures/probes/intmat-20261015/probe_hbox_torch.py \
  --compile --sizes 64,128,256 --repetitions 5 \
  --output /tank/chreatures/probes/intmat-20261015/results.json
```

This uses the existing shared memory-control slice, an 8 GiB host-memory cap,
no swap, nice 15, the existing CPU affinity policy, two Torch CPU threads and
one compilation task. Initial import/device startup was slower cold; a second
arithmetic process imported in .76 seconds and completed GPU initialization and
arithmetic about 4.55 seconds after launch. No persistent compute job remains.

## Integration path

Use a persistent Torch worker that owns its device tensors. Julia can send
physical parameters and bounded batches of advances, then request only outputs
needed for export. A process protocol or an optional PythonCall environment
avoids repeatedly importing Torch and compiling kernels. Keep all intermediate
temperature, scratch and dose arrays on the device; the measured 12.06 ms
upload/download cost for a single 256³ field exceeds the 1.29 ms resident
stencil. DLPack alone would not repair the failed Julia GPU compiler; bridging
host-resident Julia arrays through it does not create zero-copy GPU state.

A faithful writer backend must compare the original cell coordinates, periodic
stencil, source trajectory/pulsing, cooling, timestep guards, threshold-dose
quadrature and material update against Julia fixtures over multiple steps.
The current microbenchmark establishes no such full writer equivalence.

1152³ contains 1,528,823,808 voxels, so each Float32 field occupies 6.115 GB.
Five fields require 30.58 GB; three require 18.35 GB before workspace. Both exceed
this GPU's VRAM. A useful hbox implementation therefore needs active tiles,
temporal batching, out-of-core storage or another representation rather than a
whole-domain three-array allocation. GPU FFT availability supports exploring
thermal propagators, but full-volume complex workspaces need explicit budgeting;
the 256³ timing must not be extrapolated into a claim that a 1152³ FFT fits.
