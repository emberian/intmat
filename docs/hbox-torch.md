# Working hbox Torch route (15 October 2026)

**Diagnostic evidence only.** Ember requires Julia for the simulation algorithms.
The Torch probes below are retained as evidence about the working hardware,
runtime libraries and numerical comparison. They are not an approved simulation
backend or a production execution route. Further algorithm work targets Julia.

Julia recovery subsequently passed the original Interior CPU/GPU comparison
(maximum absolute field error 3.73e-9) using normal package images and the
working wheel's ROCm libraries read-only, without importing Torch. The current
Julia replay and controlled attribution are in [hbox-runtime.md](hbox-runtime.md).

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

## Implications for Julia recovery

The working wheel bundles ROCm 6.3 libraries under
`/tank/chreatures/envs/rocm-dev/lib/python3.12/site-packages/torch/lib`, and uses
`HSA_OVERRIDE_GFX_VERSION=10.3.0`. These provide a concrete runtime reference for
recovering Julia AMDGPU. The earlier proposal to use a persistent Torch worker,
PythonCall or DLPack as an algorithm backend has been withdrawn in accordance
with the Julia-only requirement. The measured transfer costs remain useful
evidence for keeping eventual Julia state resident on the GPU.

A faithful writer backend must compare the original cell coordinates, periodic
stencil, source trajectory/pulsing, cooling, timestep guards, threshold-dose
quadrature and material update against Julia fixtures over multiple steps.
The component microbenchmark above establishes no such full writer equivalence.

### Historical diagnostic writer comparison

`scripts/probe_hbox_torch_writer.py` was used to compare the opt-in three-array
`DoseThermalWriter` equations with a fused resident Torch kernel. It retains
the explicit Euler temperature update, periodic finite-difference Laplacian,
Gaussian source, cooling, nonnegative temperature clamp and right-endpoint
threshold-dose quadrature. Julia supplies exact Float32 midpoint source positions
and pulsed powers in a device-resident table, avoiding host trigonometric
differences and per-step parameter transfers. Temperature/scratch pointers swap
after each step; dose updates in place. This removes the full-volume scratch
copy without changing the next step's physical state.

Against the Julia-generated 12³, 25-step fixture, saved frames at steps 1, 7, 15
and 25 passed, including pulsing and source shutdown. Maximum absolute errors
were 2.24e-8 temperature and 1.63e-9 dose. Analytic material reconstruction
matched within 5.96e-8 nutrient and 6.98e-10 deposit. This comparison is against
the new dose-owned material representation, not the older repeated nutrient/
deposit exchange rounding.

The same resident code ran **384³ for 25 steps in .12505 seconds (5.002 ms/step)**,
using 683,225,600 bytes peak allocated VRAM. Fixture compile/check took 43.19
seconds and large-shape compile/warmup took 6.23 seconds, both excluded from the
timed advance. The benchmark reused the fixture's .003 timestep, parameters and
25 source rows; it is a short writer test, not a five-unit production specimen
or a whole-arena speedup. Evidence is in `gpu/hbox-torch-writer-verification.json`.

Generate the fixture locally with
`julia --project scripts/dose_writer_fixture.jl FRESH_OUTPUT_DIRECTORY`, copy it
and `scripts/probe_hbox_torch_writer.py` to a fresh hbox probe directory, then
use the same runtime/cache/resource-control environment as above:

```sh
/tank/chreatures/envs/rocm-dev/bin/python -u probe_hbox_torch_writer.py \
  dose-fixture/fixture.json --n 384 --steps 25 --output writer-results.json
```

1152³ contains 1,528,823,808 voxels, so each Float32 field occupies 6.115 GB.
Five fields require 30.58 GB; three require 18.35 GB before workspace. Both exceed
this GPU's VRAM. A useful hbox implementation therefore needs active tiles,
temporal batching, out-of-core storage or another representation rather than a
whole-domain three-array allocation. GPU FFT availability supports exploring
thermal propagators, but full-volume complex workspaces need explicit budgeting;
the 256³ timing must not be extrapolated into a claim that a 1152³ FFT fits.
