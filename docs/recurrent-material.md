# Recurrent material: a GPU-native 3D experiment

`RecurrentMaterial` is eight evolving channels at **every material voxel**, rather
than a sparse population of agents living over an unrelated field. Four recurrent
oscillator pairs have relative clocks 1, .43, .17, .065. Each pair reads the
previous pair; all eight channels communicate through nearest-neighbor symmetric
fluxes. A smooth spatial transition joins two constitutions with different
rotational coupling. Both constitutions share an explicitly identical eight-channel
basis; this is compatible state communication, not invented latent translation.

Fast activity writes solid by converting finite nutrient. The opposite part of
the cycle erodes solid back into nutrient. Solid persists across many fast cycles,
reduces communication permeability, shifts recurrent thresholds and strengthens
cross-timescale coupling. Nutrient also diffuses conservatively. This closes the
activity → material → communication/threshold → future activity loop. Uniform
nutrient concentration one and small seeded latent noise initialize the run; no
Gaussian world imprint or timed external waveform drives the resulting patterns.
The rotational prior and two constitutions are designed. Weights are fixed;
nothing here was trained, selected or evolved genetically.

## Execution and representation

```
julia --project=gpu scripts/run_recurrent_material.jl 128 1200 outputs/my-material
RECURRENT_PROFILE=slow-memory julia --project=gpu scripts/run_recurrent_material.jl 256 1200 outputs/my-large-material
julia --project=gpu scripts/benchmark_recurrent_material.jl metal 256 400 outputs/my-benchmark.json
```

The fused KernelAbstractions kernel updates eight recurrent channels, conservative
transport, remodeling and published activity in one launch. Latent and resource
buffers ping-pong without whole-volume copies. GPU state remains resident between
observations; there are no per-step host reductions. Driver launch allocations
still occur: the measured Metal path allocated about 3.8KB of **host** memory per
step and AMD about 17 bytes; this is not a claim of zero host allocation.

World nutrient/deposit array handles are **borrowed until the next advance**.
Consumers must fetch current `world.fields` each stage; holding an earlier array
reference observes the previous buffer. The process retains both resource buffers
and rejects a world whose resource identities were replaced. Do not compose this
process with `ActiveMedium`, which retains a deposit alias, or another owner of
these scratch buffers. Sequential GPU processes must use the same ordered backend
queue or synchronize at observation boundaries. `record_state` does synchronize.
`process_ports` records the recurrent state; material remains counted once in the
world inventory. The constructor bounds latent channels in [-1,1], maintained by the convex
update under the timestep guard. External latent edits must preserve this bound;
preflight checks shapes, backend and buffer ownership without costly GPU
reductions. The timestep guard depends on spacing, communication and writing
rates. A 256³ run with extent 256³ enlarges the physical domain at unit spacing;
it is not a refinement of the 128³ domain.

## Executed experiments, 2026-10-15

Three matched 128³ Metal runs reached t60 (1,200 steps, dt .05):

| Composition | Communication | Material threshold feedback | Writing / erosion | Final activity std | Final deposit std |
|---|---:|---:|---:|---:|---:|
| Braided | .12 | 1.2 | .25 / .06 | .435 | .040 |
| Slow memory | .28 | 2 | .5 / .035 | .491 | .172 |
| Threshold feedback off | .12 | 0 | .25 / .06 | .448 | .031 |

The off variant retains deposit-dependent transport, so it is a **threshold and
cross-timescale feedback ablation**, not removal of every material feedback law.
Its final activity differs from matched braided activity by RMS .184. Persistent
remodeling therefore changes future internal state; the result is not an optical
recoloring. Slow memory generated larger spatial variation in solid. A 256³ slow
memory run reached the same t60 with activity std .494 and deposit std .172.
All recorded latent fields remained finite. Terminal Float64 material sums drifted
by .0106/.0135 units on 128³ and .116 units on 256³, about 5–7×10⁻⁹ of the initial
inventory; Float32 reductions can instead show residuals of a whole resource unit.

A native AMDGPU run on hbox reached t150 for all three 128³ compositions. Final
activity std was .517/.452/.476 and deposit std .036/.128/.034. These are longer
runs, not a cross-hardware same-time reproducibility claim. The 8³ ten-step
CPU/GPU arithmetic gate had maximum latent difference 2.98×10⁻⁸ on AMD and
5.96×10⁻⁸ on Metal.

Warm execution measured 695 million voxel-steps/s on AMD at 128³, and 244 million
on Metal at 256³ (synchronization after each 20-step block). These are eight-channel
material updates, not the full cell ecology benchmark. Observation-inclusive
trajectories measured about 211–230 million on Metal128³, 214 million on Metal256³,
and 491–570 million on AMD128³. Performance comparisons across these sizes include
different launch and transfer overheads. A separate copy-based Metal measurement
was contended with a larger run and is not used as evidence of a copy penalty.

## Retained artifacts

- `outputs/recurrent-material-128/`: three full 3D trajectories (activity/deposit
  every 100 steps), latent checkpoints, `activity-trajectory.gif`, a slice image,
  metrics and `loaded-model-source.jl`.
- `outputs/recurrent-material-256/slow-memory/`: selected 16.8-million-voxel
  trajectory and four full state observations. Raw arrays are Float32 in Julia
  column-major order; shape 256³. `warm-benchmark.json` lives one level above.
- `outputs/recurrent-material-amd128/`: retrieved AMD metrics and exact copied
  model/runner source. Full remote trajectories remain under
  `/othersys/tmp/intmat-probe-20261014/recurrent-material-20261015/trajectory-128/`
  on hbox. The pre-existing source-dose snapshot was not modified.

The older Metal trajectories used an equivalent recurrence kernel with a fixed
.05 guard and no explicit resource-identity check; their loaded source is retained
beside the outputs. Current source adds a spacing-aware guard and ownership
preflight without changing the recurrence. Live source hashes in old state
manifests can reflect concurrent guard edits; use retained loaded-model-source.jl
for the actual experiment implementation. New runner invocations copy source
before execution and retain the exact runner too.

These experiments exhibit structured ongoing activity and persistent material
memory. They do not establish open-ended adaptation or a rich learned behavioral
repertoire. A broader family of recurrent constitutions, composition with acoustic
transport, and selection of useful material transformations are subsequent work.
