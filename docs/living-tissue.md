# Living tissue as a continuum

`studies/living_tissue.jl` is a different factorization of active material: a
spatial concentration of tissue, rather than a population of mechanically
resolved individuals. Growth, spending, excitation, persistent construction and
transport all occur on the same regular lattice. It is useful when the question
concerns fronts, communication, tissue connectivity or material memory and does
not require individual genomes, bodies or contacts.

The first 96² specimen expands smoothly and eventually loses most excitation as
its finite supply runs down. The second, explicitly seeded 192² lace specimen
forms patchy fronts, closes interrupted strands, and leaves a slower matrix
around a changing network of voids. Both are retained. The initial anatomy is
deliberately drawn; subsequent concentrations and motion follow the stated law.
Neither attractive geometry nor sustained cycling demonstrates evolution.

## State and representation

Two packed, channel-major Float32 arrays hold seven channels. Julia's first
spatial index varies fastest; each channel is a contiguous volume. Persistent
simulation storage is **56 bytes per lattice site**; 256³ requires 939.5 MB and
512³ requires 7.52 GB before captures, backend overhead or renderer storage.
There are no per-site weight matrices, contact bins, lifecycle proposals or
mechanical iterations. The current generic World has no allocated fields;
all material belongs to the process and is included by `owned_material_mass`.

| Channel | Meaning | Transport |
| --- | --- | --- |
| `precursor` N | Building stock and fuel, potential 1 | Symmetric matrix-limited diffusion |
| `biomass` B | Living continuum tissue, potential 1 | Diffusion and optional bounded precursor-directed hopping |
| `matrix` M | Persistent tissue-built structure, potential 1 | Immobile |
| `spent` W | Exhausted material, potential 0 | Immobile |
| `excitation` x | Fast bounded local state | Communication across occupied tissue |
| `inhibition` y | Delayed negative feedback | Site-bound |
| `memory` q | Slower signed activity trace | Site-bound |

This storage removes whole categories of work rather than accelerating the
same cellular law. A lattice site is not an organism. No agent speedup follows
from a voxel throughput number. The stencil is O(number of sites), including
fully occupied tissue; it can still waste work on large empty domains. Explicit
transport remains constrained by grid spacing squared. Dense fields do not solve
sparse-domain scaling, arbitrary deformations or individual inheritance.

## Transfers, work and communication

The process reuses `RecurrentMaterial`'s symmetric face primitive,
`1 / (1 + 4*barrier*(M_left+M_right))`. Precursor and biomass fluxes use this same
face value on both sides. Excitation additionally multiplies it by each side's
bounded tissue support `B/(tissue_saturation+B)`. Precursor-directed biomass
hopping uses an upwind donor and the antisymmetric signal
`tanh(taxis_sensitivity*(N_neighbor-N_here))`. Its rate is bounded by `taxis`.
The explicit bound includes both diffusion and maximum outgoing hopping.

Each step transports first, then performs local nonnegative donor-limited
transfers. The transfer fraction `dt*rate/(1+dt*rate)` is bounded below one:

1. **N → W pays activity.** The requested amount is
   `dt*metabolic_work*B*(0.25+abs(x))`. Actual spending is limited by available
   precursor. The ratio paid/requested gates neural drive, growth and construction.
2. **N → B grows tissue.** Rate depends on fuel availability, existing biomass,
   a soft density capacity and the positive excitation phase plus a small basal
   term. Growth cannot create biomass where both biomass and its incoming flux
   are exactly zero.
3. **B → N turns tissue over.** Starvation raises this recycling rate.
4. **B → M constructs matrix.** Construction follows the paid positive phase.
5. **M → N erodes matrix.** Its slower passive erosion increases in the negative
   excitation phase. Original construction history is not an irreversible scar.

`N+B+M+W` conserves material. Under the declared unit-potential convention,
`N+B+M` is the remaining chemical potential and spending increases W. Recycling
does not manufacture high-potential material from W. The reported
`potential+spent` residual is algebraically the material residual under this
simple assignment; it is **not an independent chemical validation**. This is one
abstract precursor and one spent pool, without molecular stoichiometry, heat,
ATP or a force/work correspondence. The metabolic charge is a declared funding
convention for activity. Biased tissue transport is a phenomenological
constitutive flux, without a resolved chemomechanical energy balance.

All recipes here are closed periodic specimens with **finite initial fuel and
zero renewal budget**. No implicit chemostat replenishes N. Eventual fading or
collapse is a possible result of this law, not a bound on all active matter.

Writing `s=B/(K_B+B)`, `v=s*paid*N/(K_N+N)`, the local neural rule is

```
dx/dt = v * (tanh(g*x - h*y - memory_feedback*q
                  - matrix_feedback*M/(K_N+M)) - x)
        - (1-s)*x + matrix-and-tissue-limited communication
dy/dt = (x-y)/inhibition_time
dq/dt = (x-q)/memory_time
```

The positive phase changes actual material transfers. Matrix feeds back through
both neural thresholds and transport. The slow trace is Eulerian site memory;
it is not advected with biomass or inherited by individual descendants. The
chosen timestep bound preserves physical positivity and bounded neural state
under the stated explicit stencil. Boundedness is not a general accuracy or
nonlinear stability proof. The tests include a short halved-step comparison.

## Running and observing

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=. \
  studies/living_tissue.jl experiments/art/living-tissue-lace.toml \
  outputs/FRESH-LIVING-TISSUE cpu
```

The ribbon recipe is a cheaper first specimen. `cpu`, `amd`, and `metal` select
the same Julia KernelAbstractions kernel. Accelerator execution must use the
project's verified runtime and coordinated hardware slot. Backend availability
and actual device equivalence require a successful launch; source portability
alone is insufficient evidence. A bounded probe is provided:

```sh
julia --project=. studies/living_tissue_gpu_probe.jl amd FRESH-RESULT.json 128 256
```

The runner refuses existing destinations and checks planned resident and capture
budgets before starting. Initial/final observations use generic `record_state`;
the independent restart checkpoint stores all seven native channels, effective
parameters, grid, seed and time. Restore checks the exact study source plus the
shared grid/stepping and face-law source hashes. It never silently adopts a new
law. A pinned historical run needs its historical source/environment.

Every native frame contains seven named `.f32` files and `volume.json` with
shape, extent, spacing, time, periodic boundary, endianness and Julia/F array
order. A renderer can use B and M as distinct actual geometries, x as a signed
activity control and q as an independent slow optical modulation. Optical
density, color, emission and refraction remain explicit artistic mappings.
`section.u8` is a compact middle-section replay with fixed ranges recorded in
`configuration.json`; it is a display asset, not an exact checkpoint.

`studies/living_tissue_preview.py` draws a fixed-exposure sheet of six actual
sections. It performs no simulation. Its purple tissue, amber matrix and cyan
positive excitation are visual choices. The full native fields remain available
for different optics. Captures temporarily allocate host copies and recorders
allocate export storage; the 56-byte figure covers persistent simulation buffers.

## Focused checks and interpretation

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=. test/living_tissue.jl
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=. \
  studies/living_tissue_ablation.jl experiments/art/living-tissue-lace.toml \
  outputs/FRESH-TISSUE-ABLATION
```

The narrow tests cover a separate symmetric-flux oracle, material/potential
accounting, no-fuel denial of construction, conservative directed redistribution,
generic diagnostics and recording, composition preflight, exact restart and
continuation, same-site cycling and short timestep refinement. They do not
certify a biological interpretation.

The four-arm contrast crosses matrix barrier on/off with dynamic/static material
control from identical initial states. The static gate equals the full arm's
biomass-weighted average positive excitation sampled every 100 steps. It controls
approximate mean gating, not equal realized construction or metabolic work after
trajectories diverge. Connectivity measurements count periodic grid components
of `B+M > 0.15`, and components of its complement. Those threshold-dependent
counts are not individual organisms or a proof of a material boundary.

## Observed specimens and contrasts

[Compact results and source hashes](living-tissue-results.json) retain the
measured values and exact run locations. All local runs used Julia 1.12.4 and
one CPU thread. The runner timings include native capture/probe readback and
exclude initial construction and the separate warm-up.

| Specimen | Steps / time | Actual run rate | Outcome |
| --- | --- | --- | --- |
| 96² ribbons | 4,000 / 160 | 10.13 million site-steps/s | 7 fixed-site excitation sign changes; matrix peaks at 516, then falls to 465 as activity fades |
| 192² lace | 6,400 / 80 | 4.25 million site-steps/s | 8 fixed-site sign changes; matrix grows to 2,326; mean permeability falls to .368 |

Both retain 51 seven-channel native frames. Maximum relative material drift is
5.3e-7 for ribbons and 1.2e-6 for lace. These are different spatial resolutions,
parameters, seed geometries and timesteps; their timing ratio is not a scaling
or optimization comparison. Their reviewed contact sheets remain at
`outputs/living-tissue-20261015-{ribbons,lace}/sections.png`.

The matched 128² four-arm study reaches time 80 using 3,200 steps. Its constant
material gate is .04474373, estimated from the full arm rather than chosen to
force an outcome:

| Intervention | Final biomass | Matrix | Spent | Mean permeability |
| --- | --- | --- | --- | --- |
| Full feedback | 3,568 | 1,011 | 6,117 | .378 |
| Matrix barrier off | 3,715 | 1,042 | 6,297 | 1 |
| Mean material gate | 3,457 | 1,053 | 6,456 | .272 |
| Both interventions | 3,649 | 1,112 | 6,782 | 1 |

Removing the barrier changes final biomass by RMS .0417 per site. Replacing
phase control with its approximate mean changes biomass by RMS .0245 and matrix
by .0257. This establishes that the feedbacks materially affect this trajectory;
it does not rank adaptive success. The full arm's thresholded tissue changes
from 17 occupied components to one while its complement varies from one to
eleven components. At time 80 it has two complement components. Fixed-site
sign changes occur in all four arms. The dynamic gate is therefore consequential
without being the sole cause of neural cycling.

The first contrast using an arbitrary constant .5 is also preserved under
`outputs/living-tissue-20261015-ablation/`. It constructs much more matrix and
substantially suppresses biomass; the different average requested construction
confounds that comparison. The mean-gate contrast above is the useful replacement.

The focused CPU suite passes **27 assertions**.

## Native AMD gate

The frozen `4e2889c` study completed its reserved AMD gate on the verified
RX 6750 XT runtime. Fifty complete steps of a dense randomized 16³ fixture give
maximum CPU/device error **1.19e-7** across all seven channels. The 128³ sparse
and dense cases each run three 256-step batches after ten warm steps:

| Native 128³ state | Median time / 256 complete steps | Site-steps/s |
| --- | --- | --- |
| Sparse lace | .188765 s | 2.844 billion |
| Dense tissue | .193277 s | 2.778 billion |

The timed boundary includes `step!`, process validation, the whole material and
neural kernel, buffer swaps, clock updates and final synchronization. It excludes
initialization, compilation, captures and rendering. Persistent arrays occupy
117,440,512 bytes. No cellular-arena speedup is inferred. These are three short
samples per density on a shared device, not a long-run production bake or a
cross-backend art comparison. [Native results and pin](living-tissue-amd-results.json)
retain timings, array diagnostics, source hashes and the copied manifest hash.

The initial host selector encountered Julia's late-import world-age restriction
before any kernel launched. The successful retry preloaded AMDGPU and called
the unchanged frozen probe function; no fixture, parameters or kernel changed.
The current host selector resolves its newly imported backend in the current
world. The newer active-turnover process is separate and is not covered by this
device gate.
