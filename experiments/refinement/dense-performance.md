# Native dense refinement

The implemented writing study uses a compact `World` containing only nutrient
and deposit, empty World scratch, and the unchanged ThermalWriter's three
arrays: **five Float32 volumes**. At 768³ its base state occupies **9.06 GB**,
compared with 34.43 GB for the default nineteen-array layout. This is a
writer-specific allocation contract; processes requiring additional fields
cannot use this carrier without allocating those fields.

Actual 96³ default versus compact runs for 120 steps produced exactly equal
nutrient, deposit, temperature, dose and clock. Reproduce with
`JULIA_NUM_THREADS=2 julia --project studies/verify_compact_writer.jl`.
The existing parametric World constructor makes this possible without changing
core solver code or duplicating thermal physics.

## Resolution, memory and physical time

All refinement recipes retain physical extent (2,2,2). Four times the resolution
per axis means 64 times as many cells; eight times means 512 times as many.
Explicit diffusion requires smaller timesteps proportional to spacing squared;
acoustic CFL scales with spacing. These are newly simulated fields, not
interpolated coarse volumes.

| Resolution | One Float32 field | Compact writing, 5 arrays | Default writing, 19 arrays | Default waves, 20 arrays |
|---|---:|---:|---:|---:|
| 96³ | 3.54 MB | 17.69 MB | 67.24 MB | 70.78 MB |
| 384³ | 226.49 MB | 1.13 GB | 4.30 GB | 4.53 GB |
| 768³ | 1.81 GB | 9.06 GB | 34.43 GB | 36.24 GB |

The completed 384³ writing run used the default carrier before compact allocation
was implemented; its fields and exact study source remain frozen. Writing 768
uses compact allocation, dt=.003125 below its .0033405 stability bound, and
1600 steps to physical time five. Writing 384 uses dt=.01 below .012792; wave 384
uses dt=.004 below .007518. The 768/384 comparison therefore includes temporal
as well as spatial refinement.

## Measured runs

Completed native 384³ outputs preserve raw Float32 fields, recipes, source
hashes, optical mappings and native-grid diagnostics. Rendering is separate.

| Specimen | Steps / dt / physical duration | Warm median step | Evolution time | Deposit max | Fraction deposit > .18 |
|---|---|---:|---:|---:|---:|
| writing-384 | 500 / .01 / 5 | 52.51 ms | 26.41 s | .9981 | .314% |
| waves-384 | 1500 / .004 / 6 | 137.69 ms | 207.75 s | .8163 | 8.10% |
| writing-768 | 1600 / .003125 / 5 | 409.73 ms | 651.77 s | .9981 | .312% |

The actual 384³ CPU writing pilot measured 1.072 seconds/step using two Julia
threads and one BLAS thread. Its Metal pilot measured 30.23 ms, while the full
run settled at 52.51 ms. The acoustic pilot likewise underpredicted full-run
cost (86.15 versus 137.69 ms). Compilation, allocation and export are excluded
from evolution times. Shared-host observations are not throughput promises.

The true 768³ compact Metal pilot measured 247.0 ms/step and passed finite checks.
The completed 1600-step run took 651.77 seconds (10.86 minutes), with a
409.73 ms median step and finite deposit, temperature and dose. Its three raw
fields remain native 768³ (5.44 GB), with the exact study source frozen.
At identical physical writing parameters, the completed 96³ reference's deposit
volume integral is .01292317 versus .01297516 at 384³, a .402% difference;
local maximum dose changes from .3130 to .3483. Similar integrated mass does
not imply equal resolved local structure.

## Exports and rendering

The runner reports its array memory plan before allocation, validates timestep,
and refuses existing output directories. `DENSE_PILOT_STEPS` bounds a pilot
without changing the grid; `DENSE_MEMORY_BUDGET_BYTES` sets the resident-plan
ceiling. Native Float32 output streams one field at a time: three fields occupy
679.48 MB at 384³ or 5.44 GB at 768³. Runtime/compiler/backend allocations and
staging need additional headroom beyond the base array plan.

The 768 render export uses a read-only native deposit field and surface-only
rendering, avoiding four additional 1.81 GB optical arrays and transport VDBs.
Optical carrier/material choices are display mappings; the deposit=.18 surface
comes from the actual grid. No image dimensions or surface mesh imply a changed
simulation resolution. The .018 beam width spans about 3.5 cells per standard
deviation at 384³ and 6.9 cells at 768³, versus less than one cell at 96³.

The final 768³ deposit integral is .01286222, versus .01297516 at 384³
(-.870%). Maximum dose is .348536 versus .348261. Because the stable timestep
was also refined, these differences are combined spatial/temporal observations,
not an isolated spatial convergence order. All native field sizes and finite
values were checked by bounded Float64 CPU scans.
