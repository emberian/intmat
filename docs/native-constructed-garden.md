# A funded native 3D constructed garden

The next production candidate is a fresh 64³ `DeviceConstructiveEcology`, using
[the explicit recipe](../experiments/native/constructed-garden-64.toml). It brings
the stronger developmental feedback into a paid catalytic and sorptive habitat,
with complementary raw-food patches and no ongoing soluble-food source. Its
purpose is a substantial physical history with construction, turnover and
changing cell behavior—not a benchmark fixture or a claim of adaptation.

The world has physical extent 4³, 256 randomly stratified founders, eight hidden
channels, density 128 and physical exchange radius 0.18. Founder biomass is 1.6,
stored energy 1 and each intracellular precursor 0.12. Motility speed 0.06 is slower
than the earlier native performance fixture's 0.7, allowing a cell to remain near
material it builds. Slow adaptation has strength 1.5 and time 8, the explicit
regime that broadened developmental switching in the CPU experiments. That
observation motivates the choice; it does not guarantee the same behavior under
the native law or this different 3D environment.

## Material and energy funding

Initial free-food concentrations are 0.3 each; raw-food means are 2 each. Raw A/B
have opposed smooth masks `1 ± 0.7*sin(2πx/Lx)*sin(2πz/Lz)`, applied to initial
raw inventory and later supply. Both remain available everywhere. Each raw
species receives **one total unit of mass per world-time**, normalized across
its mask, from a finite 1500-unit reservoir. These rates/budgets are world totals,
not concentrations or voxel rates. There is no soluble source or external sink.
Waste and rounding residues remain accounted material; neither restores fuel.
Raw and free precursors both have potential 4 per mass, so catalytic conversion
releases accessibility without creating chemical energy.

The unrounded initial material is 765.44. Including both finite source reservoirs
gives 3765.44, whose deliberately loose minimum-live-mass bound is 15062 cells at
minimum biomass 0.25. Capacity 16384 exceeds this bound, and the factory checks the
actual imported inventory. The guard therefore cannot act as the intended
carrying capacity. Actual source budgets last nominally until time 1500; the
requested observation horizon is 1200, or the earlier wall/artifact deadline.
The flux's ideal usable energy is 5.2 work units per time before maintenance,
motion, construction and other losses. This constrains sustainable biomass; it
is not a population target. Basal maintenance is explicitly 0.008 in this fresh
regime: the ideal flux alone would cover about 650 biomass before neural,
motor and construction costs, roughly 8% of the world volume at density 128.
The earlier 0.025 candidate would cover only about 208 biomass, below the
founder biomass 409.6. Lower maintenance gives construction and persistent
contacts more physical room without enlarging finite budgets or rewarding
colonies; it is an analytical design choice, not an empirically selected optimum.

Shared anabolism proportionally funds desired growth and construction from the
same complementary precursor/work budget. Catalyst secretion rate 0.12 costs 0.2
work per mass; catalytic rate 6, half-saturation 0.2 and decay 0.008 retain the
working construction regime. The paid catalyst also supplies binding capacity 6,
adsorption 4, desorption 0.12 and transport resistance 20. Soluble diffusivities
are 0.006/0.006, waste 0.001; raw food is immobile. The material is consequential
through access, retention and transport. No lineage-specific benefit or colony
reward is added.

## Native law and controls

Fields, bodies, geometry and controller calculations use the explicit native
Float32 law, with conservative Float64 local transfer calculations on AMD.
This new recipe explicitly selects **Float64 cellular energy accounts**: small
affordable motor costs accumulate in the battery rather than disappearing below
its Float32 resolution. There is no imposed minimum action cost. The setting is
part of the source-pinned checkpoint law and is unsupported on Metal; the
factory's default remains Float32 for recipes that omit it. Full frames and
checkpoints retain the battery type, while detached optical exports use Float32
energy values for display. It is not a bitwise CPU continuation. Native junction
exchange is synchronous and donor-limited; contacts use mutual degree-limited
proposals; mechanics uses 16 fixed, checked substeps and conservatively funded
motor work, rather than the CPU adaptive contact integrator. CPU-only
`MatrixDrag` is rejected by this factory. The existing preactivation/adaptation
law is retained.

The submitted Float32 timestep is 0.001 and must pass the coupled diffusion and
controller limit. Contact force/stiffness and travel guards remain active on
every substep. A latched failed step is never resumed or written as a complete
checkpoint. Physical support is materialized; bin width 0.25 uses the existing
conservative search based on actual radii. Radius bound 2 is a separate numerical
limit, not a prescribed body radius.

The [runner](../scripts/run_native_constructed_garden.jl) has a five-hour wall
bound, a hard 2026-10-15 20:50 UTC cutoff and a 4 GB artifact budget. It checks status
every eight steps, records scalar observations every 30 wall seconds, and retains
initial/final plus two rotating complete checkpoints. Full ledger audits occur
at checkpoints; relative accounting tolerance 0.0001 is an operational stop guard,
with actual errors and absolute tolerances retained. Failed or stopped runs keep
their previous complete checkpoints and diagnostics.

All living cells' positions, physical radii, bonds, developmental variables,
actions and catalytic readouts are captured. Chemistry is recorded as explicit
32³ block means every five model-time units and actual native 64³ fields every 50.
No state interpolation is used. Snapshots' per-cell secretion is the last
physical step, not the whole capture interval. Source, driver, profile, runtime
and checkpoint hashes accompany the recording. Resume requires the same source
and recipe and writes a new output directory.

## Production gate

No production launch is implied by this recipe. The positive-source reduction
must be integrated and pass its native gate first, followed by a bounded
whole-composition smoke and warmed throughput measurement of this exact recipe.
Earlier 64³ coupled execution measured 12.53 ms per step with the axis support;
128³ materialized execution measured roughly 30 ms after the source reduction.
Neither is a measurement of this 16384-capacity recipe. A planning range of 40–80
whole steps per second corresponds to 720–1440 model-time units in five hours at
this timestep, before captures/checkpoints. The final horizon is observed,
never inferred from this estimate.

One longer world is the default. Independent replicas are useful only if the
native ensemble adapter passes the same source/ownership/restart gates and its
measured aggregate throughput leaves a meaningful physical horizon. Smaller
performance fixtures and unintegrated reservoir approximations are not
substitutes for this recipe.

## Launched source and observation

The fresh garden launched on 2026-10-15 at approximately 16:50 UTC, from exact
commit `4b747919ec56f535f7c37bc7c666bd423316162d`. It is independent of the three
older CPU histories; their sources, states and deadlines remain unchanged.

- Runtime/source stage: `/othersys/tmp/intmat-probe-20261014/native-garden-4b747919ec56`.
- Output: `/tank/chreatures/probes/intmat-20261015/native-garden-4b747919ec56`.
- Initial Julia PID: `4182590`; supervisor PID `4182572`.
- Controls: CPUs 10–11, Julia 2 / GC 1 / BLAS 1, nice 15, 8 GiB host-memory cap,
  no swap, AMD optlevel 1. This is a host-memory cap, not a VRAM reservation.
- Bounds: model time 1200, hard 20:50 UTC, 4 GB artifacts; whichever comes first.

The exact-source native gate completed 1,100 steps at 129.33 steps/s
(7.732 ms/step), with 256 cells. Material residual was −1.892e−5 and potential
residual −4.942e−5. Strict native restoration and eight future steps were exact;
failure latching and rejection of partial checkpoints passed. Systemwide VRAM
peaked at about 1.81 GiB during that gate. This is a warmed gate measurement;
the long trajectory includes more frequent status checks and observations.

The immutable integration archive passed 18,804 main-suite assertions and 299
packed-world assertions. Local evidence is under
`outputs/integration-4b74791-20261015` and
`outputs/native-garden-launch-20261015`; the remote stage retains the gate,
source hashes, launcher, cgroup evidence, `long.log` and its eventual exit code.

Read `run.json`, append-only `metrics.jsonl`, and `checkpoint-audits.jsonl`.
Only checkpoints named in a completed manifest entry with matching size/hash
are candidates for strict same-source restoration. Never continue a failed
partial state. A resumed experiment must write a fresh directory; retain the
original output and parent checkpoint hash. No automatic reseed or rescue fuel.
The [read-only observatory](native-garden-observatory.md) consumes completed
frames and does not alter the simulation.


## Completed trajectory, 15 October 2026

The original source-pinned run stopped at model time659.992 after659,992 steps
on its unchanged completed-state material-accounting guard. It did not become
extinct: population fell from259 to43, remained near43 while stored energy
recovered, then reached93 by the final observation. There was no reseeding or
refueling intervention. Final material residual+.247363 exceeds tolerance.206432;
energy residual−.131032 remains inside.687488. The final checkpoint's exact bytes
and SHA were independently verified; no continuation was attempted. See
[numerical diagnosis](native-garden-accounting-diagnosis.md). These numerical
limits qualify the whole trajectory, including its apparent recovery.

| Observed model time | Cells | Reciprocal bonds | Largest bonded component | Mean stored energy |
|---|---:|---:|---:|---:|
| 0 | 256 | 0 | 1 | 1.000 |
| 80 | 259 | 397 | 140 | .643 |
| 170 | 52 | 4 | 2 | .000384 |
| 215 | 43 | 5 | 2 | .0374 |
| 300 | 43 | 7 | 2 | .245 |
| 400 | 43 | 11 | 5 | .814 |
| 500 | 44 | 20 | 5 | 1.611 |
| 600 | 64 | 58 | 9 | 1.978 |
| 659.992 | 93 | 100 | 17 | 1.941 |

All133 unique recorded times are preserved in the completed observatory, with
native64³ or explicitly averaged32³ material fields. The
[compact observation record](native-garden-completed-observations.json) pins
metadata and binary hashes. Components require reciprocal bonds between live
identities. Aggregation, fragmentation and renewed growth are measured features
of this history, not evidence that an adaptive multicellular organism evolved.
The earlier t100 food-access and through-t215 neural analyses are deliberately
bounded windows; neither represents the terminal outcome.
