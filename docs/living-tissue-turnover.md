# Rooms without a plan

This is a separate continuum material law, `TurnoverTissue`, composed around
`LivingTissue`. It addresses a failure of the original seeded ribbons: excitation
could accompany broadening strands without reopening the material. Here opposing
phases spend precursor to build matrix or withdraw tissue, dissolve matrix and
move biomass. Nearly uniform starts develop repeated local material turnover.
This is neither individual-cell mechanics nor a claim of biological evolution.

The runnable model is `studies/living_tissue_turnover.jl`; the bounded exporter is
`studies/living_tissue_turnover_run.jl`. The generic process, recording, validation,
material-ownership and exact-source checkpoint contracts remain usable. The
wrapper declares its tissue child, preventing accidental double scheduling.

## Material and activity convention

Eight native Float32 fields are retained: precursor N, biomass B, matrix M, spent
material W, signed excitation x, inhibition y, memory q and immobile reserve R.
N/B/M/R have declared potential 1 per material unit; W has potential 0. Signal
channels are dimensionless controller states, not additional chemical species.
The model conserves N+B+M+W+R, including finite reserve renewal R→N. Its declared
potential decreases by paid N→W work. Material conservation alone would not prove
chemical stoichiometry, thermodynamics or the energetic cost of these controllers.

The original tissue stage performs bounded precursor/biomass conversion, passive
return and metabolic work. Its construction, erosion and unfunded taxis are
disabled here. A second, double-buffered stage renews N toward a target using R,
then requests `dt*phase_work*(B+M)*(positive+negative)`. Positive and negative
phase magnitudes are `max(0, ±x-phase_threshold)`. Only available precursor is
charged; the actual paid/requested fraction gates every phase action.

Positive phases transfer B→M. Negative phases transfer M→N and B→N and fund
directed B hopping toward precursor. All reaction transfers are positive and
donor limited. Hopping uses the old donor fields and paid fractions on each
antisymmetric face, so an arriving parcel cannot be spent on another outgoing
hop in that stage. Construction/withdrawal after hopping is a stated sequential
local conversion. Grid-dependent explicit stability bounds limit hopping and
diffusion. Matrix reduces symmetric face permeability; changing matrix therefore
changes both nutrient and signal transport. Inhibition diffuses faster than
excitation; memory relaxes more slowly than inhibition.

Signals are abstract control dynamics. Their exchange and integration are not a
calibrated mechanical work calculation. The metabolic and phase charges are
explicit conventions, not inferred molecular chemistry. There is no learning,
inheritance or selection mechanism in this study.

## Physical scaling and initialization

Field values are concentrations. Rates use simulation time units; diffusivities
use length²/time, with actual voxel spacing in every stencil. `physical_*`
inventories multiply concentration sums by voxel area in 2D or volume in 3D.
The 24×24 sheets therefore have the same initial physical material inventory
7407.36 at 64² and 128². The proposed 24³ volume has spacing .1875 and finite
reserve concentration 12, not twelve units per arbitrarily refined voxel.

The refined starts share one continuous, seeded Fourier field with 64 modes of
integer radius ≤12. Initialization computes exact cell averages, including sinc
attenuation and cell-center phase, using three CPU inverse FFTs. Refinement does
not redraw independent noise. Initial means are N=.6, B=.26; N/B have 2% noise,
x has amplitude .02, and matrix/spent/inhibition/memory start at zero. No tubes,
rooms, anatomical skeleton or time-dependent forcing are drawn into this start.

The frozen-coefficient 3-channel neural dispersion study is a wavelength
hypothesis. It freezes tissue/fuel and linearizes zero signals despite the full
material state evolving. At the chosen initial B=.26, gain2.8 and excitation
diffusivity .08, axis/diagonal discrete symbols predict peak wavelength about
3.1 at both grid sizes. This is not a nonlinear organization or equilibrium proof.

## Actual 2D observations

All pilots and failed regimes remain under `outputs/living-tissue-20261015-*`.
Compact reviewed results and source hashes are in
[`living-tissue-turnover-results.json`](living-tissue-turnover-results.json).
Only actual saved Float32 fields are used for these observations.

| Study | End time | Opening / reformation events | Result |
| --- | ---: | ---: | --- |
| Initial weak-turnover 64² | 160 | 0 / 0 | Activity without strict material reopening; rejected as the target regime. |
| Stronger reopening 64² | 160 | 323 / 323 | Mostly initial transient followed by broad breathing; insufficient spatial richness. |
| Regrowth 64², dt .04 | 160 | 29467 / 26135 | Repeated turnover but visible coarse-grid effects. |
| Same 64², dt .02 | 160 | 29475 / 26244 | Time refinement preserves turnover and similar expenditure. |
| Regrowth 128², dt .02 | 160 | 128164 / 112141 | Finer early structure, repeated opening across all sites; later coarsening remains. |

Events use sampled B+M below .10 / above .24 with hysteresis. Reformation requires
a previously observed opening. These are concentration events, not organisms or
mechanical cavities. The historical 128² report contains two extra reformation
counts from its old observer initializer; the review replays native captures
with the corrected observer and preserves that original report unchanged.

The 64² time refinement changes final mean permeability .4602→.4549 and spent
physical material 1041.52→1041.11. Spatial refinement gives permeability .4951
and spent982.20, so a converged morphology or quantitative continuum limit has
not been established. At 128² the x spectrum peaks around wavelength6 at t16–64;
power beyond .7 of Nyquist is .00078/.00482/.007 at t16/32/64. The 64² run has
more high-frequency power. Both eventually coarsen toward the box scale by t160
despite substantial reserve remaining. The useful art window is a changing
transient, not a demonstrated stationary ecosystem.

The actual 128² recording is
`outputs/living-tissue-20261015-regrowth128`: `sections.png`, standalone
`playback.html`, and `activity-film.mp4` (101 recorded frames, 1216×688,10.1s).
The film compares material and signed excitation without synthetic intermediate
simulation frames. Browser interaction QA was unavailable because the shared
browser profile was occupied; generated JavaScript passed syntax checking.
Fixed-range U8 preview fields may clip B above1.2; native Float32 fields remain
the quantitative source.

## Native-volume gate and optical contract

`experiments/art/living-tissue-volume128.toml` requests native128³, physical
extent24³, dt.02, t0–80 and 41 captures every2timeunits. It uses the same physical
coefficients and initializer construction, with a genuinely 3D random spectrum.
This is a new specimen, not an extrusion or reconstructed anatomy from the sheet.
GPU execution requires the coordinated host queue. The gate first compares a
dense randomized16³ CPU/device trajectory, including zero-fuel/reserve donors,
then requires exact same-device restart and continuation before any volume bake.

The persistent simulation arrays cost64B/site (134217728bytes at128³). The
runner preflights resident/capture budgets; the requested artifact budget is
3000000000bytes including native captures, restart, generic record and event maps.
The host process is separately bounded to300s,8GiB,no swap,CPU10/11,nice15. JIT,
captures, host transfers and filesystem writes belong to the bake wall time;
they must not be presented as an isolated kernel throughput result.

Each `frames/NNNNNNN/volume.json` declares dimensions, physical extent, time,
Float32 dtype and Julia/Fortran array order. Files `precursor.f32`, `biomass.f32`,
`matrix.f32`, `spent.f32`, `excitation.f32`, `inhibition.f32`, `memory.f32` and
`reserve.f32` contain native fields. Optical phases can map B to soft body, M to
persistent boundaries, signed x to phase color, and N to surrounding resource;
q/R are available independently. Opacity is a rendering choice, not a recorded
mechanical surface. A 2D section of the volume must be identified as a section.
No native-volume result is claimed until that separate gate and bake complete.

## Checks

The focused suite `test/living_tissue_turnover.jl` passes41 tests covering donor
funding and conservative transfers, hopping, zero-fuel action suppression,
2D/3D cell-average refinement, physical inventories, process composition,
generic recording and exact restart. Separate small CPU checks exercised the
4³ volume writer and dense16³ parity-harness plumbing; CPU-vs-CPU identity is
not an accelerator parity result. Original tissue source pins and trajectories
remain separate from this law.
