# Memory representation on actual living-material signals

2026-10-16. **Keep the homeostatic prototype separate for now.** Its synthetic
scale-shift benefit did not transfer to this early-trained readout and native
observation window. This is a bounded representation result, not a general
verdict on homeostasis or on reward-trained temporal allocation.

The [recorder](../../studies/record_temporal_material.jl),
[analysis](../../studies/material_memory_transfer.jl),
[13 focused checks](../../test/material_memory_transfer.jl), and
[compact evidence/replay](../figures/material-memory-transfer.json) are separate
from the continuing learning experiment.

## The physical observation epoch

One CPU16 run used the funded slow-growth phenotype in
`experiments/native/living-temporal-construction-16.toml`: growth .1, reserve
.02, fixed-total allocation exponent1, source-mask persistence8, observation
cadence .2, action cadence1. Allocation stayed at its zero learned logit;
exploration, critic and actor learning were disabled. Native motion, recurrent
physiology, activity-dependent marks, chemistry, mutation and division continued.
The recipe's finite source/sink budgets and original physical laws were retained.

The independent `outputs/material-memory-epoch/` snapshot contains 62 files:
Project/Manifest, all Julia source, the complete temporal-policy include chain,
recorder and recipe. The run used that snapshot as its Julia project, rather
than reading concurrent study edits. Every file hash is in the compact report.
Full local observations remain at
`outputs/material-memory-observations/observations.json`, with its SHA256 also
in the report. No old run was resumed or altered.

The wrapper records immediately after the real sensing/policy-observation hook,
before that step's controller and metabolism. It does not infer the observation
clock from `advance!`'s later endpoint. Samples occur at steps divisible by10
plus each resident's first actual observation, including off-cadence births.
IDs are native identities, not reusable storage slots. First observations have
native fast=slow=current in all channels; all 4,746 ID/step pairs are unique.

The run reaches32 with 24→41 residents and17 recorded births. The final material
and potential residuals are −2.84e−14 and3.98e−13, respectively. Capturing arrays
does not change actions, physical state, private filters or birth handling.

## Causal calibration and held-out clocks

The target is each resident's **own** raw-balance, raw-total and catalyst-total
observation two physical time later. These are bounded native sensor values,
not concentrations, utility, or preferred actions. Pairs never cross identity;
missing/dead future identities are not fabricated.

Calibration uses inputs at `[2,14)` and targets at `[4,16)`: 1,466 pairs.
The frozen readout is evaluated on inputs `[20,30)` and targets `[22,32)`:
1,760 pairs. Unlabelled observations between these windows still update memory
causally. No source clock, future target, source mask or lineage metadata enters
the readouts. Source-switch labels are used only to group evaluation errors.

Four representations share one declared ridge penalty .01 and training-only
column centering/scaling. Each fits its own equal procedure once; no held-out
normalization or hyperparameter selection occurs. Raw-only uses three current
channels. Raw+fixed adds the six actual native fast/slow states. Raw+homeostatic
instead adds six normalized states from the earlier prototype. A homeostatic-only
ablation checks the cost of discarding absolute observations. The two augmented
banks have the same nine input dimensions; the raw-only baseline is smaller.

The homeostat retains its declared time20, .4/2 memory times, variance prior
`.25²/2` and gain floor .01. On a first ID, its center is the first actual
observation and normalized memories are zero. It uses the actual elapsed time
between samples, including first off-cadence newborn intervals. It sees no parent
history or future labels. This is offline representation evaluation: its toy
turnover is reported separately and is not charged to the native material run.

## Actual signal scale and transfer

Across successive eight-time phases, balance standard deviation is
.02677/.02496/.01882/.02024. Its median remains near zero. There is no synthetic
tenfold amplitude collapse. Raw-total medians are
.00997/.01421/.01364/.00951, while catalyst-total medians are
.8705/.5190/.4659/.5050. Those changing catalyst levels describe actual material
creation/retention, not an established nuisance offset to erase.

Held-out mean squared sensor errors:

| Representation | Balance | Raw total | Catalyst total |
|---|---:|---:|---:|
| Current-value persistence | .0002265 | .00006124 | .002491 |
| Raw-only fitted readout | .0002231 | .00007960 | .003201 |
| Raw + fixed temporal | .0002425 | .00007950 | .004070 |
| Raw + homeostatic temporal | .0004126 | .00033062 | .009060 |
| Homeostatic temporal only | .0005013 | .00086399 | .012228 |

The normalized representation fails both across source-switch targets and within
a source phase. Of the1,760 test pairs,330 cross the switch at24; balance error
there is .000471 for raw+homeostatic versus .000357 raw-only and .000345
raw+fixed. This does not establish a meaningful temporal-control advantage for
the small raw+fixed difference; it simply prevents attributing the normalization
failure only to within-phase samples.

The strongest failure occurs on the410 test pairs from IDs absent in training.
Catalyst error is .02314 with raw+homeostatic versus .00270 raw-only. On the1,350
seen-ID pairs, the corresponding errors are .00478 and .00335. Resetting newborn
statistics correctly creates different birth-relative coordinates from older
residents. Combined with evolving material levels, this is a plausible source
of covariate shift; it has not been isolated as the unique cause. Final
homeostatic gains span3.19–12.52, far below the gain cap100, so a saturated gain
ceiling is not the explanation here.

Representative actual histories for IDs1,13,25 and35, their homeostat states,
and frozen forecasts are included in the compact replay. The tests separately
verify that appending future observations leaves past encodings identical,
newborns do not inherit another ID's history, duplicate/backward clocks fail,
and all target pairs obey the declared identity/time boundaries.

## Consequence for the controller design

Do not install a paid native homeostat on this evidence. Preserve absolute
material observations; normalized-only features substantially degrade scarcity
and catalyst forecasts, and keeping raw values alongside them is necessary but
not sufficient to make this fitted representation transfer.

This forecast task is not the allocation objective. A temporal feature can help
choose an intervention even when it barely improves prediction of the unperturbed
sensor level. Conversely, improving that prediction would not prove better
construction. The next useful lever belongs at the observed conditional action
effect and its reward estimator, with growth/partition and current precursor
authority kept explicit—not at another unsupervised gain or timescale sweep.
