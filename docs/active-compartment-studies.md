# Funded active boundary study, 2026-10-14

A bounded five-run comparison tested whether stronger paid neural feedback changes
compartment topology rather than merely painting the existing phase boundaries.
Run `julia --project scripts/compartment_sweep.jl OUTPUT_DIRECTORY`.
The completed data live in `outputs/compartment-sweep-20261014`: `summary.json`,
`both-phase-topology.json`, per-condition metrics, exact state recordings at
`t=0,80,160`, replay JSON at `t=80,160`, and `phase-comparison.png`.

Every run used seed 11, 48×48 cells across a 24×24 periodic physical domain,
`dt=.1`, constant phase mobility, conserved mean composition −.12, identical
24-mode small initial perturbations, 96 founders, capacity 700, lifespan 80,
uptake .7, chemostat supply .023 and sink .01. The fair-share resource allocator
and catabolic fuel model used catabolic fraction .5 and energy yield 2. Neural
building, signaling and remodeling remain funded actions. Source hashes are
stored with every result; no simulation source was changed during this study.

| Condition | Remodeling gain | Signal bias | Membrane pinning | Neural remodeling rate | Secretion rate | Positive / negative components at t160 | Interface length t160 | Population |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Default feedback | .3 | .1 | .5 | .03 | .025 | 1 / 1 | 106 | 699 |
| Moderate feedback | 1 | .5 | 2 | .15 | .08 | 1 / 1 | 110 | 700 |
| Strong feedback | 4 | .5 | 2 | .15 | .08 | 2 / 1 | 110 | 692 |
| Strong parameters, neural feedback disabled | 4 | .5 | 2 | 0 | 0 | 1 / 1 | 108 | 700 |
| Strong parameters, no inhabitants | 4 | .5 | 2 | — | — | 1 / 1 | 108 | 0 |

Component counts use periodic nearest-neighbor connectivity on both signs of the
composition field. A single connected positive matrix can surround multiple
negative pores; positive connectivity alone is not a compartment count. Here
both signs were measured explicitly. Interface length counts sign-changing cell
faces times physical spacing, a resolution-dependent proxy rather than a smooth
geometric contour estimator.

The strong run maintained two substantial positive domains from t20 through t160.
At t160 the largest contained 53.95% of the positive cells, so the second component
was not a tiny isolated pixel. The negative phase remained connected. The two
controls produced identical final phase fields, while strong active feedback
changed the topology. This supports an effect of funded neural actions on boundary
evolution within this one seed and observation window. Multiple changed rates
were tested together; the experiment does not isolate remodeling from secretion,
signal bias or membrane pinning.

It does **not** show permanent arrest of coarsening. Strong interface length still
fell from 173 at t20 to 110 at t160. Its late interval t140–160 changed the sign of
111 cells versus 96 in the controls, indicating ongoing reorganization but not
sustained fragmentation. Phase variance was .525 for strong versus .636 in the
controls, consistent with membrane pinning favoring local mixing. Whether the two
domains persist for several more lifetimes needs a longer run and more seeds.

The living runs had 3,853–4,159 births and 3,257–3,556 deaths, substantial turnover
rather than a frozen population. Final material-budget residual magnitude stayed
below .00017 resource units; fuel-ledger residuals were zero in the reported
measurement. The entire sweep took about 104 seconds of simulation execution,
including first-run compilation overhead, plus snapshot preparation.

`experiments/arena/active-boundaries.toml` selects the strong settings for an
exploratory live arena. It uses the normal arena initializer: twelve larger
modes, max wave number four, mean zero. That initializer differs from this matched
study, so the profile is a candidate for exploration, not a reproduction of the
specific two-domain outcome. The study script is the reproducible comparison.
