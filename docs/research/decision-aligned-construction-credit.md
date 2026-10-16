# Decision-aligned construction credit

This separate study law gives each sampled action a return measured from its own physical decision time. The preceding [lag map](construction-lag-map.md) found delayed allocation effects and a reserve intervention whose utility effect changed sign between sampled lags 4 and 8. A longer credit horizon is therefore a concrete mechanism to test, not an assumed improvement.

For each resident and decision j, the sidecar records the 23 local actor features x_j, the three Gaussian score factors ε_j/σ², usable potential U_j, eight predictor features φ_j, and prediction b_j made before sampling ε_j. At horizon H it forms

    R_j = U(t_j + H) − U(t_j)
    ΔW = η (R_j − b_j) x_j ε_jᵀ / σ².

The native paid updater bounds each coefficient increment to ±0.01 and each learned coefficient to ±0.3. Actor and predictor changes are one affordable transaction: κ=0.001 times their total absolute coefficient change, debited from the actual cell battery into native work. An unaffordable transaction changes neither set of coefficients. The queue retains factorized features and score factors, not a full feature-by-action tensor for every decision.

The actor baseline is not recomputed from observations of the action's later consequences. Predictor training uses that cohort's stored pre-action features. The archived initial prototype also used its stored prediction for the predictor residual, reproducing the shared updater’s delayed-gradient semantics. Revision 2 instead fits the current predictor coefficients against the same stored pre-action features at maturity; the actor still subtracts its original saved prediction. This removes an unnecessary predictor-gradient delay without making the actor baseline depend on later consequences. The normalized predictor gain is 0.03 per decision versus 0.3 per old two-time window. This approximately matches nominal cadence, not the exact nonlinear predictor trajectory.

## Fixed physical comparison

All laws use the same four-resident founder fixture, original finite source budgets of 14 per channel, alternating opposed source masks, native chemistry and native work ledger. Exploration runs from time 0 through 16. Cohort targets then drain through time 24 without new scores or exploratory noise. The deterministic policy can still change when pending paid updates mature. Frozen learned/erased probes start from the same physical state at time 24 and end at 28, with continued or shifted source phase. No branch restores spent source budgets.

The H=2 and H=8 cohorts each enqueue 80 decisions per surviving resident. A prediction made before any return has matured is marked uncalibrated: it trains the predictor but does not update the actor. This is a variance-control design choice, not a requirement for an unbiased score estimate: a fixed zero baseline is also independent of the sampled action. It skips 10 decisions at H=2 and 40 at H=8. The earliest actor update is at time 4 or 16, respectively; the last pending target matures at 17.8 or 23.8. These different learning exposures are explicit. A difference between horizons is not by itself a pure horizon advantage.

A separate control wraps the original window-credit law with the same exploration cutoff at 16, deterministic continuation to 24, and frozen probes to 28. Its final exploration-window update occurs at 16. The original law, source archives and earlier 24-time exploration experiments remain distinct.

Actor rate remains 0.003: ηRΣg and ΣηRg have the same unclipped aggregate when targets agree. Dividing the actor rate by ten solely because updates occur ten times as often would change that comparison. Clipping, overlapping targets, parameter drift, calibration and paid total variation still make the algorithms different. Reports retain raw proposed increment magnitude, bounded proposed/funded increments, clipping, paid actor/predictor update counts, and queue counts.

## Boundaries

Identity changes invalidate pending private credit and reset learned actor/predictor state. No descendant receives its parent's score or unpaid return. Death remains an outcome; frozen same-resident utility is omitted for invalid identities. Failure status gates all queue and payment kernels. Queue memory is study-side credit infrastructure; its allocation/read operations are not an additional thermodynamic storage model. Paid coefficient turnover and native cellular computation retain their existing physical semantics.

Usable potential is local battery plus efficiently usable stores. It does not directly value retained extracellular capital, and a horizon can still miss later material benefits. Neither this target nor a successful pilot establishes cooperation or open-ended learning. No-raw, no-secretion, matched-no-credit and reversed-credit arms remain in the protocol; frozen/erased probes test the retained policy separately from training-time inventory differences.

## Implementation

`studies/native_cohort_construction_credit.jl` reuses the construction hooks, held-drive actuators and shared paid update kernel. `studies/native_window_construction_control.jl` provides the matched old-law schedule. The shared updater only gains an optional statically selected predictor gain in its chemistry argument; absence retains the historical 0.3 arithmetic. `test/native_cohort_construction_credit.jl` checks timing, saved pre-action data, queue draining, identity invalidation, physical closure and exact explicit/default predictor-gain equivalence. Native parity and completed pilot comparisons are recorded separately when available.

## Initial CPU pilot, retained as a separate law

One founder/noise pair (211/6201) completed all five arms for H=2, H=8 and the matched original window law. All cohort queues drained: 80 matured, zero pending/invalid; actor updates were 70 per resident at H=2 and 40 at H=8, with no coefficient increment clipping in the fed credit arms. No-credit frozen differences were exactly zero.

| Law | Mean training U, credit − no credit | Frozen ΔU continued | Frozen ΔU shifted | Paid actor+predictor work, credit |
| --- | ---: | ---: | ---: | ---: |
| Original window, act16 | +0.009349 | −0.0000071 | −0.0001346 | 0.001230 |
| Initial cohort H2 | +0.007491 | +0.0000652 | −0.0002603 | 0.002389 |
| Initial cohort H8 | −0.000755 | −0.0002088 | −0.0000284 | 0.004092 |

These results do not establish an advantage for longer credit. H8 had delayed predictor gradients and fewer, later actor updates; all laws still face a phase-sensitive physical objective. Revision 2 corrects the predictor fit, not these historical outputs. The compact initial comparison is `docs/figures/native-cohort-construction-initial.json`, with exact cohort driver/protocol/updater under `studies/outputs/native-cohort-construction-initial/`. No claim of independent replication follows from this one seed pair.

## Current-critic revision and native gate

Revision 2 completed the same one-pair CPU comparison. H2 credit exceeded no-credit training U by 0.007518, with frozen continued/shifted differences +0.0000579/−0.0002595. H8's training difference was −0.000635, with frozen differences −0.0002386/+0.0000202. The corrected predictor does not establish a longer-horizon advantage. These retained-policy effects remain small and phase-sensitive; the next useful step is a living spatial substrate, not treating another sign change as a general learning result.

The focused CPU gate passed 29 checks, including opposite actor/critic directions when the stored baseline and current predictor differ. Native CPU/AMD gates covered H2 at dt 0.02 and 0.01, fed and unfed worlds, plus H8 through time 18 and its deterministic-drain transition. Queue/calibration/update counts matched exactly; maximum checked numerical difference was 2.36e−7 and closure residual 5.68e−14. See `native-cohort-construction-current-critic.json` and `native-cohort-construction-gpu.json` in `docs/figures/`. Exact revision-2 execution sources are retained separately from the initial prototype.
