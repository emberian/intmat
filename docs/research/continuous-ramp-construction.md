# Construction that anticipates a real resource ramp

This experiment removes policy-observation occlusion from the earlier
[history-dependent construction component](history-dependent-construction.md).
Real finite raw supply creates a transient local cue. Recent change predicts
which substrate will continue arriving while a paid catalyst deposit acts.
The study is still externally initialized episodic training, not autonomous
learning in the moving/dividing garden.

Implementation: [native_ramp_construction.jl](../../studies/native_ramp_construction.jl),
using the native physics and paid update transaction through the existing
history study wrapper. No production or living-scene laws change.

## A physical crossing with hidden direction, not hidden observations

Four immobile, nongrowing residents inhabit the same 8³, unit-volume native
world as the component assay. One raw species starts at concentration .0656,
the other at .8. A real source supplies the initially scarce species at .72
mass/time, with finite budget 5.76. The opposite context swaps the species.
The source acts before sensing, so decision step 50 observes 51 grants of .02:
`.0656 + .72 × 1.02 = .8`.

Native raw-balance observations at that decision are **exactly zero** for all
four residents in both contexts; raw-total observations also match. The
continuously observed fast-minus-slow contrast is ±.09201148. No sensor is
masked, no physical state is overwritten, and no time or source-direction bit
is passed to the actor. Future supply is an explicit physical boundary
condition whose persistence can be inferred from the past ramp.

The actor chooses one Gaussian allocation logit at time 1, held for .4. Total
investment is .4, reserve fraction .2, allocation exponent 1; thus allocation
changes composition while funded construction stays constant. Return is the
resident's usable-inventory change from decision to time 8. Catalyst decay .05
gives a 20-time lifetime; diffusion .006 gives a nominal exchange-radius
diffusion time `.18²/(6×.006)=.9`. The fast/slow filter times .4/2 straddle the
one-time cue development, while supply persists across the material payoff.

This is a different operating regime from static stock. In the
[forced-action results](../figures/continuous-ramp-authority.json), supplying A
gives mean returns .4800/.6523/.6784 at allocation logits −1/0/+1; supplying B
gives .7019/.6725/.4864. Catalyst for the **forthcoming supplied channel** wins.
Persistent inflow can keep its catalyst productive; it does not merely act on
the instantaneous equal stock. The earlier static assay instead favored the
scarce channel, so its learned sign should not be transplanted here.

When the same source exhausts at the crossing, with identical observed balance
and history, the opposite preference disappears: the balanced allocation
returns .1150/.1129 in the two contexts, and both ±1 deviations do worse.
This control deliberately removes future material/energy supply (budget .7344
instead of 5.76); it establishes the predictive opportunity, not equal-budget
performance superiority.

## Exactly matched actor gain

Let `x=tanh(raw_balance/.2)` and `y=sqrt(max(0,1−x²))`. The history actor reads
`[1,x,sign_plus(trend)*y]`; the instantaneous actor reads `[1,x,y]`.
`sign_plus(0)=+1` is declared explicitly. Both vectors have squared norm 2,
including at zero trend, and both expose the same instantaneous balance.
The instantaneous representation is a deterministic function of current input;
it receives no history magnitude. Both arms run the same filter arithmetic.

At the crossing these become `[1,0,±1]` and `[1,0,1]`. Both have three paid
coefficients and identical exploration, actor bounds, clipped update rule,
critic basis and physical budget. Frozen history-sign swaps preserve the
entire physical episode and feature norm while changing only the history
coordinate. They leave the instantaneous controller exactly unchanged.
This is an information comparison against this specified instantaneous
readout, not a claim to outperform the existing 23-feature recurrent policy.

## Predictor range is a learning design parameter

The first 64-episode run retained the old predictor's unit-scale inputs and
coefficient bound ±.3. Its bias hit +.3 in every resident and every arm, although
the paired mean returns are about .65. A bias plus opposite signed feature
cannot fit both contexts above that bound. The residual offset interacts with
the actor's .01 coordinate clipping: the update can largely reflect innovation
sign instead of modest action-quality differences. The history arm failed to
improve (mean −.00539 versus zero readout), and is preserved in the
[saturated-predictor record](../figures/continuous-ramp-saturated-critic.json).

The targeted correction scales **only the shared predictor basis** by 4,
allowing a bias prediction up to 1.2 under the unchanged coefficient bound.
Actor features, rate, exploration, physical supply and reward remain unchanged.
Normalized predictor fitting and all installed coefficient changes still use
the native paid transaction. Because work prices coefficient displacement,
this parameterization also changes predictor work per unit prediction change;
costs should be compared within the corrected matched experiment, not presented
as a free efficiency improvement over the old parameterization.

The [corrected 64-episode run](../figures/continuous-ramp-construction.json)
used exploration seed 9301, sigma .5, actor rate .5, and two initial
predictor-calibration episodes. The four arms share balanced shuffled contexts
and Gaussian innovations. Frozen evaluations disable exploration and learning.
Mean gains relative to the same zero readout are:

| Controller | Training ramps | Slower held-out ramps | Training ramps, history swapped |
|---|---:|---:|---:|
| History | +.018351 | +.002808 | −.034043 |
| Instantaneous | −.001700 | −.000933 | −.001700 |
| Reversed credit | −.039018 | −.013654 | +.020004 |

The held-out ramp supplies .5 rather than .72 and starts at .29 rather than
.0656, preserving the same crossing concentration and observation clock.
All four correct-credit history coefficients are positive (.265–.342);
the predictor remains inside its coefficient bound (largest absolute
coefficient .170). Total training learning work is .005988 for history versus
.006177 for the instantaneous control and .006336 for reversed credit.
Funded construction differs by at most 1.86e−8 across these training/probe
records; maximum conservation residual is 7.82e−14.

An [independent exploration seed 9403](../figures/continuous-ramp-construction-noise9403.json)
repeated the same 64-episode protocol: history gain +.015528, instantaneous
−.000916, reversed credit −.028758, and history swap −.025502. The slower-ramp
history gain was +.002745. History work .005717 again stayed below the matched
instantaneous arm's .005800. Both runs share founder bank 211; this repeats
exploration noise, not founder populations.

The [77 focused feature tests](../../test/native_ramp_construction.jl) check
constant norm, zero-trend convention, matching instantaneous coordinates and
sign-swap invariance of the control. The
[report compactor](../../scripts/compact_ramp_construction.py) additionally
checks the actual matched training observations, actor norms, construction
amounts and physical closure. The mechanism was exercised on CPU only.

To reproduce, include the study in the verified Julia project environment and
call `NativeRampConstruction.run(profile, fresh_output; noise_seed=9301)` or
`noise_seed=9403`, where `profile` is
`experiments/native/living-patch-recruitment-16.toml`. The preserved failed
parameterization is `critic_scale=1.`. `authority(profile,fresh_output)` runs
the twelve forced-action future-supply checks. Outputs require fresh paths;
the study saves its executed driver and never migrates an existing run.

## Transfer boundary

The ramp source density .72 is 14.4 times the living recipe's volume-average
.4/8=.05 per species. The garden's alternating four-time source masks, movement,
and resource depletion can create real trends, but their predictive persistence
and action authority need measurement at the actual locations. At late times,
large accumulated raw stocks can make a small derivative irrelevant.

The useful transferable mechanism is a local temporal sign combined with
current composition, a predictable replenishment interval longer than the
construction response, and a predictor range that covers the attainable
returns. A continuously acting garden policy must train on those actual
histories, retain multiple overlapping decision scores, and reset newborn
filter state from its first real observation. This episodic crossing experiment
does not establish those properties or descendant learning.

## Explore the mechanism

The [laboratory explorer](../../site/laboratory/continuous-ramp/index.html)
connects the physical crossing to forced allocation outcomes and frozen learned
probes. Select A/B supply, continued/exhausted source, predictor range,
exploration seed, slower ramps and a frozen history-sign swap. All allocation,
return, prediction and work outcomes come from the recorded reports; the
source-only stock reference and predictor bias ceiling are explicitly analytic.
Missing prediction values in the saturated-predictor compact report remain
unavailable. The static explorer neither advances Julia nor infers a training
trajectory.

Rebuild with `python3 scripts/build_ramp_construction_viewer.py`.
[Compact browser evidence](../../benchmarks/2026-10-16-continuous-ramp-explorer.json)
records 96 exercised control combinations, exact staged-page identity and the
reviewed desktop/mobile layouts. The route adds about 64 KB including measured
data and requires no external libraries.
