# Homeostasis without changing the memory clock

2026-10-16. A bounded controller prototype, not a native ecological result.
[Julia module](../../studies/homeostatic_temporal_memory.jl),
[29 passing checks](../../test/homeostatic_temporal_memory.jl), and
[compact replay](../figures/homeostatic-temporal-memory.json).

**Mechanism:** slowly center and scale the sensory drive to two contractive
memory modes, while leaving their recurrence coefficients fixed. This lets a
frozen readout continue using recent variation when material observations change
offset and amplitude. It deliberately changes feature gain; it does not claim
more memory capacity from a larger feature norm.

## Why this mechanism rather than a criticality target

The native `dc_controller!` applies tanh to inherited weights plus marks;
`dc_regulate!` adds adaptation feedback, leaky hidden-state updates and local
activity-dependent mark changes. Its recurrent inputs also include neighbors.
The construction score learner trains a separate readout, not these recurrent
weights. Scaling all native inputs or rows would therefore change saturation,
memory, neighboring influence and action gain together.

Three primary texts were inspected directly, not inferred from citation strings:

* Schubert & Gros (2021), [Local Homeostatic Regulation of the Spectral Radius
  of Echo-State Networks](https://www.frontiersin.org/journals/computational-neuroscience/articles/10.3389/fncom.2021.587721/full),
  DOI 10.3389/fncom.2021.587721; sections 2.1–2.6, 5.5. Flow control separates
  recurrent from external drive and uses local squared activities to adapt
  recurrent scaling. Its random-matrix approximation assumes large networks
  with suitable weight statistics and weak correlations. The authors find
  systematic target-radius errors under correlated binary input. Their delayed
  XOR result is evidence for that architecture, not a guarantee for our eight
  recurrent units. Our proposal borrows the separation of drive and recurrence,
  not their radius-regulation rule or a target of one.
* Remme & Wadman (2012), [Homeostatic Scaling of Excitability in Recurrent Neural
  Networks](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002494),
  DOI 10.1371/journal.pcbi.1002494; Results, Figures 4 and 8. In their excitatory/
  inhibitory model, relative homeostatic adaptation speeds determine whether
  adaptation stabilizes or destabilizes the network. Controlling mean activity
  does not prevent instability when recurrent excitation exceeds the available
  feedback inhibition. Our native units do not have those defined cell classes;
  the transferable lesson is to analyze the adaptation loop, not import their
  numerical timescale ratios.
* Zierenberg, Wilting & Priesemann, [Homeostatic plasticity and external input
  shape neural network dynamics](https://arxiv.org/pdf/1807.01479), 2018 arXiv
  primary manuscript, sections IV and equations 8–10. In their spiking branching
  model, input strength changes the operating regime even under a fixed target
  firing rate. When the intrinsic correlation time approaches the effective
  homeostatic timescale, the quasi-stationary argument fails and bursting can
  appear. This supports keeping feedback adaptation slower than the retained
  signal as a design hypothesis; it does not establish that our material should
  burst, or that criticality maximizes its usefulness.

## Local rule and its limited guarantee

At each .2 decision tick, use the previously installed center `b` and variance
`v` to compute `d=tanh((u−b)/sqrt(max(v, .01²)))`. Then update the two modes:

```
z_fast ← r_fast z_fast + (1−r_fast)d     r_fast = exp(−.2/.4)
z_slow ← r_slow z_slow + (1−r_slow)d     r_slow = exp(−.2/2)
```

The candidate homeostat uses exponential central-moment updates, with
`a=1−exp(−.2/20)` and `δ=u−b`:

```
b_new = b + a δ
v_new = (1−a)(v + a δ²)
```

It observes only the current local signal and its retained statistics. No phase,
future value or prediction label enters this update. In the prototype, installing
the proposal costs `1e−4 (|Δb| + |Δlog(gain)|)` from a toy account initialized to
one; if insufficient, neither statistic changes. Memory continues running.
Variance-estimator bookkeeping and filter state updates have no additional
maintenance charge. This is a declared parameter-turnover model, **not native
energy accounting**, and its units must not be added to native chemistry ledgers.

For two copies receiving identical inputs with identical homeostat state, a
perturbation to a memory mode contracts by its own `r` every step, even while
the homeostat adapts. The tests verify that equality, rather than estimating a
Lyapunov exponent. It is a conditional open-loop result: differing homeostat
states add their own slow modes, and action-mediated changes to material input
can create closed-loop dynamics not covered by this guarantee. Existing memory
is not rescaled or erased when the sensory gain changes.

## What it enables in the bounded replay

The input is a period-8 bounded material-like fluctuation. Three consecutive
160-time regimes use `(offset, amplitude) = (0,.25), (.6,.025), (0,.25)`.
Each controller fits an equal-size linear readout on only the last 80 time of
the first regime, forecasting relative phase one time ahead. Readouts then stay
frozen. Fitting is an offline calibration, not a paid reward-learning claim.

The fixed control keeps the original center/scale, runs the same two memories,
and receives exactly the adaptive arm's work debits. It has its own equally
calibrated readout, so it is not penalized by using another representation's
coefficients. A zero-budget homeostat reproduces the fixed predictions exactly.

| Controller | Initial late MSE | Changed-regime late MSE | Restored-regime late MSE | Toy work |
|---|---:|---:|---:|---:|
| Slow homeostasis, time 20 | .00438 | .10677 | .00435 | .001451 |
| Fixed statistics, yoked work | .00413 | 1.24589 | .00413 | .001451 |
| No update funding | .00413 | 1.24589 | .00413 | 0 |
| Fast homeostasis, time .2 | .01269 | .14370 | .01269 | .026322 |

“Late” is the entire final 80 time of each regime, not a selected good cycle.
Fast homeostasis adapts to the shift sooner but writes much more, and is less
accurate late. It is not uniformly worse: its changed-regime first-40-time MSE
is .141 versus .614 for slow homeostasis. Slow adaptation also temporarily
hurts when the original regime returns (.452 first-40-time MSE versus .224 for
fixed statistics). These are actual response-speed/retention tradeoffs, not
evidence for a universally optimal clock.

At times 300 and 304, the instantaneous input is .6 in both cases, while its
recent directions differ. Desired forecasts are −.707/+ .707; the adaptive
controller gives −.559/+ .565. The fixed controller gives −.864/−.864. Replacing
only the retained fast/slow features swaps the adaptive forecast; erasing them
leaves approximately zero. This history intervention does not match feature
norm and is only a dependence check. The main capability is maintaining a
calibrated temporal response under changed sensory scale, not winning a
norm-matched memory comparison.

The entire replay is 2,400 scalar observations per arm and uses no GPU or
material simulation. The synthetic rhythm is predictable and the readout target
is relative phase; performance on irregular material feedback or absolute
resource prediction has not been established.

## Native integration contract, if the living task needs it

Keep unnormalized raw balance/total/catalyst observations available beside the
new temporal features: centering must not hide physically meaningful scarcity.
Start with only the balance channel. With Float32 structure-of-arrays storage,
the full three-channel version adds six statistics (24 bytes/resident) to the
existing six temporal values, for 52 bytes including the identity integer.
The present scalar prototype uses Float64 and is not a device kernel.

At an actual decision, form the normalized drive from old statistics, update
the modes, store the resulting actor features in the existing cohort queue,
and separately transact any funded homeostat proposal. Never re-normalize a
queued feature when its reward matures. On a new resident, initialize the center
from its first actual observation, reset the modes and use a declared variance
prior; do not transfer another resident's credit queue. A production study must
charge proposals through native energy/ledger helpers and checkpoint the added
state before it can claim paid continuation.

The variance floor bounds gain at 100, but can still amplify sensor noise.
Pick that floor from measured sensing variation before deployment. Abrupt
offset changes are initially indistinguishable from useful events, so this
component cannot remove them immediately without also deleting relevant
history. The useful next gate is a short recorded living-input replay with a
frozen allocation readout and the same work-yoked control—not another broad
radius or timescale search.
