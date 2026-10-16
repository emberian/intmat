# Learning a material choice from a vanished cue

2026-10-16. This is one deliberately engineered skill experiment, not autonomous
lifetime learning or a claim that the existing living garden already uses memory
well. [Driver](../../studies/native_history_construction.jl) and
[13-check mechanism test](../../test/native_history_construction.jl) reuse the
native paid construction/chemistry and the shared paid actor–critic transaction.
The [validated compact report](../figures/history-dependent-construction.json)
contains two independent exploration-seed runs and every frozen probe; the
[compactor](../../scripts/compact_history_construction.py) checks exact paired
feature-norm equality, construction totals and conservation before packaging.

## What was made learnable

The living experiment mixes reserve management, growth, investment, allocation
and movement. The lag map shows that reserve has strong, nonlinear authority,
investment often has little funded authority, and allocation has a delayed
context-dependent consequence. We isolate allocation rather than interpreting
another combined return as a learning result.

Four fixed, nondividing residents inhabit an 8³ world. The improved material's
raw diffusion 0.006 and catalyst decay 0.05 remain. Contexts contain raw A/B
concentrations `(0.8,0.08)` or `(0.08,0.8)`, equal in total; there is no ongoing
external feed. The actor sees actual local raw balance until time 0.4, then its
own observation channel is explicitly replaced with zero. Native sensing and
the physical raw fields remain unchanged. A pair of filters (times 0.4 and 2)
retains the cue. At time 1 the actor makes one Gaussian allocation decision,
held during a 0.4-time construction pulse; its return is the change in local
battery plus usable stored precursor from decision time to time 8.

This sensor occlusion is experimental apparatus, not a claim of a new chemical
law. The actor receives no time or context label. Its two inputs are bias and
`contrast/(0.02+abs(contrast))`, where contrast is half the fast-minus-slow local
raw-balance filter difference. The unsigned control replaces that second feature
with its absolute value. Thus every paired decision has exactly the same feature
norm; the deliberate difference is access to the remembered sign. Neither actor
reads the native hidden state, which could otherwise offer another route to the
cue. The critic in every arm uses the same pre-action signed history and bias;
its frozen prediction is a noise-independent baseline, not an extra actor input.

Reserve stays at 0.2 in all arms. We explicitly set allocation exponent to one
and investment to 0.4 for this task, so allocation changes catalyst composition
at constant total request. The original exponent-two law would increase total
construction when allocation became more extreme. This change uses the existing
native request helper; precursor affordability, work, catalysis, sorption and
transport still run normally. Actual total funded catalyst in the authority
probe was 0.12287994 across allocation logits −1, 0 and +1, differing only by
less than 6e−9 in representable endpoints.

The [opportunity was checked](../figures/history-construction-authority.json)
before training. In the A-rich context, mean return
was −0.01269/−0.02100/−0.03212 for logits −1/0/+1; the B-rich context gave
−0.02650/−0.01635/−0.00953. More catalyst on the scarce channel was helpful here:
the abundant channel already obtained diminishing conversion returns. This is
measured allocation economics, not an imposed reward for predicting a cue.
No-construction controls returned about −0.09358 in either context, so the
allocation task improves an already useful paid construction pulse.

## Training and causal controls

Each arm receives 128 balanced, shuffled physical episodes. The episodes restore
fresh initial physical state; only the learned readout and critic carry across
them. This reset and reinstallation are external training apparatus with fresh
physical resources, not a free ecological reset or inherited living population.
Each actual actor/critic change pays the native battery/work transaction at
turnover price 0.001. Actor rate is 0.5 for this two-feature, one-decision task,
innovation standard deviation 0.5, coefficient bounds ±2, and per-coefficient
update clip ±0.01. The first two episodes calibrate the critic without actor
updates. The rate is not a recommendation for the 23-feature living learner.

Correct signed history, correct unsigned history, no actor credit and reversed
actor credit share initial state, cue sequence and Gaussian innovations. Both
history features use identical state machinery. Construction expense is matched;
learning expenses use the same price law and are reported, not artificially
forced equal. Frozen tests use no exploration or updates. Swapping only the
remembered sign leaves physical initial state, current actor observation,
feature norm and learned coefficients unchanged.

For exploration seed 9101, all four signed-history coefficients learned the
useful sign, ending between 0.325 and 0.390. Compared with matched zero-readout
construction at time 1, mean per-resident return improvements were:

| Readout | A-rich training context | B-rich training context | Held-out amplitudes, mean |
|---|---:|---:|---:|
| Learned signed history | +0.003810 | +0.002708 | +0.002194 |
| Learned unsigned history | +0.000601 | −0.000459 | +0.000055 |
| No actor credit | 0 | 0 | 0 |
| Reversed credit | −0.004283 | −0.003283 | −0.002509 |

Held-out concentrations are `(0.6,0.12)` and its reverse. Swapping the signed
learner's remembered cue changes its training-context improvements to −0.003213
and −0.003509. The unsigned control is exactly unchanged by the swap. Extending
the cue gap by delaying the action to 1.4 preserves benefit against the zero
readout tested at that same time; comparing against the earlier-time baseline
would confound changed pulse timing. Total training work was 0.003967 for signed
history, 0.004196 unsigned, 0.000744 no-credit and 0.004282 reversed. The signed
benefit is therefore not explained by a larger feature norm or greater learning
expenditure than the unsigned control. Maximum observed conservation residual
over this training/probe set was below 5e−14.

These are small absolute utility gains on top of a useful paid construction
pulse. They establish controllable learned history-dependent allocation in this
task, not a general memory advantage across ecological regimes.

An independent exploration seed 9203 repeated the same 128-episode protocol.
Signed-history mean gains were +0.003282 on the two training contexts and
+0.002207 on held-out amplitudes, versus −0.000020 and −0.000015 for the unsigned
control. With the longer gap, matched-time gains remained +0.003217 and
+0.002146. Reversed credit again hurt. Both runs share founder bank 211: these
replicate exploration noise, not founder populations. Across both complete
records, maximum residual was 5.39e−14 and total funded construction spread was
2.18e−8. Construction amounts are thus numerically matched, while their species
composition differs. Native GPU execution has not been claimed for this new
wrapper; the helper kernels expose a backend but this protocol was CPU-only.

## What transfers, and what does not

The [same-trained-readout test](../figures/history-construction-transfer.json)
with continuous, unoccluded **constant** sensing
is decisive: both filters converge to the same value, giving zero contrast.
Mean gains fall to +0.000400 and −0.000302 in the two contexts; the remaining
effect is the learned bias. Do not install this transient decoder unchanged in
the mobile garden and call it a learned niche policy.

What transfers is the architecture: a cheap local temporal basis can control
material composition through a paid, causally validated readout. Continuous
living deployment must train on actual local histories caused by movement,
depletion, replenishment and retained material, with current material features
available alongside the temporal contrast. Its test should select same-current
observation states with different histories and use equal-norm history swaps,
not merely compare a larger feature bank against a smaller one. This experiment
provides a working component and a discriminating test for that next integration;
it does not substitute the training apparatus for the intended living world.
