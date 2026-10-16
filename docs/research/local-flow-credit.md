# Recent flow is a conditional baseline, not a better objective

The unchanged native history learner was replayed for 24 time units, with a
read-only shadow estimator. Its final observation—including every cell and
ledger entry—exactly matches the original history arm, as do paid work
(0.011672535269106082), actor updates (423), and coordinate clips (464).
The compact result and source hashes are in
[`shadow-flow-credit.json`](../figures/shadow-flow-credit.json).

For each living ID, observe inventory utility `U` every .2 time units. Add only
that parent's actual cumulative child-utility gifts, avoiding a spurious loss at
partition. Estimate its derivative with a one-time-unit exponential filter and
freeze `b = 4 * derivative` at each decision. New IDs start with zero derivative;
their first observation supplies the reference value. The shadow never modifies
the actor, critic, physical state, or existing queue. Match its stored baseline
to the original baseline, score, features, and actual matured reward by ID and
decision step. No future reward enters either pre-action baseline.

| Eligible matured windows | Count | Original residual RMS | Flow residual RMS | Original score second moment | Flow score second moment |
|---|---:|---:|---:|---:|---:|
| All training | 423 | .10973 | .14365 | .26460 | .30758 |
| Maturity time at least 12 | 327 | .09258 | .01561 | .20004 | .004872 |

The score statistic is the empirical mean squared norm of `residual * score *
features`, before gain, clipping, or payment. It is not an estimate of true
gradient variance from independent replicas. On all eligible windows, the
shadow's hypothetical coordinate clipping fraction falls from .3656 to .1395,
despite its larger second moment: a few large early errors dominate the latter.
At decision time 4, initial inventory drawdown produces predictions near −.8
for eventual returns near −.05. A derivative extrapolation is useful after that
transient, but fails across the transient itself. The late subgroup is a
diagnostic clock split, not permission to advertise the full intervention as a
variance improvement.

An earlier fixed-policy late-window replay in `local-flow-baseline.json` found
RMS .01197 versus .01323 for subtracting the preceding four-time return and
.06134 without subtraction. Its Gaussian score calculation is an analytical
fixture, not native learning evidence. The matched learner replay above
supersedes any broad inference from that easier case. Eleven helper checks cover
causality, ID initialization, partition accounting, and action-independent
baseline subtraction.

## What to build next

Do not replace the current critic with this unconditional flow estimate. First
declare the developmental objective. Inventory `U` decreases when balanced
precursors become biomass; a better baseline cannot change the expected reward
gradient or repair that preference. Native growth consumes half a unit of each
precursor and `growth_work` energy per biomass unit. Replacement accounting can
therefore use `V = U + c * mass`, where
`c = efficiency * (potential_A + potential_B)/2 + growth_work` (approximately 3
here). Substrate-only accounting gives approximately 2.6 and leaves assembly
work as a cost. These are explicit capital valuations: biomass is not stored
recoverable chemical energy in the native physics. Both choices should be
reported without fitting their coefficient to an action's preferred sign.

Conserve this chosen value through division by including child biomass in the
parent's exact partition gift. Charge later learning and maintenance normally.
Do not automatically award replacement value to extracellular catalyst: that
would reward production even where the deposit has no subsequent use. The
single-resident authority study and developmental accounting checks own that
next intervention. Revisit baseline design only on the chosen objective's
actual paired rewards, including its initial transient and births.

The recorder is `studies/shadow_flow_credit.jl`; it uses the existing native
hooks and accepts `INTMAT_TEMPORAL_SOURCE` for an immutable source epoch. Full
observations and the executed source remain in ignored local outputs. No paid
flow law, GPU run, or existing-run migration was introduced.
