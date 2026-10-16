# Joint reserve and material policy

The two-action pilot makes a concrete constraint visible: a secretion command
can lose all physical authority when native reserve-demand catabolism consumes
precursor first. A fixed lower-reserve founder restores paid construction but
hardwires the tradeoff. `studies/native_joint_construction_credit.jl` instead lets
a paid local policy choose reserve, investment and A/B allocation together.

It wraps the same construction sidecar and delegates its source schedule,
reward transaction, identity state and physical hooks. The existing two held
logits command sigmoid investment/allocation. A third held Gaussian logit
commands `actions[i,11] = tanh(u₃)` after native regulation and before native
metabolism. Thus the actual reserve target is

    energy_capacity × current_mass × (tanh(u₃) + 1) / 2.

Its initial mean `atanh(−0.6)` recovers reserve fraction 0.2. The new policy does
not directly add, remove or reserve physical inventory. Native metabolism still
limits oxidation by available stock, rate and demand, posts heat/waste/work,
pays maintenance, and only then offers remaining resources to habitat building.
Persistent reserve-policy coefficients pay the same turnover cost as the other
action coefficients. No extra resource budget is attached to the extra action.

Every policy uses the same feature-bank shape: bias; three native local material
observations; energy divided by battery capacity; two precursor stores divided
by store capacity; and the existing hidden and adaptation vectors. With h=8 this
is 23 features × 3 actions. The readout-memory ablation zeros the last 16 input
features, with the same array dimensions, Gaussian innovations, physical neural
costs and recurrent dynamics. It tests access from the learning readout to
existing recurrence; it does not remove recurrence from the native uptake and
other cellular actions. The native controller's normal local observations,
including resident age, still feed that recurrence. Neither variant receives
external feed phase, wall time, future supply or an oracle context label.

Decisions remain every 0.2 physical time, rewards every 2, actor rate 0.003,
noise scale 0.3, coefficient bounds ±0.3, increment bounds ±0.01, and turnover
price 0.001. The local eight-feature critic and first-window calibration remain
identical to the two-action study. Extra score dimensions can increase update
variance and clipping; these diagnostics must accompany any behavioral claim.
Fixed controls include matched noise/no-credit, reversed actor credit,
no-secretion and no-raw supply. Frozen learned/erased probes share physical state
and paid history, with continued/shifted future feeding. Reports distinguish
founder seeds from exploration seeds.

Implementation shares the construction run protocol through a sidecar factory;
the exact policy and protocol drivers are archived at run start. An execution
launcher can project the same kernels across independent worlds without
forking their equations. This is an experiment in temporal material control,
not yet evidence of a memory advantage or cooperation.

The [four-run comparison](../figures/native-joint-construction-comparison.json)
holds founder bank 211 fixed and pairs exploration seeds 6201/6277 with recurrent
readout enabled/ablated. All five arms and all frozen probes preserve the four
original identities. Correct-credit terminal mean U exceeds its matched
no-credit arm by 0.0300/0.0181 with recurrent readout and 0.00250/0.00115 without;
reversed credit is below no-credit in all four runs. These training outcomes
include the different material each policy built. Same-state frozen probes
remain small and phase sensitive; no-credit learned/erased differences are
exactly zero. They test a different question from accumulated training return.

This is not yet a temporal-memory advantage. Added recurrent features also
increase the effective readout update gain and the number of paid coefficients.
For seed 6201 with recurrent readout, actual final frozen investment logits move
at most 0.0332 from their zero-mark baseline, allocation at most 0.0261, but
reserve moves as much as 0.1778. Against exploration sigma 0.3 these are 0.11,
0.087 and 0.59 sigma respectively. There is no fed-arm coefficient clipping, and
total actor+critic turnover work is 0.00176 over four residents and 24 time.
The small frozen return difference therefore cannot be diagnosed simply as
uniformly negligible policy change. Reserve movement, differing feature gain,
material capital and causal lag must be separated before assigning a learning
mechanism explanation.

The next bounded intervention is the three-logit marginal lag map described in
[Learning mechanics as a material design lever](learning-mechanics-design-levers.md).
It preserves the original continuing law. Any later decision-aligned credit rule
will retain its own pre-action baseline; carrying old scores while subtracting a
new predictor conditioned on their consequences would change the credit meaning.

The [joint CPU/AMD gate](../figures/native-joint-construction-gpu.json) passes four
six-time cases: feed on/off and dt 0.02/0.01. It exercises the original host hooks
through the execution launcher, paid updates, held reserve/material commands,
and strict checkpoint restoration with exchanged masks. Maximum checked
observable/coefficient discrepancy is 1.22e−7; maximum physical residual is
1.51e−14. Eleven focused tests cover local feature ordering, recurrent-feature
ablation, held samples, identity clearing and the reserve transform. This proves
the native implementation boundary, not a behavioral advantage.
