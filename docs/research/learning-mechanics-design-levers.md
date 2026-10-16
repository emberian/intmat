# Learning mechanics as a material design lever

Research/design note, 2026-10-16. This proposes new bounded laws; it reports no new
simulation result. Primary texts were checked directly below. The starting point
is the actual [joint policy](../../studies/native_joint_construction_credit.jl),
not an assumed general learner: three held Gaussian logits choose reserve,
construction investment and catalyst allocation; a paid linear readout uses
local material, inventories and existing hidden/adaptation state. Its recurrent
representation is not trained by the construction reward. The first authority
repair—allowing precursor to survive reserve-demand oxidation—is a more direct
improvement to learnability than collecting another controller metric.

**Recommendation:** first make delayed consequences creditable across arbitrary
window boundaries; then supply a small, explicit temporal representation that
can distinguish rising from falling supply. Keep the native chemistry, affordable
work and local observations. A longer trace cannot repair an unusable actuator,
and longer memory cannot repair an objective that penalizes every useful action.

## What the literature actually licenses

| Verified primary work | Mechanism worth transferring | Boundary of the evidence |
|---|---|---|
| Bellec et al., *A solution to the learning dilemma for recurrent networks of spiking neurons* (2020), [full text](https://www.nature.com/articles/s41467-020-17236-y) | Online local derivative traces plus neuron-specific learning signals; adaptation supplies a slow eligibility component. Reward-based e-prop adds a separate reward-credit filter (Methods, equations 36–37). | The implemented learning signal approximates the unavailable recurrent gradient. An action-score trace is not the same object as a recurrent derivative trace. This is neither a theorem that all local learning works nor a reason to convert Interior to spikes. |
| Miconi, *Biologically plausible learning in recurrent neural networks reproduces neural dynamics observed during cognitive tasks* (2017), [paper](https://elifesciences.org/articles/20899.pdf) | Recurrent weights can learn memory, context selection and coordinated outputs from delayed trial rewards using local perturbation-related eligibility and supralinear increments. | Trials and reward structure are engineered. Its supralinear rule is not the Gaussian score estimator currently used here; substituting it changes the algorithm, not just numerical tuning. |
| Benna & Fusi, *Computational principles of synaptic memory consolidation* (2016), [paper](https://www.gatsby.ucl.ac.uk/~pel/tnlectures/papers/benna_fusi.pdf), [DOI](https://doi.org/10.1038/nn.4401) | Bidirectionally coupled fast/slow synaptic variables protect retained changes while preserving rapid writing. A small chain implements a broad range of retention times. | Their memory-storage benchmark assumes unstructured, uncorrelated modifications. Retention is not reward credit or task learning. A bank of independent exponential sensor filters is a simpler engineering proposal, not an implementation of their consolidation model. |
| Dohare et al., *Loss of plasticity in deep continual learning* (2024), [full text](https://www.nature.com/articles/s41586-024-07711-7) | Selective replacement of low-utility units maintains learning in the studied continual tasks; ongoing feature diversity is an intervention, not merely a diversity score. | Feedforward/deep-RL evidence does not establish that Interior's eight recurrent units have lost plasticity. Their utility is a heuristic. Resetting a recurrent unit can remove dynamics needed by other units; do not import it indiscriminately. |
| Wang et al., *Enhanced POET* (2020), [paper](https://proceedings.mlr.press/v119/wang20l/wang20l.pdf) | Maintain different environments and transfer solutions among them. Performance-based environment characterization and screened transfer reduce dependence on a single hand-designed challenge. | Outer search over bipedal tasks is not autonomous ecological evolution. The paper explicitly identifies limited environment encodings as a ceiling on discovery. |
| Nisioti & Moulin-Frier, *Dynamics of niche construction in adaptable populations evolving in diverse environments* (2023), [paper](https://arxiv.org/pdf/2305.09369) | Environmental diversity, connectivity and construction interact; niches can retain the consequences of previous inhabitants. | Here plasticity is modeled by tolerance-curve width and fitness by its density, not learned recurrent action credit. Do not cite it as a demonstrated algorithm for our delayed chemistry. |
| Natschläger, Bertschinger & Legenstein, *At the Edge of Chaos* (2004), [paper](https://proceedings.neurips.cc/paper_files/paper/2004/file/f8da71e562ff44a2bc7edf3578c593da-Paper.pdf) | Driven threshold-network computation depends on balancing separation and forgetting; a scaling rule changes that regime. | Threshold networks and their chosen tasks do not establish a universal optimal critical point for an embodied learner. Tune regimes by useful control and recovery, not by proximity to a fashionable boundary. |

These identities do not depend on the supplied Matryoshka/PaCE citation strings.
Those papers are not needed for the mechanisms recommended here. No third-party
summaries are used as evidence for the table.

## Separate four clocks and two kinds of value

Interior currently has decision cadence 0.2, return windows 2, source-mask exchange
every 4, and experimental catalyst lifetime 5. The paid pulse's first positive
private sample at 0.4–0.6 says when benefit first appears, not when its marginal
tail ends. Moreover, an action at time 1.8 has only 0.2 before its score is erased.
Thus even a return window longer than the onset lag can discard relevant credit.
This is a structural possibility, not a diagnosis of the pilot's weak outcome.

The local utility `U = battery + usable stored precursor potential` values
inventory now. A catalyst can be valuable because of future flow while reducing
U immediately. Keep physical inventory accounting distinct from an agent's
learned forecast of future return. Adding the mass of catalyst directly to U
would impose a designer-selected construction bonus and can reward useless
deposition. Predicting future *realized* inventory changes is a more defensible
route, with prediction error and horizon made explicit.

There is also a control-variate trap. The current return predictor is frozen
before the current window's exploration, appropriate for its current score sum.
If old scores persist into the next window, the next predictor's observations
may already depend on those old actions. Subtracting that predictor from the old
scores can remove precisely the self-made niche benefit we want to learn.
Do not implement continuing credit by merely deleting the trace reset in
`NativeLifetimeCredit.reward!`. Either retain each decision's pre-action baseline
with its delayed return, or introduce a separately derived value-critic/TD law.
TD errors and a window-return residual are different mathematical objects.

## Experiment 1: remove the arbitrary credit boundary

**Question:** can the joint policy learn a better investment/reserve tradeoff
when a decision receives its own subsequent return interval?

First map marginal consequences in the **current joint task**, rather than
borrowing the earlier secretion pulse's lag. From matched physical checkpoints
just before/after a source exchange and midway through each source context,
perturb one held logit by +0.15 or −0.15 for one 0.2 decision; then return both
branches to the same frozen policy and common subsequent innovations. Test each
of the three actuators separately, with an unperturbed twin. Record producer U,
funded deposition, oxidation and uptake at offsets 0, 0.2, 0.4, 0.6, 1, 2, 4, 8,
plus neighbor U and actual additions. The central difference estimates a local
response, not a universal derivative near affordability thresholds. Preserve
both signs to reveal clipping/threshold asymmetry. Report whether the incremental
U difference between offsets 2 and 8 reinforces, reverses or erases the early
effect. Use at least two founder banks and both feed contexts. This is a bounded
branch assay; it changes no continuing history. Its result selects the horizon
below and can expose a reserve/allocation interaction before changing credit.

Extend `choose_joint!` and the shared paid update through `L.policy_launch!` with
a preallocated circular queue at the decision clock. For decision j store its
score `g_j = x_j (u_j−mu_j)/sigma²`, local U, and a baseline predicted **before**
that decision's noise. After H physical time, its target is `U(j+H)−U(j)`;
apply the paid, clipped update from `g_j × (target−baseline_j)`. Learn the
predictor on matured targets. Clear private queue state on identity change.
This is delayed, truncated score credit with online policy drift; do not call
it an exact gradient of an unlimited continuing objective. Retain clipping,
turnover payments and post-payment return conventions in all arms.

Use H=2 first, then H=8 only if the short pulse/branch assay shows material
benefits still changing after 2. The H=2 comparison changes boundary placement
without changing maximum horizon. Start with four arms: existing window law,
decision-aligned law, its zero-actor-credit control, and reversed actor credit.
Match action-noise keys. Two founder banks × four independent exploration seeds
× four arms fit 32 worlds, with shared geometry and physical strides. Train for
24 time and use the existing four-time learned/erased frozen branches under
continued and shifted supply. An H=8 follow-up needs an explicitly extended,
still finite training/drain budget; do not silently stop before delayed updates
can mature. At the end of training either drain a declared no-new-score tail or
report pending scores and discard them identically across arms.

At h=8, H=8 stores 40×69 Float32 scores, roughly 11 KiB per resident plus
baseline/utility metadata. The tiny four-resident fixtures can afford that;
larger populations should test factorized storage of 23 features plus three
Gaussian innovations per decision, reconstructing the outer product only when
the target matures. Queue processing occurs at decisions, not every physics step.
Because maturation increases update frequency, compare actual paid work and
raw/clipped proposal magnitudes. Keep actor eta unchanged initially: for equal
targets without clipping or policy drift, ten separate `eta R g_j` contributions
sum to the old `eta R sum(g_j)`. Dividing actor eta by ten would suppress that
aggregate gradient. The critic's repeated regression has a different cadence;
reduce its gain accordingly and fit current coefficients to stored pre-action
features, while the actor retains its stored pre-action baseline. Sequential
clipping, policy drift and charging the sum of absolute updates rather than the
absolute sum break exact equivalence; retain a matched-gain control when needed.

Before coupled runs, use a tiny delayed-return fixture where an action at 1.8
causes reward at 2.4: verify correct sign and no update from a reward preceding
the action. Then verify restore midway through a queue. Measure whole-run
time, including setup and observations, using the existing cohort runner.

**Decision criterion:** retain the new mechanism if training-arm utility and
same-state frozen learned/erased benefit improve across exploration realizations
without relying on more external resources or hidden credit bonuses. Report
requested/funded investment, gross uptake, remaining raw feed, material state,
payments and clipping. Frozen branches test usefulness of installed coefficients;
they alone do not compare the full credit laws, whose training material differs.

## Experiment 2: give the policy usable temporal coordinates

**Question:** can cheap local memory make material control responsive to recent
history, rather than requiring the random inherited recurrence to expose the
right feature by chance?

Add three exponentially filtered versions of local raw balance, raw total and
catalyst total at physical timescales 0.4, 2 and 8:
`m_k ← exp(−Δt/τ_k)m_k + (1−exp(−Δt/τ_k))s`.
Expose current-minus-filter and fast-minus-slow differences to the same paid
three-action readout. These are nine extra state floats per resident, updated
inside `choose_joint!` once per decision. They introduce no external phase or
reward. Initialize each filter to the first local observation, so arbitrary
zero initialization does not create a fictitious history. Charge added persistent
state under an explicit study law or include the same state machinery in all
arms; describe this cost choice rather than calling the memory free biology.

Keep readout size and physical work equal across four arms: true ordered filters,
copies of current observations in those slots, native hidden/adaptation features,
and zero added features. Pad slots consistently. This is a test of accessible
temporal representation; it does not remove recurrence from native metabolism.
Use the same 32-world replication budget and the credit law selected independently
in experiment 1, avoiding a simultaneous sweep of representations and estimators.

Use finite, accounted source-mask changes with a predeclared irregular schedule
shared across a cohort, then test a held-out schedule with the same integrated
supply. Four-unit alternation alone can let age-linked recurrence imitate a clock.
For a sharper causal probe, temporarily occlude the *policy readout's* material
inputs during a short gap after a local cue, without moving physical material or
changing the native controller's inputs. Compare ordered memory to current-only
features through the gap, and test a reversed cue-to-future-supply relationship.
Specify this as an engineered partial-observation challenge, not spontaneous
ecological memory. A useful effect must survive unoccluded held-out supply too.

**Decision criterion:** improved paid return and anticipatory funded allocation
on unseen schedules, with loss of that advantage when only the added memory is
erased in a matched branch. If simple filters work, they establish a useful
temporal basis to compress into recurrent motifs or train with local recurrent
eligibility later. If they do not, lengthening neural retention is not yet the
priority; inspect actuator saturation and marginal construction value first.

The smaller implementation-ready successor is now
[paired local trend filters](local-temporal-features.md): six state floats,
three bounded contrast features, explicit descendant handling and a CPU-tested
same-current-input/opposite-history fixture. It precedes the broader three-clock
experiment above and can be integrated into the next mobile/dividing study.

## Modern learning dynamics: two useful translations, no criticality transplant

[Cohen et al. (ICLR 2021)](https://arxiv.org/pdf/2103.00065) measure full-batch
gradient descent near `η λmax(Hessian(loss)) ≈ 2`: loss can oscillate locally
while improving over longer intervals. Their SGD discussion explicitly does not
find the same fixed sharpness threshold. [Damian, Nichani & Lee (ICLR 2023)](https://arxiv.org/pdf/2209.15594)
explain a restoring effect from cubic loss terms: unstable displacement pushes
the optimization trajectory toward lower curvature. Their analysis assumes
progressive sharpening, a single unstable eigendirection and regularity bounds.
These are statements about **parameter-update dynamics on a loss**, not the
state Jacobian of a driven reservoir. Neither threshold is a recipe for physics
dt, recurrence gain, or our clipped, paid, nonstationary score learner. We have
no defensible basis for targeting its Hessian at `2/η`.

**Practical translation: constrain the behavioral size of a paid update.**
Coefficient clipping alone permits very different action changes when many
features align. For proposed readout change ΔM and current local features x,
compute `Δμ = ΔM' x`. With fixed Gaussian noise σ, the conditional pre-squash
policy KL is exactly `K = sum(Δμ²)/(2σ²)`. Scale the proposal by
`min(1, sqrt(ε/K))` before its existing affordability transaction; handle K=0
without division, retain coefficient bounds, and charge only the installed
changes. This needs one small dot product per action, no Hessian or replayed
world. It is an engineering response-budget proposal, **not** a consequence
proved by the edge-of-stability papers. The bound holds at inspected features,
not at all future states or across physical affordability thresholds. Compare
one declared ε against a fixed smaller learning rate at matched cumulative
paid work; retain it only if delayed-control benefit improves. Check predictor
and actor separately: the current normalized linear critic step already controls
single-sample curvature through its feature norm, while the actor's score update
is not gradient descent on that critic loss. This belongs after decision-aligned
credit, not as a simultaneous change to its first experiment.

[Zucchet et al., *Online learning of long-range dependencies* (NeurIPS 2023)](https://papers.nips.cc/paper_files/paper/2023/file/2184d8450c8a641f9a10c49279087c97-Paper-Conference.pdf)
offers a more direct architectural lever: independent recurrent modules make
within-layer forward sensitivities tractable and exact; their multilayer method
approximates learning signals through subsequent recurrent layers. Their
sequence-task results do not establish embodied reward learning.

**Practical translation: train a tiny temporal module rather than a dense RNN.**
After the fixed-filter experiment, replace two filters with a real two-component
mode `z_next = r R(ω) z + B s`, where `0<r<1`, R is a planar rotation and s is
local sensation. Parameterize r inside the unit interval and bound input gain.
This provides adjustable forgetting and phase with inexpensive real arithmetic.
Carry `E_next = A E + (∂A/∂θ)z + (∂B/∂θ)s` for its few parameters. Combine
readout sensitivity with Gaussian score and the delayed return estimator; pay
every installed module/readout change. At fixed parameters and given observation
history the recurrence derivatives are exact; online parameter drift and
truncated reward credit remain approximations. The stable open-loop mode does
not guarantee stable organism–material feedback. Verify sensitivities against
finite differences on a short recorded input, then compare frozen versus paid
trainable modes in the same held-out supply assay. This changes what temporal
features the material can learn, with O(number of module parameters) retained
state, rather than asking a random reservoir to happen upon useful memory.

## Beyond these two experiments

Keep a small outer archive of capable *laws*, organized by actual response:
recovery after disturbance, transport remodeling, reversible/retained structure,
temporal discrimination and energetic cost. Transfer a successful policy between
new, separately initialized supply/decay regimes; retain specialists instead of
replacing everything with the highest inventory scorer. This is a modest POET
transfer principle, not a claim of open-endedness. Preserve geometric and field
snapshots so the archive can later gain descriptors for waves, topology and optics.

If repeated reversals reveal declining ability to relearn, reserve a few
sidecar features for paid replacement/consolidation trials. Do not reset native
recurrent units on the strength of a deep-learning paper alone. Conversely,
training recurrent dynamics becomes worthwhile once readout-only learning
demonstrably hits a representation limit. Derive eligibility for Interior's
actual hidden/adaptation equations, with fixed-rank or local approximations;
do not silently rename the current action trace “e-prop.”

The artistic destination is material whose history changes its next act:
an anticipatory catalyst seam, a finite-stock repair structure that switches
strategy, or optical paths altered by remembered exposure. The same principles
apply beyond organisms: accessible controls, informative local consequence,
appropriate temporal state and diversity of discoverable laws. These experiments
build those capabilities while keeping their first claims bounded.
