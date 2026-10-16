# Give material actions something to learn from

2026-10-16 · Astra engineering design and bounded prototype.

The next useful change is **credit-bearing action plasticity**, followed by a
**competence-progress repertoire builder**. These operate at different levels:
the former can improve a cell's policy; the latter improves our ability to find
worlds worth inhabiting. Neither requires replacing the Julia physics, optimizing
for a single silhouette, or declaring one dynamical regime universally best.

## Reading that changes the implementation

The following full texts were read, including mechanisms and limitations rather
than abstracts. No new publisher PDF is retained in Git. Miconi's CC-BY XML was
retrieved through Europe PMC into temporary storage after publisher/PMC browser
challenges. Williams' publisher version was read via a public PDF mirror; the
first university copy had a broken text encoding. The author-hosted Lenia article
is the 2022 version and explicitly links a later Science Advances extension; the
claims below concern the version read.

- **Williams (1992), [REINFORCE](https://doi.org/10.1007/BF00992696)**,
  §§5–6, equations 11–13 and appendix A: an episodic eligibility accumulator is
  the sum of log-policy derivatives. For Gaussian mean parameters the score is
  `(sample−mean)/variance`. Baselines must not depend on the episode's sampled
  actions. This gives a clean engineering mechanism for delayed feedback without
  differentiating chemistry. Expected-gradient alignment is not convergence or
  low-variance learning; clipping adds a further qualification. [Readable full text](https://gwern.net/doc/reinforcement-learning/model-free/1992-williams.pdf).
- **Miconi (2017), [recurrent reward-modulated learning](https://doi.org/10.7554/eLife.20899)**,
  learning-rule and analysis sections: the published rule accumulates
  `e_ij += [r_j(t−1)(x_i−running_mean(x_i))]^3`, then applies
  `ΔJ_ij=η e_ij(R−baseline)` at trial end. The paper uses trial-type-specific
  reward baselines and clipped updates. Supralinearity suppresses the smaller
  opposite-sign relaxation contributions introduced by subtracting the running
  mean. Merely multiplying our present Oja marks by reward would not reproduce
  this mechanism. Its comparison to explicit node perturbation motivates using
  known exploratory noise first, rather than trying to infer noise from our
  deterministic adaptation dynamics. [Full-text XML](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5398889/fullTextXML).
- **Hamon et al. (2022), [Learning Sensorimotor Agency in Cellular Automata](https://developmentalsystems.org/sensorimotor-lenia/)**,
  curriculum and failure-case sections: sample reachable nearby goals, retrieve
  a prior policy, optimize it, then archive the reached outcome. Their inner loop
  differentiates through 50 Lenia updates toward a target image; obstacle
  configurations provide a recovery curriculum. Initialization matters; diversity
  can be suppressed by the target shape, and multi-creature stability is not
  guaranteed by single-creature training. Transfer the archive/curriculum
  structure; do not import nonconservative Lenia growth into our chemistry.
- **Frémaux & Gerstner (2016), [three-factor learning review](https://www.frontiersin.org/journals/neural-circuits/articles/10.3389/fncir.2015.00085/full)**,
  §4.1: eligibility bridges the gap between activity and delayed modulation.
  A reward baseline averaged across distinct tasks can contaminate Hebbian
  covariance learning. This is a useful warning about confounding resource
  abundance with action quality. It is a review, not an additional independent
  experimental demonstration. Its baseline requirements for general Hebbian
  eligibility should not be mistaken for the weaker unbiasedness condition of
  an exact zero-mean score estimator.

The 2015 arXiv predecessor of Miconi was also inspected; it uses explicit
perturbations and differs materially from the 2017 published supralinear rule.
The implementation below is **Gaussian episodic score learning**, not an
implementation of Miconi's biological rule. Criticality is not a prerequisite:
tune recurrent gain and adaptation only insofar as they improve useful actions,
recovery and behavioral repertoire. Loss-Hessian edge-of-stability arguments do
not select a critical operating point for this recurrent material.

## 1. Action credit, first

**Variables.** Preserve inherited weights `W`, hidden state `h`, adaptation `a`,
and current hidden-row marks. Add a learned action-readout offset `M[i,c,k]`,
eligibility `e[i,c,k]`, perturbation `ξ[i,k]`, and a per-cell baseline `b`.
The implemented prototype trains action rows 4, 5 and 7: uptake-A allocation,
uptake-B allocation and catabolic allocation. Existing softmax competition means
these changes also affect the unused growth-allocation share. This is a deliberately
small metabolically meaningful actuator set, not a new resource-conversion law.

For a cell `i`, pre-action observations `x`, selected action `k`, and timestep `t`:

```
ξ_ikt ~ Normal(0, σ²), independent across cells/actions/timesteps
z_ikt = Σ_c (W_ick + M_ick) x_ict + ξ_ikt
action_ikt = tanh(z_ikt)
e_ick ← e_ick + x_ict ξ_ikt / σ²
U_i = battery_i + efficiency_i Σ_food potential_food * stored_food_i
R_i = U_i(end) − U_i(start)
A_i = R_i − b_i(context)             # baseline before this trial
M_ick ← clip(M_ick + η A_i e_ick, −.3, .3)
b_i(context) ← .8 b_i(context) + .2 R_i
e_ick ← 0                           # at the episode boundary
```

This is an explicit local pre-activity × exploratory post-drive × outcome rule.
The exact Gaussian score uses the unsquashed perturbation; no extra `tanh'`
factor is needed. `Σ e` has no timestep factor: each independent action sample
contributes one score. Changing timestep changes the exploration process, so a
future temporal-resolution comparison must fix the decision cadence independently
of the physics timestep. Holding noise constant while continuing to accumulate
independent-sample scores would be wrong.

**Credit semantics.** Actions run through the actual native controller,
regulation, fair uptake, metabolism, construction, mechanics and accounting.
A score needs no derivative through diffusion, competition or saturation.
The local payoff is usable intracellular chemical potential, measured from
quantities a cell can sense. It rewards acquiring useful food and avoiding
unproductive expenditure, and charges real actions through the existing ledger.
It is deliberately myopic about reproduction, helping neighbors, and durable
construction: those are failure modes and future objectives, not hidden benefits
of this scalar. A caller can later replace the payoff with a learned local value
function and temporal-difference error, but that is a different algorithm.

All physical sources are zero in the prototype. Initial food can still be acquired
without useful policy changes. Therefore absolute positive return proves nothing:
compare frozen trained versus original readouts in matched held-out exposures.
Same exploration noise with zero credit separates noise from learning; reversed
credit checks whether the update direction matters. A common resource influx
cannot by itself explain opposite transfer directions under these controls.
This does not establish causal value for any individual synapse.

**Learning-work semantics.** This first assay is externally trained: each physical
episode starts afresh and only action offsets persist across episodes. It does
not charge computer-side training to a cell's battery, and does not pretend the
reset is conservative physical evolution. Within each episode no reward or
weight operation adds material/energy. Frozen transfer pays the ordinary native
neural/action costs. For an autonomous within-life law, add paid mark turnover:
propose `dM`, compute `L=κΣ|dM|`, pay `min(battery,L)` using `dc_energy_take`,
scale `dM` by the paid/requested fraction, and post the exact paid amount to work.
Do not create a separate unaccounted "learning energy" account. `κ` is a chosen
new material law to calibrate, not a neuroscience-derived constant. Clear credit
traces on identity changes; choose explicitly whether offspring inherit offsets,
reset them, or pay to consolidate them into slowly heritable parameters.

**Executable boundary.** [Study](../../studies/native_action_credit.jl) adds three
KernelAbstractions kernels and invokes `device_cellular_step!` unchanged. The
selected action marks are installed before the step, including noise in the
constant-bias column; the exact observations left by `dc_sense!` supply eligibility
afterward. No copied simulator schedule or alternative chemistry is maintained.
The fixture is an 8³ closed world, four actual founders, F64 inventories/battery,
no growth/division, and the native catalytic/sorptive environment. End-of-episode
identity, finiteness and physical residual gates are mandatory. Suppressing growth
isolates action credit from reproduction; it also removes an important real
tradeoff that must return in the subsequent lifetime experiment.

The harness alternates two initial food compositions, then tests frozen policies
on two unseen compositions with the same total initial food. Context-specific
baselines are initialized by an earlier deterministic rollout; the context index
is supplied by the assay, not inferred by a biological classifier. Training
histories are independent branches; this is not a sequentially depleting dish.
A deployed lifetime baseline must be computed from past local observations,
without an oracle context label or a freshly restored resource world.

**GPU-native plan and cost.** Three selected readouts and 36 inputs require two
extra Float32 tensors (offset and eligibility): `2*3*36*4 = 864 bytes/cell`,
plus noise/baseline; about 14.2 MB at capacity 16384. Installation and eligibility
are independent elementwise kernels with capacity as the first, contiguous
index. Generate noise on device using independent keyed streams; keep traces,
returns and baselines resident. Fuse installation into controller evaluation and
eligibility into regulation once correctness warrants production integration;
use a per-cell loop for paid weight updates. There is no global backward pass or
trajectory tape. The CPU study's restore/allocation and host RNG overhead is not
a GPU throughput claim. Measure complete coupled steps, including paid learning,
and batches per wall second before promoting this to production.

**Minimum next coupled test.** First establish transfer with correct/reversed/no
credit on the present small assay, across additional founders. Then run an
uninterrupted new finite-supply life with reward every fixed physical interval,
paid learning, and a baseline driven by past local observations. At equal age,
branch that learned state into frozen learned versus erased action marks in
matched changed resource patches; retain physical fields and all inventories.
A multiscale trace bank (for example 1, 4 and 16 time units) can span uptake and
slower catalytic returns, but needs separate credit tests: a long trace increases
variance and can blame old actions for other cells' work. For cooperation,
compare cell-private and neighborhood payoff rather than assuming selfish credit
will build shared infrastructure.

**Artistic consequence.** Learned uptake/catabolism can make patches acquire
history-dependent appetite. After a training pulse, cells can respond differently
to the same chemistry; their paid secretions can write this preference into
refractive scaffolds. The first art should expose the same body before/after a
changed nourishment history, with emission tied to measured flux. Training an
ornamental pattern directly is a separate, legitimate outer-design task.

## 2. A repertoire that learns what it can make

Keep a versioned archive `H={(θ,seed,d,q,trajectory)}`. `θ` contains lawful clock
ratios, controller priors and construction parameters; `d` describes outcomes;
`q` records account validity and task-specific performance. Begin with interpretable
axes—recovery time, localized transport change, phase lag and persistent
spatial differentiation—rather than an image embedding that quietly defines life.
Retain raw compact trajectories and add descriptor views when surprises appear.

For descriptor region `j`, measure recent and preceding goal errors under the
same descriptor version and evaluation budget:

```
competence(g, θ) = −distance(d(rollout(θ)), g)
progress_j = max(0, mean(old_errors_j) − mean(recent_errors_j))
P(region=j) = ε/J + (1−ε) (progress_j+δ)/Σ_l(progress_l+δ)
θ_parent = archive member nearest the sampled reachable goal
θ_candidate = bounded lawful mutation(θ_parent)
```

Record both intended and achieved goals. Keep candidates useful for other goals
even when the intended goal fails. Reserve random proposals and revisit dormant
regions; use replicated common-seed evaluation to avoid selecting noise. Advance
from short, cheap goals to recovery and compositional goals only after earlier
ones become reliably achievable. A recipe's "offspring" inherits code parameters
and a reproducible initialization; it does not silently resume or rewrite any
historical ecological run. Distinguish this outer ratchet from hereditary or
within-life change inside the material.

Derivative-free local proposals are the first implementation choice because the
full simulator includes lifecycle events, constraints and stateful accounting.
Smooth controller-only windows may use gradients, but optimizing only a
controller or a terminal image cannot certify a coupled material. The search
can later discover communication by targeting *conditional* response—an intact
bond signal must improve recovery versus a shuffled signal at matched work—then
keep several successful communication styles. Don't prescribe message meanings
before finding whether communication is useful.

GPU cost is primarily the actual rollout; vectorized independent worlds and
resident summaries matter more than fancy host archive logic. Cost `K` bounded
rollouts per candidate, report both single-world step time and candidates/hour.
Archive law validity is a hard constraint, not a reward to trade away. No special
physical energy cost is assigned to the external search; every simulated action
still pays the same local accounting law.

Minimum coupled comparison: equal wall-time and equal rollout-budget searches,
one progress-guided and one random, each from the same initial bank. Compare
held-out attainable descriptor coverage and reproducibility, not just best
survival. Challenge selected members with new perturbation positions and paired
multi-body compositions. Main failures are descriptor blindness, noisy apparent
progress, cheap near-empty "solutions," initialization lock-in, and fragile
single-body tricks. Versioned descriptors, diverse archive parents, perturbation
curricula and explicit composition trials address these with different tools.

Artistically, this produces a cabinet of lawful living-material recipes: pulses,
fronts, scars, appetite, repair and synchronized pores. The archive can expose
paths between neighboring behaviors, including unstable transitions. That is a
creative instrument, rather than a ranking of prettier endpoints.

## Bounded implementation result, 2026-10-16

Three independent four-founder banks (211,217,229), each with48 training
windows per credit/no-credit/reversed-credit arm, completed using the repaired
native metabolism and strict physical accounting. Frozen, noise-free transfer
was evaluated on both held-out food mixtures from restored physical histories.
Mean per-cell improvement in usable energy over each bank's pretraining policy:

| Founder bank | Correct credit | No credit | Reversed credit |
|---|---:|---:|---:|
|211|+0.00108251|0|−0.00231423|
|217|+0.00096332|0|−0.00196165|
|229|+0.00177640|0|−0.00173095|

Correct credit improved21/24 cell×held-out-context outcomes; reversed credit
improved2/24. This is causal evidence that the action-score bridge can improve
this short-horizon intake/work task across these founders. It is not autonomous
lifetime learning: episodes restore physical state, context-specific baselines
are supplied by the study, and learning updates in this first experiment are
external optimization. It does not establish delayed cooperation, ecological
persistence or general learning capacity. The next study explicitly charges
mark changes and retains physical state throughout training.

The [GPU gate](../figures/native-action-credit-gpu.json) verifies the shared
physical repair and all three sidecar kernels, with exact CPU/AMD sidecar parity.
Wide-inventory coupled closure passes; Float32 closure fails in two of four
contexts and remains a recorded limitation. Prior runs keep their original
source pins. The first founder bank differs from the later two only in the
study's JSON3 import spelling, which was adjusted for the pinned GPU environment.
