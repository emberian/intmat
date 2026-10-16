# Shared recurrent learning in native material

The implemented target is a trainable recurrent controller, with authored
programs retained as a [baseline](predictive-programs-baseline.md). The
user-provided ten-page *Big Blob of Compute* essay changed the direction:
capacity, training signals and varied experience should shape a flexible
observation-to-action computation, rather than forcing learned behavior into
our finite program vocabulary. The essay is a design perspective, not evidence
that a particular scale will succeed. Its private PDF remains outside Git.

## What runs

`shared_recurrent_runtime.jl` implements a shared dense tanh recurrence and
linear output heads in KernelAbstractions. `native_recurrent_agent.jl` supplies
native observations, bounded commands and partition-corrected rewards.
`recurrent_training_objective.jl`, `recurrent_training.jl` and
`train_native_recurrent.jl` supply the loss, outer update, rollout collection and
independent held-out evaluation. No gradient crosses the physical simulator.

    h_t = tanh(Wx*x_t + Wh*h_previous + b)
    outputs_t = Wy*h_t + by

The 25 inputs are six species-resolved local material observations; six body
features (relative battery, two relative stores, bounded mass, bounded age,
last-interval funded secretion per mass); and thirteen preceding native/material
commands. There is no source-phase clock, future source information, velocity,
explicit previous reward or uptake-enzyme input in this first schema. Absolute
physical quantities remain represented through the adapter's fixed scales.

Sixteen outputs comprise three action means, one value estimate and twelve
next-observation forecasts. The forecasts see h, **not the current sampled
innovation**: these are passive on-policy predictions, not a counterfactual or
action-conditioned model. Their targets are the next actual first twelve inputs.

Every .2 physical time units, independent fixed-sigma Gaussian innovations are
added to action means. Sigmoid transforms map investment to [0,.3], composition
to [0,1], and reserve fraction to [.02,.25]. Initial head biases give nominal
investment .15, balanced composition .5 and reserve .035; small random head
weights make actual initial commands slightly state dependent. These are
explicit bootstrap bounds, not claims about the final useful action space.
Native uptake, motion, adhesion and division still follow their existing
controller. The generic adapter's broader channel vocabulary is an extension
path, not evidence that this initial learner controls all thirteen channels.

Native physics pays construction and other physical work. Reward is the change
in developmental value, with actual child-value gifts correcting parent
partition and zero retained value on death. The biomass coefficient is the
previously derived replacement accounting convention, not recoverable biomass
energy. No catalyst, prediction-accuracy or survival bonus is added.

## Training and memory semantics

One fresh rollout is followed by one synchronous n-step actor-critic update with
Adam. Width 64 and eight decision steps (1.6 physical time units) are the initial
integration configuration; width 16 supports debugging. Actor, value and
prediction gradients all flow through the same recurrent weights. This is not
readout-only fitting. The actor-critic family is motivated by
[Mnih et al., 2016](https://proceedings.mlr.press/v48/mniha16.html); the synchronous
scheduler and shared predictive loss here are our implementation choices.

The current defaults are discount .995 per decision, value squared-error weight
.5, prediction squared-error weight .1, entropy weight .01, fixed Gaussian sigma
.3, and global gradient-norm clipping at 1. Forecast errors use declared positive
scales (currently ones for the normalized adapter targets). Returns and actor
advantages are frozen during differentiation. The transformed-action log
probability includes its sigmoid Jacobian. Entropy uses a reparameterized
Monte Carlo estimate of the bounded action distribution, keeping the sampled
innovation fixed; it is not mislabeled latent Gaussian entropy. Evaluation uses
mean actions without exploration.

Weights and Adam moments persist across external environment epochs. Hidden
state belongs to each stable resident identity. Incoming recurrence masks reset
newborn/reused slots and stop backward credit; outgoing return masks separately
terminate death/replacement bootstrap. Prediction targets are masked across
identity changes. Inactive observations are finite zeros and receive no loss.
Time-limit truncations bootstrap surviving residents rather than creating a
fictional death.

After an update, the driver recomputes the just-collected prefix under the new
weights to refresh carried hidden state. Its starting hidden state still
summarizes older experience under older weights: this finite-prefix
approximation is explicit. Memory may persist beyond the unroll, but gradient
credit does not. There is no off-policy replay or claim to solve arbitrary
long-delay credit.

This is **external shared training**, not paid local biological plasticity.
Central optimizer compute and shared weights are training resources. Per-cell
lifetime adaptation and paid inheritance are separate future experiments; they
must not be inferred from centralized weight updates or transient recurrence.

## Environment experience and assessment

The driver declares training seeds 211/213/215 with source intervals 6/8/10 and
held-out seeds 212/214 with intervals 7/9. These are sequential independent
environment epochs, not simultaneous world-balanced batches. Each episode has
explicit finite inventories and source budgets; reset is an external operation
that clears h and preserves learned weights. No existing continuing run is
reset. The native adapter performs one declared .02 initialization step before
the first learning observation.

At evaluation, learned and initial frozen weights receive the same held-out
environment and deterministic action rule. Reports include developmental return,
resident decisions, population, funding and physical ledger closure. A one-block
composition gate checks that the system runs; it is not evidence of learned
benefit or generalization. Default longer epochs must cross actual source
switches before their different intervals constitute useful temporal experience.

The next experiment is one measured training allocation with these varied
worlds, followed by fixed held-out evaluation. Compare curves against the
initial recurrent policy and authored-program baseline at equal environment
experience, disclosing compute differences. Then raise width to 128 and unroll
to 32 only against complete-pipeline throughput and memory measurements. Test a
matched memory-reset intervention and a separately trained feedforward control
before attributing benefits to history. Expanding action channels and experience
diversity should use this same training pipeline rather than spawning another
collection of disconnected narrow learners.

## Runtime and proof boundaries

The verified AMD runtime lacks rocBLAS/MIOpen; the implementation therefore uses
native KernelAbstractions forward, exact truncated backward and optimizer
kernels without adding a framework dependency. A separately pinned library
reference or GRU is a later alternative, not an implemented feature. CPU tests
check objective output derivatives, composed recurrent-weight derivatives,
birth/death return cuts and that actor and prediction losses both reach the
recurrent trunk. The runtime lane independently proves its backward pass and
native parity. There is no full-simulator differentiation.

Measure a complete warm rollout plus backward/update, including native physics,
sensing, masks, transfers and births. Component kernel timings alone are not
learning throughput. Report fresh resident decisions and simulated world time
per wall second, peak memory, and the amount and distribution of actual
experience. Capacity gates and short failures are engineering evidence, not a
ceiling on what this architecture can learn.

## Saved integration status

The outer objective/trainer passes 72 CPU checks, including finite-difference
output and composed recurrent-weight gradients, native-objective parity, masks,
and option validation. The [first composition result](../figures/recurrent-first-composition.json)
is width16/T8: 192 native training decisions, one Adam update, then 192 held-out
decisions for each updated/initial policy. The first cold rollout/update/prefix
refresh took 3.323 seconds. Held-out developmental totals were −3.10502 and
−3.11031 respectively, with actual funded construction and ledger residuals
below 2.2e−14. This is wiring evidence, not a learning-benefit claim.

Full native GPU objective/adapter/training parity and warm pipeline profiling
remain unfinished. A staged performance attempt was stopped before evidence
collection. No substantive training allocation or memory-reset evaluation was
launched. Preserve the ignored `outputs/recurrent-first-composition-v2` source,
parameters and report; do not overwrite it.

CPU composition restart, using a fresh output directory:

```julia
include("studies/train_native_recurrent.jl")
TrainNativeRecurrent.run(
    "experiments/native/living-temporal-construction-16.toml",
    "outputs/recurrent-composition-next";
    H=16, T=8, train_blocks=1, eval_blocks=1,
    specs=(training=[(seed=211,phase_time=6.)],
           heldout=[(seed=212,phase_time=7.)]))
```

Before a substantive run, freeze the complete executed project/source epoch,
finish the coordinated native gate, and budget from complete-pipeline timings.
Prefer fresh epochs of at most 25.6 physical time units so training crosses
source switches without unintentionally spending most experience after the
finite feed budget ends around32. Add a learned-policy memory-reset evaluation
before interpreting recurrent state as useful history. These are remaining
steps, not implemented assessment results.
