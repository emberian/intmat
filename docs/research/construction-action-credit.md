# Learning what material to build

The paid secretion experiment establishes a concrete opportunity: a resident can
spend precursor and battery on catalyst, remodel local substrate, and recover
more usable potential later. Its first positive sampled private payoff is
0.4–0.6 time in the tested geometry. The two-time credit window already spans
that lag. The next intervention therefore changes the task and actuator rather
than simply lengthening an eligibility trace.

`studies/native_construction_credit.jl` is a new, bounded study law. Four
immobile, nondividing residents inhabit an 8³ physical arena. The native cellular
controller continues to operate. A separate local two-logit policy controls the
habitat regulator's investment and catalyst A/B allocation through the shared
`dh_set_requests!` conversion. Native affordability, precursor consumption,
work, deposition, catalysis, sorption, transport and conservation guards remain
in the production schedule. There is no temporary cellular-weight cancellation.

For resident i, the actor feature vector is

    x = (1, raw_balance, raw_total, catalyst_total)
    μ_a = b_a + Σ_c M_ca x_c,       b = (−0.8, 0)
    u_a = μ_a + σ ξ_a,             σ = 0.3, ξ ~ N(0,1)
    investment = sigmoid(u_1),     allocation = sigmoid(u_2)
    E_ca += x_c (u_a − μ_a) / σ².

The three material observations use the native support sampler and its bounded
physical-scale transforms, with scales (1, 1, 0.5). The actual sampled logits are
held for 0.2 physical time, independent of physics dt. Request amounts still
respond to current resident mass at each physical step. A reward window lasts
2 time; its score is accumulated once per actual decision. This is a truncated
window score estimator: it does not claim to credit effects extending into later
windows. The measured pulse lag motivates this initial horizon but does not
prove it is optimal in the driven task.

The inherited paid lifetime update uses the local change in battery plus usable
stored precursor potential, minus a locally predicted baseline frozen before
current-window noise. Its eight predictor features are bias, local raw balance
and total, scaled usable inventory, two completed-return features, and changes
in the two raw observations. No current noise, external phase, clock, oracle
reward or neighbor inventory enters the prediction. The first completed return
calibrates local history without an actor update. Actor rate is 0.003, individual
increments are clipped to ±0.01, and coefficients to ±0.3. Predictor updates and
actor updates together pay 0.001 battery per absolute coefficient change using
the exact existing conservative energy/work transaction. Unaffordable changes
are not installed; identity reuse clears private learning. Prediction machinery
and persistent coefficient turnover are explicitly distinguished: this is a
study cost law, not a molecular implementation of computation.

The environment funds raw A and B at 0.5 world units each per time, from finite
budgets of 14 each. Opposed smooth source masks are exchanged every 4 time.
Their exact declared totals are equal, so this merely redirects future supply;
it neither teleports existing material nor changes the source quota arithmetic.
Every actual addition is posted by the native resource ledger. No sink is used.
The learner observes the resulting local material and receives no phase labels.
Initial free precursor is 0.03 of each species; initial raw substrate is zero.
Catalyst decay is 0.2, a declared new-run choice: its turnover time of 5 is longer
than the measured constructive lag and comparable to forcing. This prevents the
original slow decay time of 125 from making almost the entire assay a single
initial construction event. It is not a fitted optimum.

The predeclared pilot runs 24 time (three full feed cycles), then tests four-time
frozen learned/erased branches from the same physical checkpoint under continued
and shifted future feed. Actor and predictor learning stop during probes, noise
is zero, and both branches retain identical paid history and remaining budgets.
Matched exploration drives correct, zero and reversed actor credit; the zero
credit arm still trains and pays the predictor. Additional no-secretion and
no-raw arms retain the other study machinery. Reports include actual supplied
material/potential, work, uptake, requested versus funded secretion, material
state, clipping and local returns. Founder replication and exploration-noise
replication must be reported separately.

This tests conditional material investment and allocation, not cooperation or
open-ended evolutionary ratcheting. A successful policy might produce changing
patches of catalyst and sorptive material in response to local supply. An
unsuccessful pilot can identify inaccessible actuators, affordability limits,
credit delay or task nonstationarity without bounding the space of learning
laws. The physical task is useful even if this first credit law needs revision.

The first native gate exposed a reporting-only GPU boundary: the inherited
lifetime observer used host scalar indexing on device state. The construction
observer now explicitly snapshots physical state and copies small sidecar arrays
at reporting boundaries. No host scalar access or state transfer was added to
physical steps or controller kernels. The failed gate's source is preserved;
its partial execution is not counted as a parity pass.

## First engineering result: command authority depends on physiology

The [original pilot](../figures/native-construction-zero-authority.json) completed
all five arms with the four original identities alive. Nevertheless, funded
secretion stopped after its early transient while requested secretion continued.
Correct and reversed credit changed commands but produced identical physical
utility; every frozen learned/erased difference was zero. The primitive budget
priority explains this: reserve-demand catabolism happens before maintenance and
habitat construction. In this resource regime, the inherited reserve targets
remain above available battery, so incoming precursor is oxidized before it can
fund new material. This does not make that reserve strategy universally wrong.

A [paired authority check](../figures/native-construction-authority.json) changes
only the founder's reserve readout to a constant fraction 0.2 of battery capacity
(native action −0.6). This is a new founder phenotype, not a free lifetime update.
It receives exactly the same external resource budget and initially stored
material. At time 8 its battery is approximately 0.96 per resident, stores remain
nonzero, and essentially all cumulative construction requests have been funded.
The original phenotype has empty stores and no further funded requests. Source
additions are 8 total units in both branches; potential residual stays below
4.84e−13. A strict checkpoint restore with exchanged source masks preserves the
physical audit. Lower reserve spends more on material; it is not itself a fitness
improvement.

The [lower-reserve 24-time pilot](../figures/native-construction-fixed-reserve.json)
retains paid construction through time 24: correct credit requests 2.938 total
catalyst units and funds 1.517. Its frozen learned-minus-erased mean usable return
is +7.36e−5 under continued feed and −8.32e−5 under shifted feed. Reversed credit
has −7.04e−5 and +8.06e−5 respectively; no-credit is exactly zero. These are one
founder bank and one exploration seed, evidence of physical and credit authority,
not robust adaptive benefit. All residents survive. No-secretion retains more
usable inventory at time 24 than this construction policy; the initial policy
can overinvest, and source capture eventually limits marginal catalyst benefit.
The complete 24+4 window uses finite funding of 14 per channel, with actual
endpoint additions recorded rather than assumed exact.

A separate joint policy is the next law: native reserve priority, investment and
allocation become three held sampled actions. Its common feature bank includes
local energy/stores and existing hidden/adaptation state, with an explicit
recurrent-feature ablation. This makes the budget tradeoff learnable and permits
within-cell temporal computation. Neither joint control nor memory benefit is
claimed by the two-action pilot above.

The corrected [native gate](../figures/native-construction-credit-gpu.json)
passes all four raw-feed on/off × dt 0.02/0.01 cases through six physical time,
including paid updates and strict restoration with exchanged source masks.
Maximum checked CPU/AMD discrepancy is 6.27e−8 and physical residual 1.10e−13.
This gate uses the recorded two-action source snapshot. Later execution-only
launch adapters and joint-policy extensions require their own composition checks.
