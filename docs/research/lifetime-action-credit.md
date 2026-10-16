# Paid action credit inside one material lifetime

The [study](../../studies/native_lifetime_credit.jl) wraps the production habitat's
existing sensing and finish hooks. It uses the ordinary Julia cellular schedule,
chemistry, transport, mechanics, and accounting. No production struct is extended.
The selected uptake-A, uptake-B and catabolism rows get a private learned matrix
`M`, sampled-drive buffer, eligibility trace, previous usable energy and local
return baseline. Other controller rows retain their native dynamics.

Every0.2 physical time units, each resident samples `u=(W+M)x+ξ`, with independent
Gaussian `ξ` of standard deviation0.15, and accumulates `E += ξx/σ²`. The sampled
pre-tanh drive is held until the next decision; changing observations between
decisions do not silently change the sampled policy or add score terms. The
implementation cancels selected native weights with temporary action marks and
places the held drive on the constant input. These temporary execution marks
are distinct from the persistent learned matrix. Their cancellation writes are
not charged as material turnover: they only encode a held actuator value through
the existing hook. A production implementation should expose a direct held-drive
hook so this representation does not resemble uncharged physical plasticity.

Every2 time units, locally observed usable energy is
`U=battery+efficiency×(potential_A×store_A+potential_B×store_B)`.
The resident computes `r=U−U_previous`, `δ=r−baseline`, proposes
`M'=clip(M+clip(ηδE,−.01,.01),−.3,.3)`, then updates its baseline with
`baseline'=0.8baseline+0.2r`. It has no food-context labels, outside evaluation
signal, privileged future state, or reset during training. The initial baseline
is zero; the current version uses the first window only to establish that local
baseline, without learning an action update. Eligibility sums ten actual decisions and resets at each reward window;
this is truncated delayed credit, not a claim to solve arbitrarily long delays.
The baseline is measured before the current update. Gaussian score credit is
exact for the sampled conditional action distribution; the finite-window local
objective and continuing physical history limit the optimization interpretation.

Learning turnover costs `κ sum(abs(M'−M))`, κ=.001, measured using the actual
representable Float32 changes. An unaffordable update is rejected wholesale.
Affordable updates debit the battery conservatively and post the exact debit to
native work accounting, including the tiny endpoint-rounding overhead. The next
local return includes this learning work; the final reports include post-charge
inventory. This cost is an explicit experimental material law, not an empirical
estimate of neural biochemistry. It does not charge every computation of a
controller; it charges persistent learned changes, analogous to other model
work costs.

Slot identity is checked at sensing. A new identity erases learned coefficients,
trace and baseline and initializes its own local inventory reference. The
bounded fixture disables growth and tests fixed residents; division inheritance
and a productive birth/death ecology remain subsequent design choices.

Correct-credit, equal-noise/no-credit and reversed-credit arms begin from the
same checkpoint and run uninterrupted for24 time units. After each arm trains,
its physical state is cloned into four noise-free, frozen-policy probes:
learned/erased memory crossed with sham/conserved food displacement. The
half-period displacement moves free and raw food while preserving each species'
quantity and potential; catalyst, bound material and cell state are retained.
This separates the usefulness of learned marks from physical changes accumulated
during training. It is not a cooperation test: the local return can favor
near-term intake and disfavor secretion with remote or late benefits.

The new passes are device-shaped KernelAbstractions kernels with narrow argument
views; state is preallocated. A sensing pass holds selected actions each physics
step, while score work occurs only at decision times and paid updates only at
reward times. Wide energy arithmetic is intentional for the audited AMD/CPU
study. Metal support and complete resident-learning throughput are not yet
claimed. The focused tests exercise holding, payment, affordability and identity;
coupled observations also reject material or potential closure failures.

Artistically, the intended next instrument is matter that changes its appetite
from remembered consequences: a conserved patch displacement can reveal a
learned preference as a transient path or scar. Demonstrating that effect requires
the paired frozen/erased probes, not simply showing a moving pattern.

The preserved [uncalibrated pilot](../figures/native-lifetime-credit-uncalibrated-pilot.json)
retains its [executed source](../../studies/outputs/native-lifetime-credit-uncalibrated/executed-driver.jl)
and completed24 physical units with closure below1.25e−12, but learned marks hurt
same-state frozen transfer in founder bank211: meanΔU−.000160(sham) and
−.000266(displaced food); reverse credit had the opposite sign. The maximum mark
was .00771 after the first window and .00802 at the end, below the .01 update
cap. The local background return changed far faster than the initial baseline:
one no-credit resident returned .0879,.00963,−.00191 over the first three windows.
The baseline's effective time constant is−2/log(.8)=8.96, while the reward window
is2 and action cadence .2. This large common transient can dominate small
policy effects without any numerical instability.

The next bounded version uses one local calibration window: store the first
observed return as the baseline and clear its trace without changing learned
marks. This is a new law tested against newly matched controls, not a correction
to the pilot's result. It removes arbitrary zero-baseline credit; it may still
leave lag-induced variance in a depleting environment. A predictive local
baseline or temporally finer return/eligibility estimator would address that
mechanism more directly than increasing the learning rate.

The [native paid-learning gate](../figures/native-lifetime-credit-gpu.json)
passes on CPU and RX6750XT at physics steps .02 and .01. Both execute exactly20
decisions and2 reward windows in4 physical units. Learned-matrix discrepancy is
at most1.40e−9, paid-work discrepancy1.20e−11, and coupled potential residual
1.43e−14. This verifies the clock and physical implementation, not a learning
benefit or long-run stability claim.

All three [calibrated founder banks](../figures/native-lifetime-credit-calibrated.json)
completed the same24-unit training and4-unit frozen probes. Mean learned-minus-
erased return remained negative for correct credit: sham−.000299,−.000163,
−.000411 and displaced food−.000348,−.000281,−.000433 for banks211,217,229.
Reversed credit was positive in all six bank/probe means; no-credit was exactly
zero. Calibration alone therefore did not improve this bounded lifetime task.
The result motivates a baseline that predicts the local transient from past
observations, with its prediction frozen before the window's exploration noise.
This next comparison must also charge critic changes and report prediction error,
so a different learning outcome can be connected to the intended mechanism.

These lifetime banks share exploration seed6201. They vary founders, not the
exploration realization; repeated signs across banks are not independent evidence
about learning-noise variance. Frozen probes also encounter the depleted late
state, whereas online credit optimizes earlier short returns. A changed optimal
policy can therefore contribute to the observed loss alongside estimator noise.
The next comparison includes another matched exploration seed and measures
baseline prediction error and update magnitude, rather than attributing the
terminal outcome to one mechanism in advance.
