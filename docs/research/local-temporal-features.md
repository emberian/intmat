# A local trend sensor for the next living material

2026-10-16. Implementation-ready study component; no ecological benefit claimed.
[Helper](../../studies/local_temporal_features.jl) and
[focused tests](../../test/local_temporal_features.jl) pass 30 CPU assertions.
No native GPU gate or integrated mobile/dividing run has been performed.

The native controller already retains hidden state and adaptation; its regulator
updates those with response/adaptation clocks and activity-dependent hidden-row
marks. Construction reward currently trains a separate readout, not those
recurrent dynamics. Extra recurrent readout features can increase policy gain
without supplying useful history. The proposed sensor exposes one simple,
interpretable distinction: is local material rising or falling?

For each resident, track fast and slow exponential filters of the three existing
bounded observations: raw balance, raw total and catalyst total. At each actual
decision, update with `alpha = -expm1(-decision_interval/tau)` for tau 0.4 and 2.
Expose three features `(fast−slow)/2`. They are bounded by one when inputs are
bounded by one. Six Float32 values plus one Int32 identity cost 28 bytes per
resident. No raw trajectory, global clock, external source phase or future
information is retained. These deliberately separated timescales are starting
design choices, not fitted optima or a claim that they are biologically special.

## Integration contract

Include the module once; allocate `state(native_inputs, capacity)` in the study
sidecar, alongside its existing learned/queue state. Compute `coefficients(.2)`
at setup. Inside the decision kernel, after material sensing and before sampling
an action, call `observe!` with `C.feature(s,i,h,2:4)` supplied individually.
Append `contrast(memory,i,k)` for k=1:3 to the action readout. Store the actual
feature vector with each decision's score, as for all cohort features. Never
update the filters at the finer physics cadence while using decision coefficients.
The inline helper avoids another launch; the standalone kernel is for tests.

On a new ID, initialize both filters to its first observation, yielding zero
trend. For division, call `inherit_kernel!` after native division and before the
child's first observation. Explicit `copy_parent=true` copies normalized
information when the parent ID matches; false resets the child. The parent is
unchanged. This follows the native precedent of copying adaptation, but does
not claim information is a conserved inventory. It is separate from whether
learned coefficients are copied, and old action-credit queues must not transfer
to the child. In the first mobile run choose reset, keeping descendant policy
inheritance and sensor-history inheritance as separate design decisions.

The helper changes no energy accounting. Both history and control arms must run
the same filter state machinery and pay the same declared state-maintenance law;
coefficient updates retain the existing paid transaction. If the first study
only charges coefficient turnover, state that limitation explicitly. Do not
claim the helper itself implements paid biochemical memory. Register all three
arrays in sidecar checkpoint/projection handling before a native cohort run;
the current tests verify array-copy continuation, not the full checkpoint codec.

## One focused causal test

The executable fixture presents histories with raw balance +1 versus −1, while
raw total and catalyst total remain 0.5. Both then observe exactly `(0,.5,.5)`.
The temporal vectors have equal squared norm and opposite first components;
the same fixed readout produces equal-magnitude opposite actions. Current-only
copies cannot distinguish the pair. This establishes controllable capability,
not learned usefulness. It also verifies neutral initialization, ID replacement,
constant-input half-step agreement, inactive slots and explicit copy/reset
inheritance without parent mutation.

For the integrated test, train one paid policy, then freeze it and branch one
complete physical state. Keep native hidden/adaptation, energy, inventory,
position, learned coefficients and future innovations identical. Supply two
filter histories obtained from matched opposite-balance cue sequences with the
same final observation; compare against a neutral-history branch. Match added
feature norms, or report their discrepancy before interpreting gain. A difference
in held logits demonstrates use of history; a subsequent difference in funded
construction and realized return tests utility in that physical context. This
sidecar-state intervention is not an attainable natural history claim. Follow it
with an unmanipulated held-out supply sequence before declaring adaptive benefit.

Prefer this small mechanism over a broad clock sweep. If it supplies useful
control, its timescales can become tunable recurrent modes later; if not, retain
the result without treating the organism's entire memory as ineffective.
