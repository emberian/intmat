# Explicit clocks for constructive ecology

These studies change temporal factorization, not spatial resolution or physical
radius. They are experimental operators outside the production process and
checkpoint protocol. No live run or default integrator is changed.

## A control: slow material, fast complete cells

`studies/multirate_material_ecology.jl` advances the existing controller, uptake,
junction exchange, metabolism, paid export, mechanics, contacts and lifecycle
every delta. Fine-lattice diffusion and catalytic/reversible-sorption chemistry
advance once per H=m*delta. Body porosity is rebuilt on every fast step; cached
material permeability is multiplied back after that rebuild. All pools retain
the native Float32 representation and fine spatial grid.

Uptake and construction debit real stocks immediately on every fast step. No
unfunded request buffer promises future resource availability. Fresh enzyme is
paid immediately but its catalytic processing/decay waits for the outer event;
bound food is unavailable for uptake throughout, as in the native law. The raw
and free pools diffuse at the start of an outer interval, and chemistry finishes
at its end. This is an explicit Lie split with changed release latency.

H must satisfy the original diffusion bound; delta satisfies neural and mechanical
guards. At m=1 all compared physical, scratch, identity and counter arrays must
match the native loop exactly. Every fast step retains its original RNG/event
index. For m>1 biological decisions remain fast, but delayed chemistry can change
their outcomes; no event-history identity is promised.

Two useful local diagnostics have explicit limits. Holding old enzyme without
decay for H can overestimate it by at most exp(decay*H)-1 relative to its decayed
value. Treating fresh paid enzyme as present throughout H can add at most
kcat*H*(actual fresh enzyme mass) to the requested raw conversion, since the
saturation fraction is at most one. This is not a bound on all splitting error
or the resulting nonlinear ecological trajectory. The study records both bounds.

The matched AMD profile attributes only3.07ms diffusion and2.10ms reaction/sorption
to a38.67ms materialized-support step. Even removing those operators entirely
cannot produce a large speedup; after the roughly10.5ms closed-resource work
was removed, this subset's ideal infinite-cadence ceiling is about1.22x. A
CPU-specific improvement must not be generalized into a GPU arena acceleration.

## A broader operator: fast neural/contact dynamics, slow biology

`studies/multirate_resolved_mechanics.jl` reuses the shared native enqueue schedule
through a narrow experimental carrier. Mechanics retain the original delta/16
substep, with funded motor debits at every substep and bond capture every delta.
The default experimental variant holds the controller for H too. Its optional
`fast_controller=true` variant keeps recurrent hidden state, adaptation,
plasticity, actions, neighbour communication and physical cues at delta.

Only the eight chemical sensor channels are held between outer observations.
Native sensing/regulation kernels are split from their actual source in
`studies/multirate_cached_controller.jl`; marker checks reject incompatible source layouts. Compatible formula changes
follow the native source and require repeating the gates.
The first neural update precedes the H-sized uptake/metabolic stage; remaining
neural updates accompany the fast mechanical intervals. This performs exactly
m neural updates of delta per H, rather than double-advancing the network.

Sensing/support refresh, fair uptake, junction exchange, metabolism, catalytic
planning/export, diffusion, sorption and lifecycle operate at H. Uptake is a
funded zero-order policy evaluated at the outer event, not an integral of all
intermediate fast actions. Growth/construction share the existing budget and
motor work always debits actual remaining energy. There is no speculative
exchange that can spend a donor twice.

This preserves a resolved contact/motor clock while deliberately introducing a
slow biological response. Birth, death and genotype mutation occur at H events;
transient threshold crossings can differ, so there is no blanket bounded-delay
claim. Persistent contact capture is checked every delta, although adhesion
actions can differ because the network's chemical information is held. Existing
force/displacement guards remain active. At ratio1 both controller variants must
match every native array exactly.

Unlike slowing material alone, this broader split can remove repeated support
construction, chemical sampling and global allocation work. The pre-reduction
profile spends about4.7ms on mechanics and about1ms on the controller, while much
of the remaining work could move to H. That supports investigating a several-fold
whole-workload gain, not claiming one before a controlled native GPU run.

The prototypes use exact integer clock ratios and report physical duration
separately from step count. A failure leaves a partial operator state; there is
no rollback, generic composition or restart support. Only closed reservoirs
are presently nominated for comparison. Source-driven budgets and an escrowed
integral of fast control requests require a separately declared operator, not
silently reusing these zero-order policies.

## Working coupled CPU comparison

The material-only control passes8 assertions, including exact ratio1 arrays,
conservation and clock alignment. The resolved-mechanics variants pass14 coupled
assertions, including exact ratio1 arrays with the fast recurrent path, plus a
forced simultaneous death/division. The daughter inherits from surviving ID2,
starts with zero secretion observation, and preserves the closed material/work
budget. Six additional assertions cover actual bonded neural/mechanical coupling
and funded movement. Two invalid-duration guards are also checked, for30 total
observed assertions across the study gates.

The stronger operator was compared twice on native32³,32 founders/capacity128,
delta=.001,16 mechanical substeps and fixed duration .256. All32 identities stay
unchanged in these timing trajectories. They exercise the complete coupled loop,
but are not mature contact-rich ecological populations. Initialization, audits
and final downloads are outside the timers; target-shape warmups precede them.

| Outer ratio | Native-loop speed ratio | Ratio1-carrier speed ratio | Food A/B relative L2 | Bound field relative L2 | Hidden-state relative L2 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2 | 1.68–2.17x | 1.34–2.04x | .420% / .389% | 2.55% | .0327% |
| 4 | 4.04–4.14x | 3.23–3.88x | 1.32% / 1.34% | 8.72% | .106% |
| 8 | 6.30–6.73x | 5.38–5.91x | 3.49% / 3.55% | 24.9% | .311% |

The exact ratio1 carrier itself varied from1.07 to1.25x the native timer despite
identical arrays, so the second speed column is an additional control for timing
noise/host effects. These CPU measurements establish a useful complete operator,
not an AMD/Metal throughput claim. The measured GPU stage profile suggests a
smaller potential gain, roughly2–3x at a high ratio, requiring an actual native
whole-loop benchmark.

The reference bound inventory is .0264464; ratios2/4/8 yield .0268518/.0277286/
.0299336, differences of1.53%/4.85%/13.2%. At ratio4 maximum periodic position
error is .000233, around0.13% of exchange radius. Material residual is4.32e-6;
potential/work/heat residual is-1.37e-4. No ecological-equivalence conclusion is
drawn from these short cases. Ratio2 is the conservative next performance/accuracy
probe; ratio4 is a useful explicit two-clock material experiment. Ratio8 materially
changes bound chemistry and should be treated as a different response regime.

[Recorded comparisons](multirate-ecology-results.json) include all timing samples,
source pins, material-only and slow-controller controls, pool errors and actual
physical duration. Raw observations and source copies are retained under
`outputs/multirate-cpu-prototype-20261015`.

## What integration would require

The fast neural variant holds a chemical policy rather than integrating fast
flux requests. A future accumulated-request operator should collect positive
per-cell rates over fast steps, then allocate actual donor-limited grants fairly
at the outer event. It must specify whether those grants become a new owned
feeding escrow or immediately available intracellular stock; simply averaging
actions is not equivalent because enzyme allocation and kinetics are nonlinear.
Both alternatives need an explicit event/availability phase and conservation
ledger. The present prototype does not pretend that extension already exists.

A production version should extract shared sensing, neural regulation and request
operators rather than retain experimental source slicing. Its strict checkpoint
must include both clocks, phase, held inputs, controller-derivation source pin and
biological RNG/event index. Saving only the underlying native carrier would lose
these laws. Native GPU failure/ownership/composition and whole-step parity gates
remain necessary. No GPU slot has been requested or used for this study.

## Native AMD gate with Float64 battery

The bounded 64³ closed-world gate used 256 founders, capacity 4096, materialized
physical support and the native fine lattice throughout. Its actual funded
mass/minimum-mass headroom bound was below capacity. Each trajectory covered
256 fast steps of Float32(.001), physical duration 0.256000012; neural updates
remained at that fast clock and mechanics at its sixteenth. Biological/material
updates occurred every 1, 2, or 4 fast steps. Float64 battery storage was selected
explicitly. This is a temporal approximation, not production acceleration.

Three alternating warm whole-loop pairs gave median wall times 2.296 s native,
1.768 s ratio 2 and1.508 s ratio 4: **1.30× and1.52×** by ratio of medians. Simulated
time per wall second was 0.1115, 0.1448 and 0.1697. Timings include all coupled
stages and identical status checks every 8 fast steps; construction and detached
final observations are outside. One initial ratio4 sample was 1.841 s, versus
1.506/1.508 s on the following pairs, so the individual samples remain recorded.

This gain has a measured cost. Food-species field L2 errors were 1.14%/1.09%
(ratio 2) and 3.83%/3.72% (ratio 4) against native at the same physical time.
Bound-food field L2 errors were 2.70%/9.80%; their absolute L1 errors were
0.00463/0.01627 against native bound inventory 0.20898. Catalyst field L2 errors
were 0.143%/0.397%. These are short-time errors, not a mature ecological
trajectory or long-time equivalence result. All 256 identities survived this
performance fixture; separate prefix gates force simultaneous death/division
and actual reciprocal contact communication.

Ratio 1 matched all native state and scratch arrays exactly on AMD in ordinary,
lifecycle and contact gates. Ratio 4 CPU/AMD maximum floating differences were
1.19e-7,2.38e-7 and5.96e-8 respectively. Cross-backend comparison decodes atomic
catalyst UInt32 storage as Float32 and excludes nonphysical linked-list order
and cached-radius bits; integer identities remain exact. Final full-trajectory
material residuals were 6.81e-5/3.48e-5/1.80e-5, and potential/work/heat closure
residuals 4.06e-4/2.09e-4/1.08e-4 over approximately 3426 potential units.

The immutable numerical base was 0159111, with only the study/fixture/driver
overlays recorded in `docs/multirate-amd-results.json`. It predates the newer
production lifecycle integration. The pinned AMD runtime used the process-local
GPUCompiler opt_level 1 workaround, not a different simulation implementation.
The complete invocation took 110.83 s including compilation and gates under a
300 s, 8 GiB/no-swap, CPU 10/11 budget. Wrapper checkpoint/restart remains unsupported.
