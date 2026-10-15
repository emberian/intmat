# Directional transport instead of a smaller biological timestep

This unloaded CPU study replaces only transport in a complete device-cellular
constructive ecology. It does not change the production process or live garden.
The biological/controller/lifecycle clock, paid construction, catalytic release,
closed-form sorption, and 16 mechanical substeps are shared across candidates.
The explicit control subcycles **only** diffusion to its original CFL bound.
Thus a transport speedup cannot come from omitting controller or cell work.

The experimental step retains the native nontransport limit
`min(0.25, response_time/2)` plus a stricter study neural/adaptation guard.
The fixture is explicitly closed and refuses nonzero source/sink rates.

## Operators and physical contracts

For an edge joining masses x,y, the native symmetric arithmetic face coefficient
is retained. Its exact two-site flow preserves x+y and multiplies x−y by
`exp(-2 D face dt / h²)`. Disjoint periodic matchings compose without atomics;
a forward/backward palindrome gives second-order splitting for a *held*
permeability field. This prototype explicitly refuses odd periodic dimensions.
Positivity does not ensure useful large-step accuracy: a fixed number of local
matchings cannot represent arbitrarily wide diffusive spreading.

The stronger candidate solves each periodic coordinate line implicitly:
`(I - dt L_axis) u_new = u_old`. Face coefficients are the same heterogeneous
symmetric coefficients as the native stencil. Each line matrix has positive
diagonal, nonpositive off-diagonal and row/column sums one. Its inverse is
nonnegative and preserves total mass in real arithmetic; successive coordinate
solves retain these properties. The implementation uses `SymTridiagonal` LDLᵀ
and a cyclic rank-one correction, with a separate two-voxel periodic case that
retains both parallel faces. Odd line lengths are supported. Solves use Float 64
workspace and return Float 32 fields, so actual rounding residuals are measured.
No unrecorded clamping or discarded mass is used.

These line sweeps are **first order in time**, both from backward Euler and from
noncommuting coordinate operators. Heterogeneous material makes cross-direction
commutators nonzero; it is not a separable exact heat solver. Refining the number
of full directional cycles reduces this error. Permeability is frozen over each
transport interval in both the subcycled explicit control and candidates. This
comparison does not validate holding a rapidly moving/material-changing medium
for an arbitrary interval: the unchanged outer biology/material clock remains a
separate accuracy constraint. Catalytic reaction is still the existing bounded
native update and is not granted a new accuracy guarantee.

## Scope of execution

`studies/pair_transport_ecology.jl` supplies both operators and a complete coupled
carrier. `test/pair_transport_ecology.jl` checks periodic edge flow, zero barriers,
mass/positivity, Fourier convergence and independent dense cyclic-matrix oracles
for lengths 2,3,7,8 and time intervals 0.001,1,100. Full coupled smoke steps include
uptake, metabolism, construction, chemistry, movement and lifecycle.

The line implementation is deliberately a CPU reference and refuses device
arrays. A GPU implementation would require batched cyclic tridiagonal solves
(e.g. parallel cyclic reduction) and separate complete-loop parity/timing gates.
No GPU throughput is inferred from the CPU result. The existing implicit phase
solver uses FFT diagonalization with constant mobility, and implicit contact
studies use a different mechanical CG/Newton problem; neither supplies a reusable
heterogeneous periodic line solver here.

## Initial coupled accuracy frontier

The complete closed fixture has extent 4, 8 founders, capacity 64, inherited neural
controllers, paid enzyme construction, reversible binding and heterogeneous
body/material permeability. A declared ±50% sinusoidal initial food/raw profile
makes transport active immediately. Rates, radius, geometry and initial masses
are shared. The test regime uses diffusion 20, outer interval 0.01, 4 outer steps
(physical time 0.04); this is intentionally a faster-solute regime than the live
garden's diffusion 0.001. At 24³, explicit transport requires 48 substeps per outer
interval, per raw/free pool. All candidates retain four biological updates and
16 mechanical substeps per biological update (64 across this whole four-update
trajectory).

At 24³, one line cycle gave free-food L2 difference 3.73%, raw 2.24%, bound 1.23%
against the transport-subcycled explicit control. Four cycles reduced these to
0.949%,0.634%,0.230%. Pair matching at one cycle was fast but had 22.9% free-food
and 21.7% raw error; even four cycles retained 7.16%/7.45%, so it is not the useful
candidate in this regime. All variants remained finite and retained 8 cells.
These are differences from the chosen explicit discretization, not a rigorous
bound against the exact nonlinear ecology. The full explicit control itself uses
the native CFL 0.45; independent Fourier/dense-matrix tests address operator
correctness but do not replace a future full-trajectory reference refinement.

Local one-thread timings varied substantially with other desktop work, so those
raw samples are retained without a reliable performance claim. The pinned hbox
CPU gate is separate. No mature ecological equivalence, source-fed behavior,
long-run stability, artifact quality, GPU performance or checkpoint support is
claimed by this short closed trajectory.

## Pinned one-thread CPU result

A separate immutable df33fa7 source archive plus this study ran on hbox CPU 8/9
with one Julia/GC/BLAS thread,8GiB/no-swap/nice 15, and no GPU imports. All 52
focused tests passed. Three warmed 24³ complete trajectories took 0.37974,
0.37986,0.38016 seconds for the explicit control. Four line cycles took 0.05732,
0.05746,0.05745 seconds: **6.61× by ratio of medians**, with the approximately 1%
free-food and 0.63% raw-food differences above. One line cycle took 0.02307,
0.02362,0.02308 seconds (16.46×) with the larger 3.73%/2.24% differences. The
conservative recommendation is the four-cycle accuracy point, not the largest
speedup. All raw samples and exact source/runtime hashes are recorded in
`docs/directional-transport-results.json` and its referenced local artifacts.

The four-cycle material closure residual was 9.90e-6 and potential/work/heat
residual 5.93e-5; explicit-control values were 1.97e-5 and 1.18e-4. These residuals
are measured Float 32 endpoint rounding, not exact bitwise mass conservation.
The first combined CPU launcher completed 24³ but failed when reincluding a
module for 64³ due to an unqualified exported name; the standalone 64³ invocation
uses the same immutable numerical source. The current driver qualifies the
study module explicitly to avoid that repeated-include ambiguity.

At 64³ with the same physical extent/rates/time and eight founders, explicit
transport requires 342 substeps per biological interval. One warmed complete
trajectory took 50.44 seconds, versus 1.158 seconds for four line cycles:
**43.55× in this single CPU sample**. Useful physical time per wall second was
0.000793 versus 0.03453. Free-food L2 difference was 1.75%, raw-food 0.592%,
bound-food 1.16%. The material residual was 9.85e-6 and potential closure
5.91e-5 (explicit control 3.71e-5 and 2.23e-4). One line cycle was faster again
but had 6.27% free-food error; it is not the recommended accuracy point.

This is a substantial representation-of-time gain in the declared fast-solute,
transport-dominated regime, not a claimed speedup for the diffusion 0.001 garden.
The short reference follows the same biological clock, so the gain does not
buy more controller, birth or contact resolution. The 64³ run has only one warm
sample; it establishes feasibility/scaling, not a timing distribution. Combined
remote CPU process wall time was 131.10 seconds within the approved 180-second
budget, and both CPU slots were released. A GPU batched cyclic-line solver is
the concrete next implementation step if this accuracy tradeoff fits an intended
material; neither a new GPU benchmark nor a production law change is included.
