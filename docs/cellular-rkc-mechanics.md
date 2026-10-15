# Stabilized explicit contact mechanics (unloaded study)

This is a separate numerical experiment, not the default integrator or a
migration of an existing run. `studies/cellular_rkc_mechanics.jl` implements the
second-order damped Chebyshev RKC recurrence, with damping 2/13. Its negative-real
linear stability interval grows quadratically with stage count. The derivation
follows [Sommeijer, Shampine and Verwer, section 2](https://ir.cwi.nl/pub/4764/04764D.pdf).
A linear stability interval is not a nonlinear contact guarantee: compressed
spheres have transverse negative curvature. Every stage therefore checks travel
against its conservative pair-cache skin, and each completed passive solve must
pass the actual mechanical-potential decrease gate. Motor motion uses symmetric
splitting with the original remaining-outer-interval funding horizon, explicit
work reservation, and endpoint potential accounting.

The actual middle 3D fixture has 2,393 cells and 6,242 bonds. Its contact/bond graph
is one connected component, so component clocks cannot isolate a cheap quiet
majority. Over one physical interval 0.1, the default needs 62 substeps. Against
the retained 1,024-substep explicit reference, default position RMS is 3.519e-5,
maximum 2.919e-4, and work error -0.005131. RKC with 16 macro intervals and 48
force evaluations gives RMS 1.085e-5, maximum 2.044e-4, and work error +0.004509.
Both retain exactly the reference bond set and pass the mechanical work ledger.
The reference's 512-to-1024 RMS difference is 2.268e-6. Fewer RKC intervals are
cheaper but less accurate; the 16-interval candidate was selected before timing.

The first two implementations lost badly: roughly 543 ms and 483 ms against
245 ms default, due to source-wrapper broadphase rebuilding and boxed loop
state. The cached funded wrapper removes those rebuilds and reuses production
contact-workspace bond columns. These changes leave all five tested RKC physical
states bitwise identical, including funds and bonds. Under bounded CPU8,9,
2 Julia threads/1 GC/1 BLAS, 8 GiB/no swap/nice15 on hbox, five alternating warmed
pairs gave median default 231.157 ms (16.255 MB allocated), candidate 122.940 ms
(33.917 MB). Median paired speedup is 1.879. This is a **mechanics-only** fixture
result, not a complete-step or ecological speedup.

Evidence is retained in `outputs/cellular-rkc-study-20261015/host-cached-62fce99`
and immutable host directory `/othersys/tmp/intmat-rkc-cached-62fce99-20261015`.
Earlier failures remain in the corresponding original/typed directories. Core
source was pinned to 62fce99; each probe records exact study sources and input
hashes. No existing trajectory or pin was changed.

`studies/cellular_rkc_ecology.jl` is an opt-in source-derived physiology adapter:
it checks loaded source hashes and replaces only the mechanics call. It supports
unit drag, zero extra bond range, and retains the ordinary biology schedule.
It rejects the ordinary Arena checkpoint adapter, whose source fingerprint does
not describe unloaded study code. Diagnostics must unwrap the physical models
with `rkc_physical_processes`; generic Arena accounting cannot infer this new
wrapper's ownership. The new `scripts/study_cellular_rkc_ecology.jl` compares a
fresh, common seeded dense closed initial condition across default, RKC16 and
fine512 mechanics, without splitting the physiological clock. The complete-trajectory
measurements below characterize this adapter. Stage refinement can change contact
history, funded work, births and future controllers; no exact-trajectory claim is
made for the new integrator.

## First complete-step gate

A fresh 192-founder closed constructive ecology in a 5³ periodic physical box
and 10³ transport grid ran ten common dt=0.1 steps. The controlled source pin
was 0a1c6a8 plus hashed study overlays, retained at
`/othersys/tmp/intmat-rkc-full-0a1c6a8-20261015`. Single warmed whole-trajectory
samples were default 210.457 ms, RKC16 191.747 ms, fine512 1,906.933 ms. This is
only a preliminary 1.098× result, not the 1.879× mechanical-fixture gain.
RKC allocates 42.825 MB versus default 19.293 MB over the interval.

Against mechanically fine512 with the same physiological outer clock, default
position RMS/max were 0.01468/0.1520, RKC 0.003354/0.03084. Paid-work error was
-0.05673 default, +0.03100 RKC. Bond Jaccard was 0.98927 default, 0.98217 RKC:
the candidate is not uniformly closer. Both retained all 192 stable IDs with
no births/deaths during this short interval, finite state, positive funds and
closed material/chemical-energy residuals below 3.1e-7. This is not a turnover or
long-term ecological equivalence gate. Repeated interleaved timing and finer
whole-trajectory reference convergence remain necessary before a performance
conclusion. The initial condition was constructed independently and identically
for every arm, not selected from a candidate trajectory.

## Refined complete-step result

The second immutable probe keeps the same core source and initial state, adds
fine256/512/1024 and RKC8/16/32, and alternates five warmed default/RKC16 timing
pairs. Median paired speedup is **1.1114×**; median times are 213.362 ms and
191.979 ms. The candidate allocates 2.22× as much memory in this case. The source,
all terminal observations and full logs remain at
`/othersys/tmp/intmat-rkc-refined-0a1c6a8-20261015`; compact evidence is in
[cellular-rkc-ecology-results.json](cellular-rkc-ecology-results.json).

All rows cover the same physical time and physiology schedule. The table gives
errors against fine1024; timings outside the paired default/RKC16 comparison
are single warmed observations.

| Method | Time (ms) | Position RMS | Max displacement error | Paid-work error | Resource relative L2 | Bond Jaccard |
|---|---:|---:|---:|---:|---:|---:|
| Default | 213.36 (median) | 0.01518 | 0.1564 | -0.05910 | 0.4502% | 0.98927 |
| RKC8 | 132.15 | 0.003682 | 0.03413 | +0.09347 | 0.3286% | 0.98217 |
| RKC16 | 191.98 (median) | 0.002873 | 0.02530 | +0.02862 | 0.3199% | 0.98217 |
| RKC32 | 280.13 | 0.001226 | 0.008711 | +0.02327 | 0.2599% | 0.98927 |
| Fine256 | see JSON | 0.001680 | 0.01717 | -0.007131 | 0.04397% | 1 |
| Fine512 | 1,909.09 | 0.000540 | 0.005613 | -0.002375 | 0.01435% | 1 |

Fine512→1024 is a useful uncertainty measurement, not a proof of the continuum
limit near contact-history discontinuities. RKC16 improves positions, work and
resource-field errors but slightly worsens bond-set agreement. RKC8 pays more
work than the reference despite better positions than default; RKC32 restores
default bond agreement but costs more. This is an accuracy/cost Pareto result,
not uniform dominance. All seven arms retained 192 IDs without births/deaths;
closed mass/energy residuals stayed below 7.1e-7, with finite positions and positive
funds. The separate 55-assertion composition suite exercises actual turnover,
secretion, exact unchanged-mechanics oracle behavior and checkpoint rejection;
that suite is not a long-term ecological equivalence test.

No production solver or numerical default was changed. A stronger integration
gate would need repeated crowded trajectories with turnover and an explicit
acceptable error target for mechanical work, chemistry and contact history.

The explicit RKC stage has no nonlinear solve residual: its private adapter
returns `NaN` in that inapplicable diagnostic slot. Acceptance uses stage travel
and true endpoint potential, not this slot. Physical position/fund finiteness and
the material/work ledgers are checked separately.
