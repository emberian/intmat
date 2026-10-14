# Active arena experiments

The experiment program compares the complete arena against three matched ablations: medium feedback disabled, within-lifetime controller plasticity disabled, and inherited mutation disabled. Every comparison uses the same initial random field and founders for each seed. Resources enter continuously; births, deaths, and habitat changes should therefore be measured after the initial transient, rather than treating a crowded final frame as evidence of a living ecology.

Initial composition is a conserved mean plus small random Fourier modes. This seeds competing spatial structures throughout the periodic domain. Cell spacing is one; the coarse phase model uses bulk 1, gradient coefficient 0.7, mobility 1, stabilization 2, and a step of 0.1–0.2. At 32² use modes through wave number 4; at 16³ use 2 so the shorter domain still seeds unstable, resolved wavelengths. A positive phase region is a diagnostic habitat proxy, not a definition of an organism.

The recorded measurements answer distinct questions:

- Late births and deaths establish whether turnover continues after startup. Population change and blocked births reveal collapse or capacity saturation.
- Surviving founders and the exponential Shannon entropy of founder frequencies distinguish numerous descendants from diverse ancestry. Founder diversity cannot increase without immigration, so its loss is expected and its persistence is informative.
- Genotype variance measures inherited controller differences; plastic weight drift measures changes within a lifetime. Disabling mutation does not remove founder differences. Disabling plasticity does not remove recurrent memory.
- Positive-phase occupancy, periodic connected components, largest-component share, and component-size variation measure habitat structure. They do not establish membrane individuality or replication.
- Resource input/output ledgers, signal variance, phase mean, and finite state checks distinguish maintained dynamics from numerical drift or unaccounted flux.

A short suite can expose collapse, saturation, inert controllers, and broken coupling. It cannot establish indefinite novelty or open-ended evolution. A useful follow-up is a longer replicate study, novel resource regimes, and independent behavioral probes, after a coarse configuration demonstrably sustains turnover.

Run `julia --project=. scripts/arena_experiments.jl` from the repository root. The program writes its seed-level measurements and comparison summaries under `outputs/arena-experiments-fairshare`; exact settings and duration belong with those measurements.

## What the first comparisons exposed

The original controller granted new chemical work on every uptake while retaining all uptake as reusable structural matter. Death returned that matter to dissolved nutrient. The first suite therefore measured an externally refueled recycling loop, rather than energetics limited by supplied feedstock. Its complete records are explicitly tagged `legacy-uptake-work-recycled-material` under `outputs/arena-experiments`.

That distinction changed the design. At step 0.1, uptake 0.7 and lifespan 50, both 32² replicates had continuing birth/death turnover at time 1000, but both sat at the 240-agent capacity. Seed 11 retained four founders; seed 29 retained one. They accumulated roughly a million blocked birth attempts each. A finer 48² grid with extent 24 also filled its 600-agent capacity by time 160. Reducing resource input below the observed viability threshold caused extinction. These observations motivated explicit catabolism: a fraction of uptake becomes inert waste, and only that consumed fuel supplies chemical work. Recycled structural matter can no longer provide unlimited work.

The first catabolic suite writes `catabolic-sequential-allocation-v1` results separately under `outputs/arena-experiments-catabolic`. These are distinct models, not a parameter comparison. A further correction allocates pre-step cell nutrient and erosion stock proportionally among colocated agents, eliminating the earlier advantage from feeding in persistent-ID order; this final model is tagged `catabolic-fuel-fairshare-v2` under `outputs/arena-experiments-fairshare`. The material budget includes inert waste; the energetic budget includes initial agent energy and `energy_yield × catabolized fuel`. Continuing evolution and biological realism remain questions for the measurements.

## Final fair allocation measurements

The final mechanism was measured on 32² cells with spacing 1, step 0.1, 32 founders, capacity 600, uptake rate 0.7, lifespan 80, source 0.023 and sink 0.01. Both runs reached time 1000, over twelve lifespans. The late window is time 500–1000.

| Seed | Late births / deaths | Final population | Late population range | Surviving founders | Blocked births |
| --- | --- | --- | --- | --- | --- |
| 11 | 8489 / 8429 | 383 | 295–391 | 1 | 0 |
| 29 | 9097 / 8934 | 398 | 221–407 | 3 | 0 |

This establishes persistent turnover below the configured capacity for these two seeds. Founder diversity contracts strongly. The positive habitat proxy becomes one connected component in each run. These histories did not count negative-phase components, so they cannot establish how many enclosed pores remain within a connected positive matrix, or whether compartment individuality is maintained. Future captures record both phase signs. The medium's paid building, erosion, signaling, and remodeling mechanisms are available for stronger boundary studies.

At matched time 160, reducing the step from 0.1 to 0.05 in seed 11 gave final populations 82 and 90 with three surviving founders each; phase variance changed from 0.6481 to 0.6420. Event counts and intermediate populations differ substantially, so this is a gross viability check rather than trajectory convergence. The periodic conserved phase solver has separate refinement tests.

A 16³ smoke run at time 160 retained 14 founders, with effective founder count 6.88 and 519 inhabited grid cells. Its population reached capacity 600 and produced 279963 blocked birth attempts. Its resource supply per inhabitant differs sharply from the 2D study; it verifies a working three-dimensional composition rather than superior ecology.

Across these records the largest absolute material ledger error was below 0.001, and catabolized-fuel work matched its declared yield exactly. Energy ledger roundoff remained small relative to accumulated work. All sampled states were finite. Compact tracked settings, source hashes, final metrics, timings, and summaries are in `docs/arena-experiment-results.json`; full local sample histories remain in the corresponding output directory.

For the matched feedback ablation, seed 11 at time 1000 had population 349, two surviving founders, and late births/deaths 7635/7567, with no blocked births. The complete model had population 383 and one surviving founder. Thus these data establish usable bidirectional mechanisms, but do not show that feedback improves diversity or is required for turnover. The ablation suppresses secretion/erosion, signal emission and phase remodeling while retaining sensing, transport, controller memory, and the same paid work budget; effort allocated to a suppressed output is still dissipated. Plasticity and inherited mutation remain enabled.

The arena constructor supports configurable mode count and per-mode amplitude, and honors `phase.mean`. The coarse lab and active boundaries profiles now seed twenty-four modes with coefficient standard deviation 0.08 around mean −0.12. The experiment harness instead normalizes its twenty-four coefficients by the square root of mode count and uses a separate seeded founder positions stream. These profiles therefore transfer operating parameters rather than replaying the measured histories exactly. Restoring an existing checkpoint preserves its saved initial condition and current state; it does not rerun the updated initializer.
