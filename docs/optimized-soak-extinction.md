# Optimized continuation: population loss at the finite-funding boundary

Read-only inspection on 15 October 2026 found no living cells in the optimized
broad-condition continuation. The Julia process was still running. Its last
populated sampled observation and first empty observation bound the transition:

| Observation | Step | Model time | Living cells | Cumulative precursor supply A / B |
| --- | ---: | ---: | ---: | ---: |
| Last populated sample | 5,204,150 | 520,415.007755 | 214 | 999,193.588 / 999,193.270 |
| First empty sample | 5,220,815 | 522,081.507780 | 0 | 1,000,000.002061 / 1,000,000.000359 |

These samples were recorded at 14:07:32.843 and 14:07:52.835 UTC. Every subsequent
sample available at inspection remained empty, through model time 997,741.014867.
Final cumulative counters were 42,106 births and 42,170 deaths, consistent with
the original 64-founder population. This is population extinction in the
recorded model, not termination of its runner.

Each precursor has a declared finite source budget of one million mass units.
The nominal rate density `.0025` over area `32 × 24` supplies 1.92 mass units per
precursor per model time, placing the nominal funding boundary at 520,833.333.
That lies inside the observed transition interval. The supplied counters then
plateau; their milliscale excess over the nominal totals is retained rather than
rounded away in the evidence.

Inspection of exact local copies of immutable remote archives confirms that the
source funding was effectively spent. The last populated archive, at time
513,798.507656, held 217 cells and about 13,510 mass units of remaining funding
per precursor. The first empty archive, at time 624,070.509299, retained source
budgets of only `2.2197e-42` and `2.2477e-42`. Those amounts are positive, not
literally zero. Dissolved precursor masses `.182003` and `.206381` also remained;
the observation does not mean every molecule of chemical fuel vanished.

The exact ordering of funding completion and the final deaths is unresolved
within the sampling interval. The evidence locates this population loss at its
declared finite-funding boundary; it does not establish a limit on possible
organization or persistence under other conditions.

Both archived states were preserved locally, SHA256-checked against the remote
files, and loaded without advancing under the exact pinned implementation.
All 31 loaded physics hashes matched; resource/grid aliases survived strict
restoration. The original parent and its preserved copy both retain SHA256
`d7f81c4cbbf21b05ce1918a142764eb0d22e1d6b68565c012c66cd759a6f01de`.
All 115 files in the remote source manifest also matched. Nothing remote was
written, stopped, restarted, refueled or reseeded.

[Selected observations, source maps and checkpoint hashes](optimized-soak-extinction.json)
accompany this note. Binary evidence remains under
`outputs/optimized-soak-extinction-20261015`; the original run and parent remain
under `/othersys/tmp/intmat-soaks-optimized-20261015`. The original broad and
mosaic soaks are separate trajectories, as described in
[the soak protocol](cellular-soaks.md).
