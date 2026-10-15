# Cellular overnight explorations

Two funded, resumable conditions give the new cellular substrate a longer horizon.
They use the same Julia implementation as the local arena. The cellular lifecycle,
contacts and resource pool currently run on CPU; a GPU on the host does not make
these jobs GPU simulations.

| Profile | Physical environment | Development and junctions |
| --- | --- | --- |
| `experiments/arena/cellular-soak-broad.toml` (seed 401) | 32 × 24 periodic domain, one broad complementary resource stripe; diffusivity 0.06 | Response 0.8, adaptation 6, junction conductance 0.8; plastic decay time 50 |
| `experiments/arena/cellular-soak-mosaic.toml` (seed 402) | Same domain, complementary cosine-product patches with wavelengths 16 × 8; diffusivity 0.02 | Response 1.2, adaptation 24, junction conductance 0.25; plastic decay time 200 |

Both start with 64 founders, physical exchange radius 1, spacing 0.5 and the same
metabolic machinery. `niche_waves = [kx, ky, kz]` specifies integer cycles around
the periodic domain; positive axes multiply their cosines. The default `[1,0,0]`
preserves the original stripe. Masks are complementary and each has mean one for
these resolved recipes. Contrast 0.9 leaves both resources locally available.

The broad landscape allows movement and connected bodies to encounter a long
chemical gradient. The mosaic offers shorter, intersecting gradients, slower
external mixing, and developmental memory that lasts across more encounters.
Slower junction exchange leaves room for intracellular differences to persist.
These are analytical design choices, not claims that a particular organization
will evolve. Different seeds and several changed parameters mean the pair is an
exploration of two conditions, not a causal ablation or fitness ranking.

## Supply and population headroom

Initial precursor concentration is 0.45 per species. Continuing supply is 0.0025
per species per unit area per unit time: 3.84 total mass per model time. At chemical
potential 4 and conversion efficiency 0.65, that provides at most 9.984 work per
time. Basal maintenance alone at 0.035 per biomass implies about 285 biomass in a
steady supply-only balance, or about 1,141 cells if every cell were at the minimum
living mass 0.25. Neural operation, motion and synthesis also consume work.

This calculation motivates a guard of 4,096 cells; it does not predict carrying
capacity. Initial and accumulated resource stores can support overshoots, bodies
can become large, and a time-varying population need not be in steady balance.
Inspect actual population, biomass, stored resource, births/deaths and blocked
births. Sustained blocking means the computational guard is affecting the run;
it is not evidence of higher fitness. There is no fixed-age death.

Each precursor has a finite source budget of one million mass, enough for roughly
520,833 model time at the nominal source rate. Waste export is budgeted separately
and counted by the open-system ledger. It prevents inert waste accumulating
indefinitely; the model has no waste toxicity reward or recycling into fuel.

## Execution and retention

The intended hbox allocation is two jobs, each with two Julia threads, one GC
thread, one BLAS thread, an 8 GiB memory limit, no swap, low CPU priority and a
12-hour wall bound. Use disjoint allowed CPU pairs 2–3 and 4–5 under the existing
shared build controls. Each job has a fresh directory and immutable source,
project and profile snapshot; the older running arenas remain separate.

The runtime lane owns deployment and records the actual host commands and paths.
Retain 20-second longitudinal metrics, 90-second rotating restart checkpoints and
10-minute full-state archives, with a 4 GB artifact budget per run. End by saving
a checkpoint; resume with the same pinned source and Julia environment. A budget
limit should pause collection or stop cleanly, never delete another run's files.
If a population becomes empty, preserve the terminal state; do not silently
reseed it and join that trajectory to the original run.

Useful first readings are spatial chemistry, body and bond organization,
fragmentation/reconnection, hereditary and developmental diversity, age and
lineage histories, resource turnover, and actual junction flux. Longitudinal
checkpoints allow later assays without changing the overnight trajectory. A
negative or quiet outcome concerns these founders and conditions; it cannot
establish that the substrate lacks the corresponding potential.

## Running hbox jobs (15 October 2026)

Both 12-hour jobs run from the immutable source at
`/othersys/tmp/intmat-soaks-20261015/source`. Its thirty Julia source files match
commit `606ca2e`; the copied profiles, driver and launcher have independent SHA256
fingerprints in `source/provenance.json` and `environment.json`. The resolved
manifest retains the portable `Interior` path `.` and uses Julia 1.12.4. All new
package writes use the isolated `depot/`, with the existing verified probe depot
available as a read fallback.

| Condition | Output | Julia PID at launch | CPU affinity |
| --- | --- | --- | --- |
| Broad | `/othersys/tmp/intmat-soaks-20261015/broad` | 3736521 | 2,3 |
| Mosaic | `/othersys/tmp/intmat-soaks-20261015/mosaic` | 3736525 | 4,5 |

Both actual cgroups enforce 8 GiB memory, zero swap and nice 15. The driver saves
at its 12-hour wall deadline, with an external interrupt/kill guard 15 minutes
later. `status.json` gives the actual deadline and progress; `checkpoint.json`
names the latest completed restart slot. `metrics-SESSION-NNNN.jsonl` chunks keep
longitudinal metrics. `state-STEP-SESSION.jls` files retain complete future state
every ten minutes. `latest.json` is the bounded current observation. The 4 GB
artifact cap reserves room for final checkpoints and stops without deleting
history. A `STOP` file in a run directory requests a checkpointed stop between
steps; rename that operator-created file before intentionally resuming.

A short real hbox soak completed at step 415, then an independent Julia process
strictly restored it, checked resource/grid aliases and advanced to 416. This
validates restart machinery rather than any long-term ecological conclusion.

To resume an interrupted **broad** job on hbox, first confirm its PID is absent,
then read the current pointer and remaining wall budget:

```sh
root=/othersys/tmp/intmat-soaks-20261015
checkpoint=$(python3 -c 'import json;from pathlib import Path;p=Path("/othersys/tmp/intmat-soaks-20261015/broad");print(p/json.loads((p/"checkpoint.json").read_text())["file"])')
remaining=$(python3 -c 'import json,time;from pathlib import Path;p=Path("/othersys/tmp/intmat-soaks-20261015/broad/status.json");print(max(0,int(json.loads(p.read_text())["deadline_unix"]-time.time())))')
test "$remaining" -gt 0 && SOAK_DEPOT="$root/depot" "$root/source/scripts/launch_cellular_soak.sh" "$root/source" "$root/source/experiments/arena/cellular-soak-broad.toml" "$root/broad" 2,3 "$remaining" "$checkpoint"
```

For mosaic, replace both occurrences of `broad` in the Python paths, the profile
and output names, and use CPUs `4,5`. The launcher accepts an explicit seconds
argument for a separately authorized extension. Strict source validation remains
mandatory; there is no automatic source-migration or reseeding path. A live
writer is refused by the output lock. The existing local arenas and their
observers are separate and were left running.

### Recovery boundary for the currently pinned jobs

The two jobs launched above retain their original driver hash
`cba3c8f41102aee0ebcd8dd4a2aa774dbb70377d0961bc57f37c92682e452aad`.
Their normal wall-deadline and `STOP` shutdown paths save between completed steps.
That first driver also attempts a final checkpoint after an exception. If a job
fails or is interrupted during a physics step, do **not** automatically resume
its final checkpoint pointer: inspect the log and use a rotating checkpoint or
full-state archive completed before the error. Preserve the questionable state
for diagnosis. This caveat does not require restarting a normally running job.

The repository driver has since been hardened for future launches. It rejects a
resume profile that differs from the stored configuration, and tracks complete
step boundaries. A mid-step failure leaves the prior checkpoint and observation
untouched and writes `unsafe_in_memory_state=true` failure metadata; it does not
serialize the partially mutated state. The standalone fault regression is
`julia --project test/cellular_soak_faults.jl` (it deliberately replaces a method
inside that test process, so it is excluded from the shared test harness).
