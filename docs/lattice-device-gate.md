# Native AMD gate for the fixed-voxel ecology

Prepared 2026-10-15. **AMD execution is pending an allocated window after the production long garden.** Local CPU plumbing is tested; it is not accelerator evidence. The gate uses the reviewed finite-genotype lattice law, not the moving-agent ecology, turnover/packing law, or experimental packed-bank variant. The model SHA256 is explicitly checked before running: `e6d6cab2b0624965a5e7ec66b7f9f9ce44a01f9ded5653db660de61368cde9b5` (reviewed baseline in commit `cfb9544`). Its three kernels and G×14×11 bank layout remain unchanged.

## Predeclared acceptance

`studies/lattice_device_gate.jl` writes `protocol.json`, exact gate/model source copies, Julia version, and loaded shared-core hashes **before** loading or launching AMD. It requests a real `AMDGPU.ROCBackend`, checks runtime availability, and disables scalar device indexing. It does not substitute CPU for an unavailable accelerator. Any assertion or compile/runtime error prevents a passing result and remains visible in the launch log and partial output.

1. On 2×2×2 sites, forced contact fixtures independently exercise paid division with a one-bit inherited mutation and paid identity-preserving hopping. Each must actually occur with positive recorded work; body/energy ledgers must remain valid. CPU and target identity/event arrays must match exactly.
2. Real complete CPU and target steps begin from the same source-pinned checkpoint on 8×8×1 and 4×4×4 domains, seed 29, four founder codes from a 64-code bank, initial occupancy 0.25. Advance 200 steps of dt=0.1 and compare every ten steps. Both fixtures must actually exhibit births, mutations, hopping, matrix inventory and paid remodeling. Neither an all-empty nor inert pass counts.
3. Cross-backend integer arrays are **exact**: genotype, founder lineage, generation, events, IDs and parents. Ticks and world time must match exactly. Every compared floating element must satisfy `abs(target−CPU) ≤ 1e-4 + 2e-4*abs(CPU)`. Environment, diffusion scratch, cells, hidden state, actions, compensated heat, paid work and the shared weight bank are compared separately. Tolerances are fixed before AMD outcomes; failed identities cannot be relabeled as numerical success.
4. Independently on each backend, all measured physical stores must be finite and nonnegative; body+matrix capacity excess ≤2e-6; hidden magnitude ≤1.00001. Material and accounted-energy drift must each stay within `2e-5*max(initial inventory,1)`. These are declared numerical allowances, not biological resources. CPU/target field errors and each ledger are retained.
5. Record both states through generic native recording. Target checkpoint restoration must reproduce every persistent field exactly; scratch is intentionally recomputed. Seventeen further target/target-restored complete steps must then match **all** fields including scratch exactly. Replay is same-backend; CPU/AMD floating replay is governed by item 3, not falsely claimed bit-identical.
6. Reject aliased environment/scratch ownership, a negative dt without time/tick/cell advance, overwriting an existing checkpoint, a mismatched law-source hash, and mismatched loaded shared-source hashes. These are actual attempted invalid operations, not merely checks for validator functions.

The output includes per-field errors at every observation, event/work evidence, native recordings, initial/target checkpoints, intentional corrupt-checkpoint fixtures and source pins. Only after all correctness fixtures pass is `gate.json` written with `passed=true`. A later throughput failure does not invalidate an already recorded numerical gate, but it prevents claiming the requested sustained benchmark succeeded.

## Sustained whole-step benchmark, after parity

The optional `--benchmark` path is restricted to AMD. It measures full 32³ and 64³ worlds (32,768 and 262,144 sites), dt=0.1, seed 29, 4,096 shared codes/eight initial founder codes. Initial occupancy is 0.65, initial body mass 0.7, precursor 2 and finite reserve 1,024 per food/site. Benchmark parameters use construction 0.006 and matrix decay 0.02 to retain an active dense workload; all other baseline parameters remain unchanged. This is an explicitly declared workload, not the original ecological pilot or an adaptation experiment.

Each specimen warms for 128 complete steps. Three timed repeats start from identical copies of that warmed persistent state. Each repeat runs batches of 128 complete steps until at least three synchronized wall seconds, capped at 128 batches (16,384 steps). Synchronize before and after each batch. Timing includes normal validation, diffusion, metabolism, contact events and bookkeeping; it excludes compilation, state construction/reset, host observation and recording. No kernel-only/controller-only rate is substituted. The duration cap is reported: if the batch limit prevents three seconds, `reached_minimum_duration=false` must be treated as insufficient sustained-duration evidence, not silently accepted.

At warm-up and every timed-repeat endpoint, at least half the sites must contain living bodies and the physical checks above must pass. Population, births/deaths, hopping, mutations, matrix, paid work, ledgers, elapsed seconds, steps, resident bytes and throughput are retained for each sample. Empty or failing populations abort instead of producing an impressive inactive rate. Each repeated measurement covers a finite-resource trajectory; the reset between repeats is an experimental operation outside timing, not an unrecorded energy source inside the model.

At 64³ the baseline arrays are approximately 48 MB plus a 2.5 MB bank; the host warm-state copy and transient host observations remain comfortably within an 8 GB job budget. One grid is benchmarked at a time. The code never selects larger grids or an unbounded timing loop. The enclosing launcher must still enforce a hard timeout, since a compile/runtime failure can exceed a model-step estimate.

## Launch when the GPU owner grants a window

Use the existing pinned Julia/AMDGPU environment and reviewed source copy on hbox; do not instantiate, precompile, replace shared environments or overlap the production long garden. No use of CPU8–9 or CPU10–11 is authorized by this plan. Reserve CPU20–21, at most two Julia threads, one precompile task, BLAS1, memory8 GB, swap0, nice15, and a 15-minute wall timeout through the designated swarm resource wrapper. Confirm its active controls and accelerator ownership before launching; no command below has been executed on AMD.

```sh
JULIA_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1 \
JULIA_PKG_PRECOMPILE_AUTO=0 SWARM_MEM_MAX=8G \
/home/hbox/bin/swarm-build taskset -c 20,21 nice -n 15 \
timeout --signal=TERM --kill-after=15s 900 \
"$LATTICE_JULIA" --startup-file=no --project="$LATTICE_AMD_PROJECT" \
studies/lattice_device_gate.jl amd "$LATTICE_FRESH_OUTPUT" --benchmark
```

`LATTICE_JULIA`, `LATTICE_AMD_PROJECT` and the fresh output path are launch-specific paths to the existing configured runtime and pinned project. Preserve its configured ROCm/depot environment. Verify the wrapper actually provides the requested no-swap and memory controls rather than assuming environment variable names imply enforcement. The gate only adds a study file and outputs; no production law or arena is migrated.

CPU plumbing command (small fixtures only):

```sh
JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 julia --project=. \
studies/lattice_device_gate.jl cpu outputs/my-lattice-device-cpu
```

GPU compiler incompatibility, excessive cross-device divergence, a ledger/identity failure, missing event coverage, or inadequate active sustained throughput are useful gate outcomes. Preserve their sources/logs and repair a separately versioned successor; do not weaken tolerances after seeing failures or quietly rerun a changed law under the same provenance. Adaptation, organized multicellularity and useful morphology remain separate questions even if this numerical/performance gate passes.

## Prepared-source verification

Final local CPU plumbing passed with the exact prepared gate source SHA256 `de586b0f57aa153ee78d30b562cce296568d0d7b61cdb61bd5ac1f81d74e1119`. Evidence: `outputs/alife-20261015/lattice-device-cpu-frozen/{protocol.json,gate.json}` and archived sources/native recordings/checkpoints. Both shapes passed all twenty observation comparisons, the paid-contact cases and strict continuation/failure cases. The two earlier CPU smoke directories remain historical. No AMD compilation, parity or performance outcome is available yet.
