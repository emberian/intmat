# Native workbench control, version1

A fresh resident Julia worker owns two initially identical physical states. The
`off` branch leaves native sensory slots16–18 zero; `sensing` uses the profile's
explicit receptor channels/scales. This is a declared law edit at branch creation.
The worker starts paused and never advances without a finite `step` command. It
does not resume historical checkpoints, execute browser-supplied code, load arbitrary code or mutate a
running step or run a background soak.

```sh
DEVICE_BACKEND=cpu WORKBENCH_SIDE=8 julia --project=. \
  scripts/native_workbench_worker.jl \
  experiments/native/constructed-garden-material-sensing-64.toml FRESH_DIRECTORY
```

Native AMD uses the existing pinned environment and host resource wrapper in
`hbox-runtime.md`, including the compiler workaround when needed. The script's
network interface is intentionally empty. A separate loopback, same-origin HTTP
gateway writes only one in-flight `command.json` atomically. The worker reads at
most4096bytes. Expose observation JSON/binary files through the gateway, never
checkpoint `.jls` files. Public Pages can display recordings without this worker.

## Request and acknowledgement

```json
{"id":"a-unique-command-id","action":"step","branch":"both","steps":100,"requested_at":1792110000}
```

IDs contain1–80 ASCII letters, digits, underscores or hyphens. Allowed keys are
`id`, `action`, `branch`, `steps`, `scales`, `rates`, `requested_at`. Branch is an existing `off`, `sensing`, `fork1` or `fork2`; `both` (default)
always means the original `off`/`sensing` pair. `steps` is a positive integer only for the step action; booleans,
unknown fields and arbitrary expressions are rejected. The `receptor_scales` action
also accepts `scales`: exactly three positive finite numbers representable as
positive finite Float32 concentrations, with `branch: a single receptor-enabled branch`. Other actions
reject that field.

| Action | Effect at a complete step boundary |
|---|---|
| `describe` | Acknowledge current branches; physical configuration is in `worker.json` |
| `observe` | Publish immutable actual blockmean fields and cell state |
| `step` | Advance selected branches by the requested bounded count; publish actual observations afterward |
| `pause` | Acknowledge paused state; all branches are already paused between finite commands |
| `checkpoint` | Write a strict complete native checkpoint with source and content hash |
| `receptor_scales` | Apply three half-response scales to the sensing branch, preserving physical state and recording a timed intervention |
| `material_rates` | Change catalyst decay, food release and raw diffusion in an explicit fork, preserving state and the timestep |
| `fork` | Strictly clone one current complete branch into an independent generated child, preserving the immutable parent checkpoint |
| `stop` | Stop accepting new work and let the worker exit |

`control.json` contains `id`, `ok`, `accepted`, `status` (`completed`, `rejected`
or `failed`), `message`, `sequence`, and actual branch clocks and pointers. An
acknowledgement appears after completion or failure, never merely after a request
was queued. Duplicate identical IDs return the original acknowledgement without
repeating steps. Reusing an ID with different content is rejected. This guarantee
is scoped to the resident session; a crashed worker cannot be restarted over its
output directory. `commands.jsonl` preserves requests and acknowledgements.

Preflight errors do not change physics. A runtime failure may leave a partial
state; that branch is flagged failed and cannot be observed or checkpointed. A
wall deadline can interrupt a requested pair between complete steps; the
acknowledgement reports actual clocks, not the requested target. Previously
published observations and checkpoints remain immutable. Pause does not promise
interrupting a command in mid-kernel: keep commands small enough for interaction.

## Bounds and publication

The CLI publishes `worker.json` with `status: "starting"` before loading Interior
and compiling the native factory. This initial manifest has no branch state or
limits yet; the gateway and browser reject/disable commands until `ready`. Its
fresh directory can be adopted only with the startup token from the same process,
before any state exists. Initialization counts toward the wall deadline; the idle
clock starts after initial observations are complete. Initialization failures
caught by the CLI publish a stopped reason and preserve any existing evidence.

After initialization, `worker.json` reports `status` (`ready`, `running`, `stopped`), `stop_reason`, source
pins, recipe, receptor parameters, step size, branches and finite `limits`. Driver,
factory and observation-helper hashes are frozen when the worker loads, so later
source edits cannot relabel a resident process. Limits are:
`max_steps` per branch4096, `max_steps_per_command`256, `max_model_time`64,
`wall_seconds`900, `idle_seconds`120, `artifact_bytes`512000000 and
`max_commands`128. Display side defaults8 and must divide the simulation grid.
The Julia constructor accepts other explicit positive finite limits for tests or
a separately bounded launch; browser commands cannot increase them. Device and
host memory controls remain the launcher's responsibility.

Each branch has `branches/NAME/observatory.json`, schema
`native-garden-observatory-1`. Entries point to immutable
`frames/00001/frame.json` and `state.bin`, with the existing
`native-garden-observatory-frame-1` schema. Arrays are column-major with explicit
shape/type/offset and little/big endianness; the binary has a SHA256 hash. All nine
existing habitat channels are Float32 blockmean concentration observations,
including dimensionless permeability. Cells, bonds, identities, hidden state,
actions and final-step secretion are separate arrays; `material_inputs` adds the
three actual receptor values. Display interpolation never becomes simulation
resolution. The existing native observatory renderer can consume these indices.

Strict checkpoints live in `checkpoints/`, and observation export is never used
as restart state. The worker keeps native state resident, transfers only display
blockmeans plus cell records, and synchronizes observation after complete steps.
There is no racing GPU write during export. Artifact reservations precede new
files; command count, step/model-time caps, idle and wall deadlines bound the
session independently. Stopping does not automatically create a final checkpoint:
request one explicitly at a desired complete boundary.

## Exercised vertical slice

The [16³ miniature profile](../experiments/native/constructed-garden-workbench-16.toml)
uses four founders and matching per-volume concentrations/source funding; it
changes domain size and uses a study timestep.01. A CPU integration test covers
paired state identity, bounded steps, rejected commands, duplicate acknowledgements,
strict checkpoint restoration, actual viewer-compatible arrays, file hashes,
artifact and total-step limits, pause and stop. A separate blockmean reference
checks both inventory precisions against direct concentration averages.

The native RX 6750 XT protocol gate executed nine commands: describe, paired16-step
advance, checkpoint, pause, observation, sensing-only8-step advance, checkpoint,
rejected over-budget257-step request, stop. The final clocks are16 and24 steps;
six immutable8³ observations have verified binary hashes. Final checkpoint
material error≤1.78e−15 and potential≤3.11e−15. The first paired advance took12.67s
including compilation; the warm8-step command acknowledged in.151s including
polling/export. These are protocol latencies, not sustained throughput measures.
No worker remains running from this test.

[Compact native evidence](../benchmarks/2026-10-16-native-workbench-protocol.json)
records the exact source pins and bounded host controls. Full protocol report,
logs and source list are in `outputs/native-workbench-20261016`; the immutable
host stage is `/othersys/tmp/intmat-probe-20261014/workbench-v2-20261016`. The
first failed stage is retained. It exposed a pre-existing capture helper issue:
Float32 tuple accumulators widened heterogeneously when reading Float64 inventory
fields, producing unsupported dynamic tuple operations on AMD. Uniform Float64
accumulation now rounds only the final Float32 display observations. This repair
changes observation rounding, never physical simulation state. The gate pins its
own source, independent of later physiological changes in concurrent work.

The native v2 gate predates the metadata-only loaded-hash correction; its exact
executed driver is retained as `outputs/native-workbench-20261016/worker-v2-executed.jl`.
The repaired manifest's four loaded-pin assertions pass independently. Native
physics and command semantics are unchanged by this correction.


## Browser control gate

The loopback gateway serves the paired workbench and each branch through the
existing native observatory viewer. Embedded headers are compact; status and
errors remain visible. Standalone observatories retain their original layout.

```sh
python3 scripts/serve_native_workbench.py FRESH_DIRECTORY --port 8807
```

An isolated browser controlled a fresh CPU pair from the frozen source in
`outputs/native-workbench-ui-source-v2-20261016`: advance both 16 steps, advance
sensing 8 more, checkpoint sensing, stop. Actual viewer clocks ended at 16/24
steps with four cells each; each 8³ field export contained 4,608 values, and each
hidden-state array contained 32 values. Checkpoint downloads returned HTTP404,
257-step requests HTTP400, and there were no browser or observer errors. This
checks end-to-end control and export, not ecological usefulness.

Full source hashes, command acknowledgements and checkpoints remain in the source
stage and `outputs/native-workbench-ui-v2-20261016`. Browser results and reviewed
screenshots are in `outputs/native-workbench-ui-review-20261016`. The worker was
stopped; port 8807 serves the preserved observations.

`test_native_workbench_browser.cjs` exercises a fresh live pair.
`test_native_workbench_display.cjs` is read-only: it uses an explicitly synthetic
starting-status fixture, then displays the actual stopped session. This checks
disabled controls and compact embedded status without another simulation.

The final isolated Julia module gate passed 179 checks; six HTTP gateway checks
and five lightweight startup ownership/failure checks also pass.
[Compact browser evidence](../benchmarks/2026-10-16-native-workbench-browser.json)
keeps the actual run source pins distinct from the later startup/layout-only gates.


## Paused receptor intervention

```json
{"id":"scale-study-1","action":"receptor_scales","branch":"sensing","scales":[0.5,0.03,0.006]}
```

The three channel identities remain fixed at creation. The miniature profile uses
raw total, bound total and catalyst total concentrations. Each sensor response is
`c / (scale + c)`; a smaller scale increases response at the same concentration.
The worker validates all three values and constructs a same-type habitat wrapper
with the existing cells, fields, inventories and controller state. One wrapper
assignment applies the tuple. No simulation step or receptor recomputation occurs
until the next explicit advance; recorded inputs remain those of the last step.

Each application preserves `-before` and `-after` complete checkpoints and actual
observations at the same physical time. Their physical arrays agree; effective
parameters differ. Checkpoint `workbench` metadata and the command journal carry
the intervention ID, step, time, old and new scales; the strict native payload
carries the effective parameters used for restoration. A full artifact reserve is
checked before publication. If publication fails after application, the failed
acknowledgement and current manifest report the applied configuration; this is
an atomic parameter tuple change, not a transactional filesystem commit.

Observation identities use metadata hashes so a configuration-only change is not
mistaken for an unchanged frame. Optional `observation_event` distinguishes same-
time before/after frames in the existing viewer. Its stamp shows the intervention
phase; historical exports without this field retain their existing behavior.

For a **new** resident interactive session, launch-time environment variables
`WORKBENCH_WALL_SECONDS=3600 WORKBENCH_IDLE_SECONDS=900` can keep Julia warm across
many short commands. They change only that fresh worker's finite bounds; elapsed
initialization counts toward wall time, and no simulation advances while idle.
The gateway has no deadline-extension command. Existing closed runs are unchanged.


### Exercised scale change

The CPU browser gate changed `[4, .02, .004]` to `[.5, .03, .006]` at step16,
time.16. Before/after binary observations are identical, while their metadata
identities and effective parameters differ. Both same-time frames were reachable
with the viewer's best-resolution mode, and the intervention phase was visible.
An additional eight sensing steps then completed; the pair stopped at16/24 steps.
Strict restoration retains the new parameters. The isolated Julia suite passed
204 checks; the HTTP gateway passed7. The native algorithm and GPU kernels were
not changed or duplicated by this adapter.

[Compact intervention evidence](../benchmarks/2026-10-16-native-receptor-intervention.json)
pins the fresh source, event, binary identities and completed checkpoint. Full
browser data and reviewed images are in `outputs/native-workbench-receptors-review-20261016`.
The worker is stopped; port8808 serves those preserved actual observations.


## Fork a current completed state

```json
{"id":"fork-study-1","action":"fork","branch":"sensing"}
```

This uses the selected branch's **current completed state**, not the viewer's
recorded timeline position. Child names are generated internally as `fork1` and
`fork2`; four total branches is a hard maximum. No arbitrary path, external
checkpoint, source migration or failed partial state is accepted.

The worker saves and verifies a complete parent checkpoint, then uses the native
strict source-aware restore into independent arrays. Time, world seed, attempted
steps, complete controller state, inventories, paid source funding and original
accounting baseline are retained. Parent provenance and inherited intervention
history are included in the child's manifest and checkpoint. The parent checkpoint
is immutable; either live future can subsequently advance without modifying it.

Children share the existing session wall deadline, idle policy, artifact cap and
command cap. The absolute step and model-time limits are unchanged; forking does
not buy additional lifetime. A child from a receptor-disabled parent stays
receptor-disabled. Receptor-enabled children support the existing scale command;
its channels and concrete native parameter types remain fixed.

The browser selects an existing branch for commands and independent left/right
branches for comparison. A successful fork selects its parent on the left and its
child on the right. Commands are gated by the worker's advertised actions, so
older frozen paired workers remain usable without gaining unsupported forks.


The fresh CPU browser gate forked sensing at step 4, advanced parent and child
independently to step8, and obtained identical actual exported arrays. It applied
new receptor scales only to the child, selected the two futures in the browser,
created the fourth branch and verified refusal of a fifth. The parent checkpoint
hash stayed unchanged. All workers from this gate stopped normally. The Julia
composition gate passed233 checks, including exact full-state future parity,
changed-child divergence after equal-time advances, unchanged parent during a
child advance, inherited absolute step bounds and rejection of failed state.
Five additional ownership checks verified that every restored cell, environment,
habitat and feedstock array, plus the observation buffer, is independent of its
parent. The original wall deadline was unchanged across the browser session.
Eight HTTP tests passed. A read-only old-worker ready-status fixture verified
that existing frozen paired workers do not expose the new fork action.

[Compact fork evidence](../benchmarks/2026-10-16-native-workbench-fork.json) records
the source pins and checkpoint hash. Full browser results and the reviewed
comparison are in `outputs/native-workbench-fork-review-20261016`; port8810 serves
the stopped actual session. This is a control/independence result, not an ecology
or long-horizon organization claim.

## Inspect a numeric resident ID

The optional **Inside one cell** panel reads existing completed exports. It
follows the left/right comparison branches, offers only exact shared capture
times, and keeps same-time before/after intervention phases selectable. The
newest event at a shared time is the initial default. If the latest branch
captures have different times, the panel explicitly labels its older shared
comparison; it never substitutes a nearest frame or interpolates state.

It displays energy, biomass, radius, metabolites and the three material receptor
readings, with expandable hidden state, adaptation, controller commands, catalytic
expression and last-step secretion. Differences are right minus left. Effective
scales and intervention phase come from each frame's metadata. Controller commands
are not measured forces, flows or successful births. Inputs remain those of the
last completed physical step, so an immediate receptor-scale change can correctly
show changed scales and zero state differences.

A numeric ID must be resident in both selected captures. Missing or dead residents
are reported without zero filling; incompatible ancestry metadata suppresses
all differences. Even matching numeric IDs and ancestry do not establish that
descendants born after independent forks are the same individual. The UI calls
cross-branch correspondence unverified. Missing older optional fields remain
unavailable, rather than fabricated.

Optional energy dots cover up to eight exact shared capture times, using the
newest event at each time. Dots are actual observations, with gaps for missing
pairs and no connecting/interpolated curve. Full binaries are hash-checked, capped
at32MiB each, then discarded after extracting compact cell arrays. The cache is
bounded by20 captures and16MiB of cell records. This is read-only browser work;
there are no simulation, checkpoint or publication writes.

Nineteen decoding/comparison checks cover column-major arrays, exact times,
missing IDs/fields, ancestry mismatch, empty populations, nonfinite values and
invalid binary metadata. Nine HTTP tests pass. The read-only browser gate used
preserved `sensing`/`fork1` histories at steps4/8, selected both intervention
phases, plotted four actual energy dots, and verified explicit absent-ID and
empty-history fixtures. It made zero write requests and produced no JS errors.
Reviewed evidence is in `outputs/native-cell-inspector-review-20261016`.

[Compact inspector evidence](../benchmarks/2026-10-16-native-cell-inspector.json)
links the compared immutable frame hashes to the read-only review.

## Living construction adapter (bounded local session)

A fresh worker can select the existing moving/dividing, paid cohort-policy law:

```sh
JULIA_NUM_THREADS=2 DEVICE_BACKEND=cpu WORKBENCH_LIVING=1 WORKBENCH_SIDE=4 \
  julia --project=. scripts/native_workbench_worker.jl \
  experiments/native/living-patch-recruitment-16.toml FRESH_DIRECTORY
python3 scripts/serve_native_workbench.py FRESH_DIRECTORY --port 8811
```

This optional mode passes the existing `NativeLivingConstruction.LivingHabitat`
into the native step. It adds no numerical implementation. The original finite
session limits, paused command boundaries, loopback gateway and immutable export
protocol apply. Use short finite advances; the default 64 model-time cap allows
separate continued experiments within the original deadline. The initial pair
shares policy initialization and differs in declared material receptor access.
Fork a current completed branch, then change its receptor scales to compare a
causal change in what the learned policy observes. This edit preserves existing
learned weights, critic, hidden/adaptation state, delayed cohorts and private
history; it only affects future physical sensing steps.

Native physical checkpoints alone are insufficient for this law. Living-mode
checkpoints therefore also contain a source-pinned sidecar: all actor/critic and
cohort arrays, ID and invalidation records, cumulative paid work, pending births,
birth event arrays, original feed masks, policy options, private-memory mode,
inheritance law and policy step/noise clock. Fork restore allocates every sidecar
array independently and binds the new wrapper to the restored native habitat.
Native checkpoint state retains actual current feed masks, finite budgets and
all physical/controller arrays. No external restore or source migration is
accepted. The original branch clocks, source funding and session deadline apply
to the child.

Observations additionally export learned weights, critic, cumulative learning
work, pending cohort IDs and enqueue/maturity counts for active residents. These
are actual completed-state arrays; the existing viewer may ignore optional
arrays. `living` frame metadata and worker branch metadata identify policy
options and clocks. Complete checkpoint files remain inaccessible over HTTP.

The effective `config.living_policy` table is always explicit in the worker
manifest. Its optional recipe overrides are:

```toml
[living_policy]
inheritance_fraction = 0.5
horizon = 2.0
actor_rate = 0.003
noise_seed = 6201
sigma = 0.3
```

These are this workbench adapter's declared defaults, including the half-copy
inheritance used by the copy-law scene. They differ from the standalone study
runner's default zero inheritance. The policy acts throughout the finite session
model-time allowance; birth-transfer adjustment and private recurrent features
remain enabled. The stored checkpoint carries all effective policy options.
Changing this table requires a fresh session; receptor edits within a resident
session retain its declared inheritance and learning law.

### Living HTTP and viewer gate

A frozen-source four-founder 16³ fixture exercised the actual loopback gateway:
advance sensing two steps, fork it, change only the child's receptor scales,
advance the child two steps, then stop. The parent remained at step2 and the
child reached step 4; the parent's exported binary and immutable complete
checkpoint hashes remained unchanged. The receptor edit preserved identical
before/after binaries at the same time. The original off branch stayed at step0.

The embedded native viewer decoded all six optional living arrays directly from
the served, hash-verified binary: learned weights `[4,23,3]`, critic `[4,8]`,
Float64 paid work `[4]`, pending cohort IDs `[4,10]`, and enqueue/maturity counts
`[4]`. Values were finite, metadata lengths agreed, and the browser reported no
page errors. Complete checkpoint downloads returned HTTP404. This short fixture
checks the user-facing instrument; it establishes no ecological outcome or paid
learning progress. The owned worker and gateway both stopped. Existing local
servers were untouched.

[Compact living HTTP evidence](../benchmarks/2026-10-16-native-living-workbench-http.json)
records the frozen source, finite 600s wall/180s idle limits, commands, clocks and
hashes. Full receipts, actual checkpoint/observation state and reviewed comparison
remain in `outputs/native-living-workbench-http-20261016`.

## Fork-only material timescales and raw transport

A fresh worker advertises `material_rates`. Select a current complete branch,
**fork it**, then edit the generated child. The original `off` and `sensing`
branches reject physical rate edits, as does `both`. Changes are acknowledged
at a complete step boundary and affect the next finite advance.

```json
{"id":"material-study-1","action":"material_rates","branch":"fork1","rates":{"decay_rate":0.05,"desorption_rate":0.12,"raw_diffusion":0.006}}
```

The three required fields form a typed Float32 law vocabulary:

| Field | Meaning | Accepted range |
|---|---|---|
| `decay_rate` | Catalyst's exponential decay rate; lifetime `1/k` | 0…10 per model time |
| `desorption_rate` | Bound-food release rate; release clock `1/k` without rebinding | 0…10 per model time |
| `raw_diffusion` | Same diffusion coefficient for raw A and raw B | 0…min(1, diffusion stability ceiling) in length²/time |

Zero stops the named kinetic process. Catalyst loss can still reduce binding
capacity and release bound food through the existing accounted collapse law, even
with zero desorption. Positive inputs must remain positive finite
Float32 numbers. Booleans, negative/nonfinite values, underflow, missing or extra
fields, and rates attached to another action are rejected. The finite rate caps
bound this interactive vocabulary; they are not claims about all possible
material laws. Adsorption, binding capacity, chemical potentials, source budgets,
cell capacity and array schema are fixed. The diffusion edit changes both raw
coefficients together; a profile's other feed parameter roles remain unchanged.

For actual grid spacing, the raw ceiling is conservatively rounded down from
`.45 / (dt × sum(1/spacing²))`, capped at 1. The worker also calls the native
`check_dt` on the candidate coupled process, covering the cellular bound and
raw diffusion at the unchanged submitted Float32 timestep. A rate edit never
silently reduces `dt`, advances a step, redistributes a field, or changes an
inventory or funding ledger. Decay and sorption already use the existing native
exponential and closed-form finite-capacity laws.

The adapter replaces same-type feed/habitat parameter wrappers around the
existing arrays and rebinds the existing living wrapper. It does not reconstruct
or reset the actor, critic, cohort queues, private recurrent/adaptation state, noise clock,
birth records or physiological state. Complete before/after checkpoints and
observations carry the event, old/new rates, unchanged `dt`, native stability
limit and original fork provenance. Before/after binaries at the edit time are
identical; metadata and effective parameters differ. Native checkpoints include
both habitat rates and feed diffusion, so a later strict fork inherits the
edited law and independent complete state. Source identities remain pinned to
the actual code; parameter changes are declared provenance, not source edits.

The local interface displays proposed catalyst/release times and nominal
`R²/(6D)` transport time over the actual exchange radius. The latter is an
analytic diffusion scale, not a measured arrival time. Its physical comparison
table uses actual completed observations: total catalyst, bound food, cumulative
raw conversion and cumulative total work. Differences appear only when the
selected latest captures have exactly equal model times. Unequal clocks show an
explicit mismatch instead. The existing volume viewers retain the corresponding
actual fields and cells. Frozen older workers do not acquire new controls.

The CPU composition gate passes 48 material-control checks; the existing
workbench/receptor suite passes 238 and the gateway suite 10. A real loopback
browser forked a four-founder 8³ living fixture at step 4, applied rates 2 / 1 / .03,
then advanced both futures to step 12. It verified identical edit-time binaries,
unchanged parent checkpoint, actual unequal-clock refusal and equal-clock
material differences, native viewer decoding, event labels and HTTP rejection
of unsafe/original-branch edits. Desktop/mobile views were reviewed at the
actual final 12/12 steps; the final review made zero write requests.

The same frozen physical adapter passed the bounded RX 6750 XT gate. At shared
time .24 the catalyst totals were .035610599 versus .028642580; these are short
control-fixture results, not ecological claims. Final material residual was 0
and maximum potential residual 1.45e−14. The edited strict checkpoint restored
and retained complete private state; the parent's immutable checkpoint hash
stayed unchanged. Both actual sessions stopped; the local gateway and read-only review server exited.
The device gate used CPU 10–11, two Julia threads, 8 GiB and a 900s timeout.

[Compact material-control evidence](../benchmarks/2026-10-16-native-material-workbench.json)
links actual source pins, finite bounds and retained receipts in
`outputs/native-material-workbench-20261016`. Reproduce the host gate with
`scripts/check_native_material_workbench.jl PROFILE FRESH_OUTPUT`, or exercise
a fresh local worker through `scripts/test_native_material_workbench_browser.cjs
LOOPBACK_URL REVIEW_DIRECTORY`. Neither script changes an existing recording.
