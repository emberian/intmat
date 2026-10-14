# A running active-compartment arena

The arena couples externally funded nutrient exchange, heritable recurrent organisms,
a remodeling and signaling medium, and implicit conserved phase separation. Its
timestep is simulation time; UI and checkpoint cadences use wall-clock seconds.
The default profile has 96² unit-sized cells and 1,000 founders. The optional 32³
profile uses the same physical cell spacing. Changing resolution and changing
physical extent are independent decisions; profile coefficients are not rescaled
silently.

```sh
julia --project bin/interior.jl arena --continuous --output outputs/arena/active-2d
julia --project bin/interior.jl arena --profile experiments/arena/active-3d.toml --steps 100 --output outputs/arena/active-3d
julia --project bin/interior.jl arena --resume outputs/arena/active-2d/checkpoint-1.jls --continuous --output outputs/arena/active-2d
```

`latest.json` contains downsampled central field slices, a bounded organism sample,
lineage and action metrics, full population size, and material budgets. It is an
observation; `checkpoint-N.jls` stores all future-determining state, including
controller weights, plastic marks, organism IDs, deterministic mutation counters,
full 3D fields, coupling aliases, simulation clock and source ledgers. Checkpoints
use Julia Serialization and are trusted local files, requiring the exact Julia
version and loaded-module source hashes by default. Fingerprints are captured when
Interior loads, so editing files while a process runs cannot relabel its compiled
implementation. They are not a portable archival format.
`restore_arena(path; allow_source_change=true)` explicitly permits an implementation
migration, without promising identical subsequent evolution. `checkpoint.json`
names the latest completed slot. Ctrl-C and the stop command save before exit.

A fresh output directory must be empty. The runner creates `arena-owner.json`
and acquires an atomic `.arena-writer` lock recording its PID. Existing owned
output requires explicit `resume=true` (the CLI forwards `--resume`). A second
live writer is refused. A recognized stale lock is recovered only after its PID
is confirmed absent; unknown lock contents require manual inspection.

The runner owns only `latest.json`, `control.json`, `checkpoint.json`, and its fixed
ring of checkpoint slots within the chosen output directory. These paths are
atomically replaced. Choose a dedicated directory; it never removes other files
or grows an unbounded snapshot history. Three slots are retained by default.
To archive a particular state, copy its checkpoint to an independently named path.
The existing `record_state` format remains useful for portable field inspection.

A local UI/server can atomically write `command.json`:

```json
{"id":"pause-2026-10-14-1","action":"pause"}
```

Supported actions are `pause`, `resume`, `stop`, `checkpoint`, and `set`. Each command
needs a unique ID; a failed command is acknowledged and consumes its ID. The
acknowledgment appears in `control.json`. Setting accepts only `dt`, `source_rate`,
`sink_rate`, `snapshot_seconds`, and `checkpoint_seconds`, for example:

```json
{"id":"feed-2","action":"set","values":{"source_rate":0.015}}
```

Commands are processed between complete physics steps, including while paused.
Negative/nonfinite settings, unsupported keys and unstable timesteps are rejected.
A resumed checkpoint retains its paused/stopped flags; issue resume for a paused
arena. The CLI clears the stopped flag when explicitly resuming; programmatic restart
sets `arena.stopped=false` before `run_arena!`.

The chemostat supplies inventory per cell and removes dissolved nutrient with an
exponential outflow fraction. Both ledgers record actual Float32 inventory changes,
including rounding. The displayed material residual compares grid nutrient,
deposit and organism storage/waste with initial inventory plus supply minus
outflow. Chemical work input and dissipated heat are a separate ecology ledger.
Founder storage is included in initial inventory. No automatic reseeding occurs
on extinction: a fresh run is a distinct experiment with its own seed and ledger.

The arena constructor seeds low Fourier modes rather than the phase process's
standalone white noise initializer. `initial_modes` is the maximum absolute wave
number, `initial_mode_count` defaults to 12, and `initial_amplitude` defaults to
0.12 (the standard deviation of each mode's coefficient, not total field RMS).
The constructor recenters that field to `[phase].mean`, default zero; `[phase].noise`
pertains to the standalone constructor and is superseded by the arena mode field.
Existing defaults preserve the original arena initialization exactly.

An implementation migration requires an explicit CLI opt-in:

```sh
julia --project bin/interior.jl arena --resume outputs/arena/active-2d/checkpoint-1.jls --allow-source-change --continuous --output outputs/arena/active-2d
```

This permits changed source hashes while keeping the checkpoint's saved fields,
controllers, coupling aliases, clocks and ledgers. The profile's initializer is
not rerun during restoration. Source validation remains strict without this flag;
Julia-version and checkpoint-schema checks remain strict with it.

The opt-in `model = "metabolic"` profile composes a `ResourcePool`,
`MetabolicEcology`, `BoundaryMaintenance`, and phase separation. The boundary
wrapper advances its referenced `ActiveMedium` exactly once. The original neural
constructor remains the default. Start the richer arena alongside an existing
run in a fresh output directory:

```sh
julia --project bin/interior.jl arena --profile experiments/arena/metabolic-2d.toml --continuous --output outputs/arena/metabolic-2d
```

The profile uses 256 founders on 64² cells with capacity 3,000. Two independent
chemical substrates begin in complementary smooth periodic niches and feed through
the same masks; these imposed opportunities do not establish emergent specialization.
Sources and sinks have finite total mass budgets, separate from their per-cell-time
rates. Exhausting a source budget stops its feed. Chemical potentials fund work;
structural nutrient and deposit provide material without becoming fresh chemical fuel.
The profile's lower-potential third species begins empty and can receive exported byproduct. Species, potentials,
initial levels, diffusion, source/sink rates and budgets are configured in
`[resources]`; `spatial_niches=false` selects uniform initial levels and masks.

Snapshots publish every resource species, separate boundary `film`, pool mass and
potential ledgers under `metrics.resources`, and boundary contributed material,
spent work, and recycled material. Organism `metabolism` contains the inherited
species allocation simplex, assimilation, efficiency and cooperation traits.
`boundary_exposure` is local neighboring film, an exposure measure rather than a
fitness result. The material residual includes world fields, organism cargo/waste,
chemical species and film once each despite their shared references. Exact local
restart preserves those shared references and all finite source budgets; source
hash checks cover every Julia module, including the added processes.

Metabolic snapshots report `chemical_energy_accounted` (remaining chemical
potential, stored organism work, dissipated controller heat and metabolic heat)
and `chemical_energy_error` relative to initial energy plus external supply minus
outflow. Boundary work is already dissipated controller heat. The inherited
`byproduct_fraction` exports a share of unassimilated material to a lower-potential
species, permitting measured local reuse without restoring spent potential.
Resource and film slices use material per cell; organism boundary exposure is the
mean over its cell and adjacent cells. Live scalar chemostat feed controls apply
to the original model; metabolic per-species feed rates and finite budgets are
set in the profile.

`experiments/arena/metabolic-3d.toml` is an optional 32³ prototype with physical
extent 16³ (cell spacing 0.5), 1,024 founders and capacity 12,000. It retains the
2D profile's chemical rules, rates and finite budgets. A bounded CPU run on
2026-10-14 completed 64 steps to time 6.400000095 in 13.30 seconds measured around
`run_arena!`, including its first-call compilation and snapshot/checkpoint writes;
package loading and construction were outside that timing. Population reached
1,536, all inspected chemical/film/phase/controller arrays were finite, material
residual was −0.00218 and energy residual −0.00243. Actual central slices, full
checkpoint and source fingerprints are in
`outputs/arena/metabolic-3d-smoke-20261014`. This checks 3D composition and restart
plumbing over a short interval; it does not establish long-term diversity or
persistent compartments.

```sh
julia --project bin/interior.jl arena --profile experiments/arena/metabolic-3d.toml --steps 64 --output outputs/arena/my-metabolic-3d
```

## Develop beside a running arena

Loaded Julia modules keep evolving under their captured source hashes. Editing the checkout does not hot-reload that process. The two local running studies retain source snapshots at `outputs/arena/runtime-source-a2dd363` (original) and `outputs/arena/runtime-source-8765769` (metabolic). Both snapshots strictly restored their corresponding live checkpoints after later checkout edits. To recover a run, invoke its pinned `examples/arena.jl` with that snapshot as `--project`, its profile, the owned output directory, and the checkpoint named in `checkpoint.json`. Do not start a second writer while the existing process is alive. Reviewed migrations remain an explicit opt-in; ordinary development is not permission to reinterpret stored state.

New source snapshots must include `bin/` and `experiments/` alongside `src/`, `Project.toml` and `Manifest.toml`; include `studies/` when preserving study commands. Launch a new snapshot through its own `bin/interior.jl arena` with that snapshot as `--project`. The older snapshots above predate this layout and retain their original entry points.
