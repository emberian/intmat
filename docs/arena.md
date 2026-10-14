# A running active-compartment arena

The arena couples externally funded nutrient exchange, heritable recurrent organisms,
a remodeling and signaling medium, and implicit conserved phase separation. Its
timestep is simulation time; UI and checkpoint cadences use wall-clock seconds.
The default profile has 96² unit-sized cells and 1,000 founders. The optional 32³
profile uses the same physical cell spacing. Changing resolution and changing
physical extent are independent decisions; profile coefficients are not rescaled
silently.

```sh
julia --project examples/arena.jl --continuous --output runs/arena
julia --project examples/arena.jl --profile examples/arena/active-3d.toml --steps 100 --output runs/arena-3d
julia --project examples/arena.jl --resume runs/arena/checkpoint-1.jls --continuous --output runs/arena
```

`latest.json` contains downsampled central field slices, a bounded organism sample,
lineage and action metrics, full population size, and material budgets. It is an
observation; `checkpoint-N.jls` stores all future-determining state, including
controller weights, plastic marks, organism IDs, deterministic mutation counters,
full 3D fields, coupling aliases, simulation clock and source ledgers. Checkpoints
use Julia Serialization and are trusted local files, requiring the exact Julia
version and source hashes by default. They are not a portable archival format.
`restore_arena(path; allow_source_change=true)` explicitly permits an implementation
migration, without promising identical subsequent evolution. `checkpoint.json`
names the latest completed slot. Ctrl-C and the stop command save before exit.

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
