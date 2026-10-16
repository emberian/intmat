# Per-cell program observatory

The observation-only recorder in `scripts/cell_program_observatory.jl` accepts
plain completed-state frames from a Julia mechanism. It never steps, chooses,
learns, or forecasts. A capture has actual `time`, `step`, and `cells`; each cell
has its stable numeric `id`. Additional named fields remain as supplied, with
unknown or unrecorded quantities absent or null rather than reconstructed.

The first mechanism contract is `NativePredictivePrograms`: species-resolved
raw/food/catalyst receptor readings; fast/slow raw memory; a forecast of the next
raw receptor reading; a named program and actual investment/allocation/reserve;
cumulative requested/funded secretion; and cumulative paid prediction updates.
Forecasts are receptor forecasts, not utility predictions. Program catalogs and
units belong in metadata. A program phase is shown only if actually exported.

Forecast origin/target step stamps link a prediction to an actual later reading.
Without an exact recorded target, the viewer reports that the target capture is
unavailable. It does not substitute a nearby frame or attach a delayed error to
the currently displayed decision. Counter differences require the same resident
in adjacent recorded captures and nondecreasing cumulative values; these describe
the interval, not the causal effect of one chosen program.

```julia
include("scripts/cell_program_observatory.jl")
using .CellProgramObservatory
recorder = Recorder(metadata; ids=[1,2], max_frames=128, max_cells=8,
                    max_bytes=2_000_000)
capture!(recorder, actual_observation)
write_observatory(recorder, "fresh-program-observatory.json")
```

ID selection is explicit and stays fixed. Empty selection records all current
residents subject to the cell bound. No slot is relabeled as a departed resident,
and no missing resident is zero filled. Reaching a declared frame/cell/byte bound
stops further recording with an explicit reason; it does not alter the mechanism
or its clock. The recorder copies plain observations so later simulation writes
cannot change recorded samples. Separate run/branch identity belongs in metadata;
numeric IDs alone never imply identity across independently evolving branches.

The published baseline packages three existing eight-unit runs from
`outputs/predictive-programs-baseline-v2`, frozen against
`outputs/predictive-programs-epoch-v2`. Each recording contains 40 captures from
0.02 through 7.82; the measured endpoint report reaches eight units. Predictive
and current-reading selectors have exactly matching captured actions, funding,
and resident outcomes. This demonstrates no prediction-control benefit. The
finite programs are authored scaffolding; the reusable observation contract can
also carry future learned-controller state.

Build the compact, read-only page without running Julia:

```sh
python3 scripts/build_cell_program_observatory.py \
  --recording outputs/predictive-programs-baseline-v2/predictive-observatory.json \
  --recording outputs/predictive-programs-baseline-v2/persistence-observatory.json \
  --recording outputs/predictive-programs-baseline-v2/balanced-observatory.json \
  --run-report outputs/predictive-programs-baseline-v2/report.json \
  --destination site/laboratory/cell-observatory
```

The optional endpoint report supplies the duration caption, never a synthetic
endpoint capture. Exact forecast matching uses `sensed_step`, not the completed
frame step; the sense hook precedes frame completion. The expandable reading
guide separates model scores, receptor forecasts, physical funding and paid
updates. Source records and hashes remain downloadable. The laboratory links to
this page; replay controls never start or modify a worker.

Validation: 21 Julia recorder checks, four Python packaging checks, the synthetic
exact-clock browser fixture, and 72 selections across all three actual runs.
Desktop and 390-pixel mobile browser checks found no page errors, write requests
or horizontal page overflow. The synthetic fixture is a decoding test only and
is not packaged as an experimental outcome. Compact receipt:
`benchmarks/2026-10-16-cell-program-observatory.json`.
