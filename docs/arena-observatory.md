# Arena observational recorder

`scripts/record_arena_observations.py` reads published HTTP snapshots without controlling the simulation. Each distinct observed step adds one compact JSONL metrics row per source. Every five minutes it also saves an immutable, uniquely named gzip of the complete published JSON. These are sampled observations, not complete simulation history or restart checkpoints. Earlier extinction times and trait trajectories cannot be recreated from them.

The live run started on 2026-10-15 UTC is in `outputs/observatory/20261015/live-054735-437165/`. Its `manifest.json` records PID, sources and parameters. It observes original arena port 8780 and metabolic arena port 8781 every 20 seconds. The metabolic arena was already extinct when observation began. Its simulation remains available.

Each metrics row preserves the server's full `metrics` object, model time, observed wall time, source generation time, actual population, captured agent count, captured scalar trait summaries, and captured genome signature counts/diversity/entropy. Agent-derived summaries describe only the published agents. `agents_complete` identifies truncation; the known snapshot default display cap is 3000, explicitly labeled as a source-code default because the HTTP schema does not advertise its effective cap.

`status.json` is atomically replaced with current health, latest step, age of the published snapshot, time of last distinct step, per-source errors, archive names and charged bytes. `errors.jsonl` preserves endpoint and parsing failures even after recovery. An unchanged snapshot does not duplicate metrics or archives. Stale snapshots remain observable through timestamps; an unavailable source does not discard the other source's observations.

The default budget is 6,000,000,000 bytes of new observations per recorder run, with 64 KiB reserved for status. At that limit the recorder exits with `budget_exhausted`; it never deletes old data. HTTP reads are limited to 64 MiB per snapshot. Signal SIGTERM or SIGINT requests a clean exit (an active HTTP read can take up to its timeout). Status and manifest record the process ID; stopping the recorder does not stop either arena.

Start a separate run using a fresh output directory:

```sh
python3 scripts/record_arena_observations.py \
  --source original=http://127.0.0.1:8780/latest.json \
  --source metabolic=http://127.0.0.1:8781/latest.json \
  --output outputs/observatory/NEW-RUN \
  --interval 20 --archive-interval 300 --budget-bytes 6000000000
```

Existing output directories are refused. `--polls 2 --interval 1` provides a bounded smoke run. Run status timestamps use UTC Unix seconds; local process logs are adjacent to the live output directory. No checkpoint copies are made.
