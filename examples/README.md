# Compatibility entry points

The supported application interface is `julia --project bin/interior.jl --help`.
Runnable programs live in `bin/cli/`, research programs in `studies/`, and editable
recipes in `experiments/{arena,materials,bigruns}/`.

These old Julia paths forward to the canonical programs. TOML paths are relative
symlinks to the single recipe source; edit the canonical recipe. They remain for
existing commands, tests and recorded provenance. New work should use the unified
CLI and canonical recipe paths. Existing output directories and checkpoints have
not moved.
