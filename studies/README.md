# Studies

Run an investigation with `julia --project bin/interior.jl study NAME`.
These programs combine simulation laws, analysis and artistic observation; they
are separate from the application commands and TOML recipes. Inspect a program's
parameters before starting a long run. Most accept an output path as their first
argument. `verify` writes bounded scene verification results; `compositions`
prints coupling comparisons. Default outputs resolve from the repository root.
