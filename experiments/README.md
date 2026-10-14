# Experiment recipes

TOML files configure the supported `bin/interior.jl` interface:

- `arena/`: persistent living arenas, with checkpoint and observation cadences.
- `materials/`: bounded material simulations with optional particle controllers.
- `bigruns/`: explicit refinement experiments with stability and resource budgets.

Recipes are data; reusable simulation laws belong in `src/`, executable research
and artistic investigations in `studies/`, and command orchestration in `bin/cli/`.
