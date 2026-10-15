# A constructed catalytic habitat

`CatalyticHabitat` gives growing cells another way to participate in their
material environment: they construct extracellular catalysts that make food
accessible. It composes with the existing cellular bodies, developmental
controllers, bonds, resource transport and finite work reservoirs. It is a new
synthetic ecology, not an assertion that the earlier substrate was incapable of
organization.

The habitat contains two inaccessible feedstocks and two enzyme materials.
Enzyme A releases soluble precursor A from feedstock A; enzyme B releases B.
Cells still need both soluble precursors to build biomass. Any nearby cell can
import the released food, including cells that did not contribute catalysts.
Persistent construction can therefore provide a local benefit, create dependence
on somebody else's chemistry, be exploited, or be abandoned.

Construction consumes equal amounts of both intracellular precursors and an
explicit amount of stored work. The precursor chemical potential is dissipated;
the enzymes have zero potential in the chemical-work ledger. Catalytic release
moves equal mass downhill from feedstock to soluble food, with any potential
difference dissipated. Catalysts are immobile material fields, decay into inert
waste, and are never recycled into fuel. Sub-ulp transfer residue is retained as
explicit inert mass, with its lost chemical potential accounted as heat. Actual
field changes determine the transfers.

This gives construction a chemical function and environmental memory. The
current enzyme material is not also a membrane, spring, obstacle or optical
organism boundary. Cells keep their existing soft-body mechanics and conservative
junction exchange. A rendering can show the catalysts as material, but should not
invent those additional physical functions.

## Regulation and inheritance

Each cell has a separate three-output secretion readout. It senses the existing
inherited recurrent state, local feedstock and enzyme concentrations, internal
precursors, stored work, body size and age. One output controls investment; two
competing logits divide catalytic construction capacity between A and B.
The readout genotype is inherited and sparsely mutated at division using stable
cell IDs. No developmental state is reset by the habitat. Regulation feeds back
through material, energy and the local food field into the original controller.

The default exponent of two on allocation is an explicit synthetic
capacity/specialization tradeoff: dividing capacity across both catalysts has
lower throughput than specializing. Material and work are charged only for
actual construction. There are no imposed producer/consumer classes, ancestry
gates, bonuses for group size, or privileged access to one's own secretions.
Founder readout biases make varied investments and allocations available; they
are priors, not evolved specialization.

Three processes own chemical fields: the feedstock pool, the soluble pool and
the habitat's enzyme material. The cells and habitat borrow the pools. Typical
order is `(feedstock, soluble, cells, habitat)`: cells import and develop, then
secrete; released food becomes available on the next cellular step. Optional
existing phase/medium stages can remain in the composition. Schedule the habitat
on every cellular step so new descendants inherit from the recorded parent.
Distinct inventory arrays are required; aliases across pools or world fields
are rejected before mutation. Arena ledgers follow borrowed owners once and sum
supply/removal potential across both pools.

## Running and seeing it

```sh
julia --project=. --threads=2 --gcthreads=1 bin/interior.jl arena \
  --profile experiments/arena/catalytic-mosaic.toml \
  --output outputs/my-catalytic-arena --continuous
julia --project=. --threads=2 --gcthreads=1 studies/catalytic_habitat.jl \
  experiments/arena/catalytic-mosaic.toml outputs/my-catalytic-study 3000
```

Use fresh output directories. The study records three matched conditions:
regulated construction, catalysis switched off with construction still paid,
and a fixed secretion readout. These diagnose machinery and exploratory
trajectories; they do not measure evolved adaptation or rank the possible
architectures. Complete restart states and actual observation frames are retained.
Pin source before long runs; existing cellular soaks are unchanged.

Snapshots contain `enzyme_a`/`enzyme_b` (constructed mass per physical volume),
`release_a`/`release_b` (current soluble-food release per volume per time), and
`bound_a`/`bound_b` alongside the ordinary soluble chemistry, physical cells and
bonds. Per-cell `catalysis` reports investment, allocation, actual last secretion
and a separate genotype signature. These expose the actual causal material
rather than an illustrative coating.

The first low-turnover exploration (catalytic rate 0.35, soluble bootstrap
concentration 0.12) constructed catalyst but released little food before its
founders dwindled. That motivates a regime with faster catalytic turnover (4),
a modest larger bootstrap (0.3), and secretion capacity 0.12. The initial run is
retained as evidence about that regime. Neither its decline nor a later viable
trajectory settles whether enduring collective organization can evolve here.

## First active exploration

A source-pinned, 3,000-step run of the faster-turnover mosaic recipe reached model
time 300. All arms used the same founders and seed; only the named mechanism
setting changed.

| Condition | Living cells at t=300 | Cumulative births / deaths | Food unlocked |
| --- | ---: | ---: | ---: |
| Regulated | 18 | 67 / 81 | 284.9 |
| Catalysis off, construction still paid | 0 | 29 / 61 | 0 |
| Fixed secretion readout | 14 | 56 / 74 | 210.9 |

The regulated trajectory's largest recorded catalyst inventory occurred at t=80:
78 living cells, 114 bonds and 3.238 units of constructed enzyme material. Its
later population and material inventory declined; this is not a demonstrated
self-maintaining collective. The observations establish a functioning paid
construction–food-release feedback. They do not establish evolved adaptation,
cooperation, causal superiority of regulation, or long-term persistence. There
was no population-cap blocking. Maximum absolute accounting residuals across
these arms were about 5.5e-6 material and 1.2e-5 chemical energy.

The complete run is under `outputs/catalytic-habitat-20261015-active`, with its
immutable Julia implementation in the adjacent `-active-source` directory.
`report.json` contains all configurations, source fingerprints and endpoint
metrics; each arm retains actual sampled snapshots and three complete
checkpoints. The earlier low-turnover and interrupted numerical-regression runs
are also retained. A tiny negative precursor roundoff in the latter prompted
strict donor-bounded deposition, not a physiological change.

`catalytic-habitat-diagnostic.mp4` shows 121 actual recorded states, without time
interpolation. The left field is constructed catalyst, the right accessible
food; amber/teal distinguish A/B, with fixed logarithmic exposure throughout.
Recorded disks and bonds overlay both fields. `regulated-material.png` uses the
existing cellular renderer on the t=80 state: physical planar positions/radii
are measured; lens height and lighting are optical staging. Its background shows
soluble chemistry, while the diagnostic explicitly shows the catalyst fields.

[Watch the recorded catalyst and food fields](https://emberian.github.io/intmat/animations/#catalytic-habitat-20261015).

## Complete-step cost

The recorded t80 state (78 cells, 114 bonds) took a median **1.445 ms** per complete
step on Apple M2, Julia 1.12.4, two threads and one GC thread. Forty warmed
fresh-copy samples allocated 587,600 bytes per step; copies were outside timing.
The habitat stage took 0.197 ms (13.7% of the full-step median), cellular ecology
0.980 ms, tuple validation 0.182 ms and both pools 0.079 ms. The next useful targets
are cellular allocations and repeated validation; secretion arithmetic is a
smaller share. This measurement uses the experiment's own source pin, independently
of later working-tree optimization, and does not establish throughput at larger
populations. [Samples and source fingerprints](../benchmarks/2026-10-15-catalytic-t80.json).
