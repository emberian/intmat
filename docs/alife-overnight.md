# What the overnight arenas actually retained — 15 October 2026

There is sustained, mutating, resource-funded artificial life here: organisms reproduce and die, inherited controllers vary, and living populations remodel a persistent medium. The stronger claim—evolution produced a new functional organization, persistent ecological division of labor, or adaptive cooperation—has not yet been tested. The spatial phase structures and complementary chemical source masks have substantial imposed physics; attractive bands are not themselves evidence of evolved organization.

## Saved observations

Exact JSON copies of both live observations are preserved in `outputs/alife-20261015/{active,metabolic}-observed.json`. `summary.json` was generated with `scripts/analyze_arena_history.py`, which accepts observations, replay arrays, and the finite study's sampled results. Live process state was neither paused nor modified.

| Observation | Simulation time | Population | Reproduction / ancestry evidence |
| --- | ---: | ---: | --- |
| Original arena, saved live | 61,189.50 | 1,332 | 3.48 million cumulative births; maximum living generation 1,659; one surviving founder lineage |
| Metabolic 96², earliest replay frame | 183.50 | 2,811 | Maximum living generation 22; 48 surviving founder lineages |
| Metabolic 96², last replay frame | 247.50 | 3,115 | Maximum living generation 22; 41 surviving founder lineages |
| Metabolic 96², saved live | 411,704.01 | 0 | 743,520 births / 744,096 deaths, including the 576 initial founders |

The original arena's 1,332 observed inhabitants have 1,332 distinct inherited-controller signatures despite sharing one founder. Founder loss does **not** mean genotype identity. The metabolic replay similarly contains 2,811 distinct signatures among 2,811 organisms initially and 3,000 among its last bounded sample of 3,115. Signatures are truncated SHA256 hashes of inherited controller weights (`src/arena.jl`), excluding lifetime plastic state. Distinct signatures establish variation, not functional diversity, adaptation, or reproductive isolation.

The 16-frame metabolic replay covers only t=183.5–247.5. Its largest **sampled** population is 3,228 at t=228; that is not the overnight peak. Mean inherited efficiency changes from .65208 to .67266; mean cooperation changes from .23516 to .24017. The profile initializes those traits at .6 and .25. Mutation, selection, founder sorting, and drift can all contribute to these mean shifts. No common-garden fitness comparison separates them.

The arena retains three rotating checkpoints and one atomically replaced observation, with no longitudinal overnight journal. The current slots are all long after extinction. The exact extinction time, whole-run maximum population, and late trait trajectories cannot be recovered from these files. The current empty population's `max_generation=0` is an empty-set diagnostic, not evidence that generations never occurred. Future studies need append-only scalar history and periodic independent lineage/trait samples; merely leaving a live viewer open does not preserve them.

## The metabolic collapse

Both million-unit substrate budgets are exhausted. Remaining substrate masses are approximately 110.51 and 258.17 units, with zero remaining byproduct. Thus the model has undergone major finite-fuel turnover and eventual collapse; it was never an indefinitely fed steady-state experiment. Remaining structural nutrient (456,566 units) and waste (1.41 million units) are not fresh chemical potential. Residual fuel may be inaccessible or uneconomic to surviving phenotypes; exhaustion of the external supply is not proof that every final death was caused solely by inventory exhaustion.

Cumulative lower-potential byproduct uptake is 4,600.82 units, compared with 1,004,497.49 and 1,004,349.85 units of the two primary substrates. Recycling executes, but byproduct represents only .229% of measured uptake mass. Paid film received 6,578.56 units of material/work and later decayed/recycled. Movement consumed 1.121 million work units; remodeling consumed 1.167 million. Those are real modeled expenditures, not demonstrations of their benefit.

The saved material residual is −88.85 relative to approximately two million externally supplied units (−.00444%); the chemical energy residual is −9.32 relative to seven million supplied work units (−.000133%). Long Float32 integration has accumulated measured numerical drift. Neither ledger reports creation of enough energy to explain persistence.

The separately controlled, smaller 15-run affordability experiment retains its full trajectory. All arms reached the imposed capacity of 180, then went extinct with finite 2,000-unit supplies per substrate. First sampled extinction times were:

| Arm | Seed 11 | Seed 29 | Seed 47 |
| --- | ---: | ---: | ---: |
| niches | 592 | 568 | 840 |
| uniform | 1,048 | 464 | 432 |
| withhold | 400 | 520 | 416 |
| no_mutation | 384 | 472 | 440 |
| no_crossfeed | 504 | 448 | 504 |

The sample interval is eight simulation units. These comparisons show variable finite-budget ecology; three seeds do not establish an advantage to film, exported byproduct, or niche geometry. The large overnight arena has different physical extent, initialization and budgets and should not be pooled with those arms.

## A stronger next experiment

The next claim worth trying to establish is **heritable adaptation to resource geometry**, before adaptive cooperation. Run a separately named, bounded, continuously funded selection study with a source budget explicitly larger than the maximum possible feed during its declared horizon. Keep the live historical arenas intact. Use paired seeds for complementary niches versus uniform sources, each with mutation enabled versus disabled. Retain append-only scalar samples, genotype/trait distributions and periodic full checkpoints. Track supply depletion, births/deaths, population, resource uptake, founder effective count, controller and metabolic variance, patch allocation, work expenditures and spatial mixing.

After at least ten founder lifespans of sustained replacement, transplant randomly sampled descendants and the exact initial founders into fresh equal-resource common gardens. Reset age, energy, storage, plastic state and positions to matched distributions; preserve inherited controller weights and metabolic traits. Compare equal-number mixed populations for five lifespans with mutation disabled during the assay. Reciprocal home/away gardens separate general improvement from adaptation to geometry; randomized seed/position assignments and independent repeats reduce founder/place effects. Primary outcome: change in descendant-to-founder log frequency, with extinction handled as an outcome rather than dropped. Secondary outcome: resource-derived reproductive work per input chemical potential. A reproducible positive descendant advantage, especially a home-versus-away interaction absent in no-mutation controls, would be substantially stronger evidence than mean trait drift or unique genotype counts.

Film cooperation should then receive its own mixed cooperator/free-rider invasion assay and a matched intervention removing film's physical benefit while preserving its cost. Whole-population withholding only shows the system without maintenance; it cannot establish resistance to free riders. Likewise, no-inhabitant phase controls are necessary before calling a physical pattern organism-built.
