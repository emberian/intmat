# Finite metabolic comparisons

`scripts/metabolic_experiments.jl` supplies matched seeded comparisons of two chemical substrates and one lower-potential byproduct, heritable enzyme allocation, assimilation, efficiency, and funded cooperative film. This is an opt-in CPU experiment: the new resources and metabolism do not yet run on Metal. GPU execution of existing dense fields or batched controllers is a separate capability, and cannot establish acceleration of this coupled experiment.

Two seeded founder populations begin with alternating 75/15/10 and 15/75/10 enzyme allocations. Every enzyme column sums to one; increased specialization removes capacity from the other substrate. Efficiency has a throughput cost. Chemical uptake supplies both structural storage and bounded chemical work. The film consumes the contributors' own storage and work, decays, and slows local transport. The founders' initial specialization and cooperation are imposed, so persistence alone is not evidence they evolved.

The default four arms share seeds, positions, controller genotype, initial resources, total source rates, and finite source budgets:

| Arm | Intervention | Question |
| --- | --- | --- |
| niches | complementary half-domain sources | Do local resource regimes associate with different allocations? |
| uniform | spatially uniform sources of equal total mass | Does the niche contrast require spatial segregation? |
| withhold | maintenance rate zero | What happens when costly film maintenance is removed? |
| no_mutation | controller and metabolic mutation zero | Do changes require novel heritable variants? |

The mutation control still inherits founder traits. It is not a no-inheritance intervention. The withholding control removes the entire population's maintenance; it is not a mixed cooperator/free-rider invasion test. An optional `no_crossfeed` arm resets export fractions to zero at every step, preventing exported byproduct even under mutation. Byproduct export transfers fuel mass and remaining chemical potential, reducing the exporter's immediately available work. Cumulative byproduct uptake demonstrates recycling activity rather than a fitness advantage. Film exposure measures access to a physical public effect, not a demonstrated fitness benefit. Patch occupancy is sampled residence rather than individual dwell time; positive phase occupancy is not tracked compartment identity.

Implementation fingerprints come from `Interior.arena_sources()`, identifying the sources loaded into the running module. Editing shared files does not replace physics inside an already running Julia process. Partial checkpoints require matching loaded sources; an explicit `METABOLIC_ALLOW_ARENA_METADATA_MIGRATION=1` exception permits only an `arena.jl` metadata difference reviewed to leave the physical process unchanged, preserving the prior hashes in the result.

Samples retain allocation, efficiency, and cooperation means and variances, controller variance, founder effective count, births and deaths, cumulative uptake by each substrate, patch population and allocation, phase occupancy, film exposure, cumulative maintenance material/work, and incoming/stored/dissipated work plus controller work components. These distinguish trait variation from founder survival and expose resource and maintenance costs. A compact summary retains every sample; individual JSON files make completed comparisons restartable. Atomic serialized checkpoints every sample preserve partially completed comparisons and reject changed implementation hashes. Existing result files are skipped rather than overwritten; use a fresh output directory for changed settings.

Run a bounded pilot:

```sh
METABOLIC_STEPS=400 julia --project=. scripts/metabolic_experiments.jl outputs/metabolic-pilot-20261014
```

Run a longer, explicit comparison after reviewing the pilot:

```sh
METABOLIC_STEPS=12000 METABOLIC_SEEDS=11,29,47 METABOLIC_CONDITIONS=niches,uniform,withhold,no_mutation,no_crossfeed METABOLIC_UPTAKE=1.4 METABOLIC_LIFESPAN=160 julia --project=. scripts/metabolic_experiments.jl outputs/metabolic-long-20261014 > outputs/metabolic-long-20261014.log 2>&1 &
```

The source supply is finite (2,000 units per substrate). Long runs can exhaust it and shift into starvation; compare actual supplied mass and late birth/death rates, not merely terminal population. Lifespan is 80 simulation units; 400 steps of 0.2 reach only one founder lifespan. Small finite pilots establish that mechanisms execute and produce measurable trajectories. They cannot establish stable adaptive diversification, selection for cooperation, cross-feeding, open-ended evolution, or GPU scalability. Longer multi-seed trajectories, mixed-strategy invasion experiments, and resource-budget sensitivity remain necessary before making those claims.

## Measured finite pilot, 14 October 2026

The frozen pilot ran 400 steps in each of five arms for seeds 11 and 29 (24×24×1 grid; 24 founders; capacity 180). Every retained sample was finite; all ten runs recorded unchanged physics source hashes. The largest absolute final material discrepancy was 0.000167 units. The initial tracked study settings and source hashes are in `docs/metabolic-experiment-results.json`.

| Arm | Final population, seed 11 / 29 | Final cumulative byproduct uptake, seed 11 / 29 |
| --- | --- | --- |
| niches | 9 / 0 | 4.70 / 2.88 |
| uniform | 10 / 2 | 3.73 / 3.31 |
| withhold | 25 / 12 | 9.60 / 5.88 |
| no_mutation | 22 / 0 | 6.91 / 3.54 |
| no_crossfeed | 25 / 8 | 0 / 0 |

This parameter regime lost many agents and extinguished two runs by one founder lifespan. Withholding and disabling export left larger terminal populations in both seeds, consistent with immediate costs of film and exported chemical potential. These two-seed terminal comparisons do not establish a selective effect, a persistent population regime, or a benefit to cooperation. The byproduct recycling mechanism executed; adaptive cross-feeding has not been established. A longer, separately labeled resource-affordability study increases uptake capacity to 1.4 and lifespan to 160 while retaining matching within that study; its results must not be pooled with this pilot as one unchanged regime.

The affordability regime's 400-step, three-seed niche check retained populations 39, 24, and 41 (seeds 11, 29, 47) and late births exceeded deaths. A longer checkpointed run was launched with a 600-second wall budget; partial progress is in `outputs/metabolic-long-affordability-20261014.log`. Seed 11 first grew above 100 agents beyond one founder lifespan, then fell to one by time 584 after finite external supply was exhausted. A harness reduction at zero population interrupted sampling; its checkpoint was preserved and resumed with the corrected empty-population diagnostic in an isolated physics snapshot of commit `274ce2e`. The controller had accepted 6,660.38 work units and dissipated/discharged 6,660.13 by the last pre-extinction sample. Extinction did not consume every remaining fuel species: the longer uniform arm retained approximately 274 substrate-B and 267 byproduct units. Accessibility, evolved allocation, controller uptake decisions, and continuing work costs therefore require diagnosis alongside finite supply. In that arm, cumulative remodeling work (3,904) exceeded movement (1,441), signaling (1,043), reproduction (402), and film maintenance (135). These trajectories do not establish a self-sustaining evolutionary regime.

The resumed longer study completed all 15 arms (three seeds × five arms), each for 12,000 steps / 2,400 simulation units, within its wall budget. All arms reached the population capacity of 180 and eventually went extinct. Every sampled state remained finite. The largest absolute material-ledger discrepancy across these much longer trajectories was 0.0515 units (approximately 0.0011% of total supplied plus initial chemical mass); the largest work-ledger discrepancy was 0.00716 units. Float32 transport and repeated exchanges accumulate small ledger drift, so these values are measured tolerances rather than exact conservation claims.

| Arm | First sampled extinction time, seeds 11 / 29 / 47 |
| --- | --- |
| niches | 592 / 568 / 840 |
| uniform | 1048 / 464 / 432 |
| withhold | 400 / 520 / 416 |
| no_mutation | 384 / 472 / 440 |
| no_crossfeed | 504 / 448 / 504 |

Sampling is every eight simulation units, so these are upper bounds within one sample interval. The same finite budget can support multi-generation expansion and eventual collapse; terminal extinction alone would hide this trajectory. Across three seeds, changing source geometry or removing film has variable effects. These results support mechanism execution and motivate supply/work/allocation studies; they do not establish evolution of adaptive cooperation or long-term diversity. The compact tracked JSON now includes all longer summaries, extinction timings, ledger bounds, and loaded implementation hashes; full longitudinal samples remain in the output directory.
