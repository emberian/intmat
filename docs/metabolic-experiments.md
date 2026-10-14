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

Samples retain allocation, efficiency, and cooperation means and variances, controller variance, founder effective count, births and deaths, cumulative uptake by each substrate, patch population and allocation, phase occupancy, film exposure, and cumulative maintenance material/work. These distinguish trait variation from founder survival and expose resource and maintenance costs. A compact summary retains every sample; individual JSON files make completed comparisons restartable. Atomic serialized checkpoints every sample preserve partially completed comparisons and reject changed implementation hashes. Existing result files are skipped rather than overwritten; use a fresh output directory for changed settings.

Run a bounded pilot:

```sh
METABOLIC_STEPS=400 julia --project=. scripts/metabolic_experiments.jl outputs/metabolic-pilot-20261014
```

Run a longer, explicit comparison after reviewing the pilot:

```sh
METABOLIC_STEPS=12000 METABOLIC_SEEDS=11,29,47 julia --project=. scripts/metabolic_experiments.jl outputs/metabolic-long-20261014 > outputs/metabolic-long-20261014.log 2>&1 &
```

The source supply is finite (2,000 units per substrate). Long runs can exhaust it and shift into starvation; compare actual supplied mass and late birth/death rates, not merely terminal population. Lifespan is 80 simulation units; 400 steps of 0.2 reach only one founder lifespan. Small finite pilots establish that mechanisms execute and produce measurable trajectories. They cannot establish stable adaptive diversification, selection for cooperation, cross-feeding, open-ended evolution, or GPU scalability. Longer multi-seed trajectories, mixed-strategy invasion experiments, and resource-budget sensitivity remain necessary before making those claims.
