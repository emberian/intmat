# Developmental regimes on the repaired source

Stronger delayed feedback makes repeated developmental switching widespread and
persistent in both completed seeds on the shared repaired source. Variation also
reaches controller actions and metabolic allocation. This changes the prevalence
of modes the baseline can already support; it does not establish that stronger
inhibition is necessary for switching, or that switching always improves ecology.

The matched regimes change adaptation strength from `.3` to `1.5` **and** its
timescale from `6` to `8`. Configuration checks find no other differences within
each seed, and loaded physics hashes match exactly. All runs start from fresh
founders after the finite-precision transfer repair. They are separate from the
[original-source individual examples](constructed-food-development.md), whose
identities and figure remain unchanged.

The table summarizes time ≥1000, retaining identities with at least 50 recorded
observations. Standard deviation is computed within identity and channel, then
combined across channels; the table shows its median across qualifying IDs.

| Seed / regime | Population at 3000 | Qualifying IDs with hidden sign changes | Median hidden std | Median action std |
| --- | ---: | ---: | ---: | ---: |
| 211 / baseline | 0 | 2 / 71 | .00146 | .0293 |
| 211 / stronger feedback | 138 | 105 / 117 | .788 | .223 |
| 307 / baseline | 176 | 1 / 91 | .00179 | .0254 |
| 307 / stronger feedback | 123 | 123 / 126 | .761 | .208 |

These dense traces sample up to the 32 oldest living IDs, not random cells or a
fixed cohort. Qualifying IDs are selected survivors rather than independent
evolutionary replicates. Seed 211's baseline extinction also limits its
observation exposure. Within-ID variation can include environmental forcing and
division partition, not only autonomous developmental activity.

The difference persists late. Between times 1600 and 3000, the stronger seed-307
trace records 39,078 hidden-channel sign changes across 44,742 consecutive ID
pairs; its baseline records none across 44,789 pairs. Full-population snapshots
every 25 time units independently capture 24,488 channel changes across 8,793
same-ID pairs in the stronger regime, versus none across 10,167 baseline pairs.
Of the stronger changes, 19,081 occur outside the oldest-32 sample at both
observations. Seed 211's stronger run likewise varies in both sampling schemes.
The coarse cadence can alias fast cycles; these counts establish recorded
changes, not precise periods or autonomous oscillation.

Controller variation is not confined to hidden state. Median within-ID standard
deviation of derived metabolic shares changes from `.0135` to `.0849` in seed 211
and `.00976` to `.0824` in seed 307. Those shares are regulatory outputs, not
realized material or work fluxes. Their propagation into sustained construction
and spatial organization requires the ecological measurements in the
[dependence study](constructed-food-dependence.md).

The endpoint outcomes resist a universal benefit claim: stronger feedback retains
seed 211, whose baseline becomes extinct, while seed 307 ends with fewer cells
under stronger feedback. Founder ancestry also changes without a single richness
ranking: seed 307 retains eight founder lineages with effective count 1.95 in the
baseline, versus five with effective count 2.31 under stronger feedback. These are
ancestry statistics, not genotype diversity. Two seeds and this combined regime
intervention support a descriptive dynamical result, not an adaptation assay or
isolated causal attribution to gain alone.

[Compact results, actual configurations and source hashes](constructed-food-fixed-development.json)
retain denominators and both observation schemes. Full traces remain in the
two `outputs/constructed-food*-fixed-20261015` roots; analysis helpers and detailed
readouts are in `outputs/constructed-food-fixed-analysis-20261015`.
