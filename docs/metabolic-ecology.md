# Metabolic ecology

`MetabolicEcology` composes a `ResourcePool`, the existing `NeuralEcology`
controller and an optional `BoundaryMaintenance`. It preserves the old neural
state format. Agents have inherited enzyme allocation simplexes, assimilation,
chemical efficiency, cooperative investment and partial oxidation export. All scalar traits stay in
[0,1]; mutation renormalizes the enzyme budget. Birth inherits the parent's
traits and neural genotype, with seeded mutation. Hidden state and plastic marks
reset at birth. Contribution histories belong to living individual IDs, not descendants; death removes individual entries while the boundary retains cumulative totals.

An enzyme share invested in one substrate cannot be used on another. Efficiency
also reduces catalytic throughput by `1 - 0.7 * efficiency`, so higher energy
yield costs uptake capacity. Assimilation directs matter toward structural cargo
and away from fuel work. Neural uptake actions and substrate availability gate
the shared capacity. Species compete proportionally for pre-step stocks.

The first two neural local sensors and gradients represent species one and two;
remaining species add an allocation-weighted cue to channel one. These replace
the old nutrient and deposit cues. Phase, signal, permeability and social cues
retain their original slots. Outputs are computed from this pre-exchange snapshot
and passed unchanged to the neural lifecycle. Chemical gains also contribute to
lifetime plasticity reward without adding energy twice.

Withdrawn chemical potential is split into useful work and chemical dissipation.
Consumed matter becomes inert structural storage, inert waste or an exported lower-potential third species. A heritable export fraction diverts nonassimilated donor material to that species. Export only occurs from higher-potential species one and two; the exported potential is subtracted before computing useful work and chemical heat. Neighbors can consume exported fuel on subsequent steps. Neither dead biomass
in World nutrient nor eroded deposit supplies fresh fuel. Export cannot increase total chemical potential, and species three cannot export itself.
Boundary investment uses the construction output and inherited cooperative
fraction, pays actual material and work, and competes with movement, signaling,
construction and reproduction. Shared film changes permeability for neighbors;
selection for cooperation is an experimental outcome, not an assumed result.

Mass includes pool fields, World structural fields, agents, waste and shared film.
Chemical potential removed from the pool equals useful work added, retained byproduct potential and chemical
heat. The ordinary controller energy ledger then tracks work and birth heat.
These tests establish implementation balances and trait tradeoffs; they do not
establish increasing complexity or evolutionary superiority.
