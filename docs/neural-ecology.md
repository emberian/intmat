# Neural ecology in a continuous arena

`NeuralEcology` adds a population of resource-limited recurrent agents to the
continuous fields. It is an explicitly designed artificial-life substrate; a
richer controller and ecological feedback do not establish that open-ended
organisms, intelligence, or evolutionary novelty have emerged.

```julia
world = World(Grid(32,32; extent=(32,32,1)))
phase = ImplicitPhaseSeparation(world)
medium = ActiveMedium(world, phase)
population = NeuralEcology(Float32[4 12; 4 12; 0 0]; medium)
step!(world, (population, medium, phase), 0.01)
```

Agents occupy floating-point periodic **cell coordinates**, with three rows even
in 2D. The default controller has eight recurrent hidden units. Inherited
`weights` are separate from lifetime `plastic_marks`: expression uses their sum,
reward-modulated Hebbian plasticity changes only the marks, and offspring start
with zero marks and zero hidden state. Children inherit the parent's genotype
plus seeded mutation, lineage and generation metadata. Birth offsets are seeded
as well. Sorting by persistent IDs before observation and exchange makes input
array reordering irrelevant to the trajectory and mutation sequence.

Each step snapshots senses and controller outputs before agents exchange resources.
Agents sharing a cell receive scarce nutrient in proportion to their requested
uptake. Contested membrane erosion uses the same proportional allocation, with
its demand determined by each agent's paid action budget. Float32 shares round
down before the resource grid is debited, leaving tiny roundoff residuals rather
than favoring a selected ID. Allocation uses pre-stage stock: material returned
at death and freshly deposited membrane become available to other agents in the
next step. Persistent-ID ordering is retained for reproducible birth events, not
for feeding priority. The 32 sensor inputs
follow hidden state in the controller input vector:

| Input indices (after hidden units) | Meaning |
| --- | --- |
| 1–6 | nutrient, membrane deposit, signal, phase composition, interface, permeability |
| 7–9 | stored material, abstract energy, age/lifespan |
| 10–27 | central xyz differences for each of the six local channels, in channel order |
| 28 | nearby agent count, squashed as count/(1+count) |
| 29–31 | mean periodic neighbor displacement xyz, divided by neighborhood radius |
| 32 | constant bias |

Resource and spatial inputs use `x/(1+abs(x))` except normalized age and neighbor
displacement. Periodic occupied-cell bins restrict neighbor lookup to a local
cube; distance filtering gives a spherical neighborhood. Wrapped cells are
visited once, including in tiny domains. z gradients, neighbors, movement and
birth offsets operate in 3D, and vanish on a one-cell z axis.

Nine outputs follow hidden targets: xyz movement, uptake, deposition, erosion,
signal emission, division, signed chemical remodeling. Hidden state relaxes toward
a tanh target with `response_time`. Founder feeding/division and food-gradient
movement biases are evolvable initial conditions, not fixed policies. The other
outputs and recurrent weights start with reproducible random variation.

Uptake removes available nutrient and splits it between catabolism and structural
assimilation. `catabolic_fraction` (default 0.5) becomes **inert waste** in the
material ledger; the remainder becomes agent storage. Only catabolized fuel
credits abstract chemical work, at `energy_yield` (default 2) per unit fuel.
`catabolic_fraction=0` supplies no work, and `=1` supplies no new structural matter.
These are global chemistry parameters, separate from the inherited controller.

`energy_input` includes founders' initial work and this harvested chemical work.
Basal metabolism, active work, division overhead, and energy remaining at death
become `heat`. Thus stored energy plus heat equals energy input, up to
floating-point accumulation. Membrane deposition debits storage; erosion credits
it; division splits stored matter between parent and child; death returns it to
the final occupied nutrient cell. Nutrient, membrane, agent storage and inert
waste form the material ledger. The waste reservoir is globally accounted, not a
spatially transported or biologically reusable resource. Signal and phase-bias
fields are information/forcing variables rather than additional matter reservoirs.

This chemistry closes a recycling loophole: returning structural material at death
can recover matter, but every subsequent harvest spends a fuel fraction into
unusable waste. In a closed arena, additional work is bounded by `energy_yield`
times the initial usable material. Sustained turnover must eventually depend on
fresh feedstock; repeated death/reuptake cannot generate unlimited work from the
same matter. Actual external nutrient additions should therefore be recorded in
an arena material-input ledger as well as in the ecology's harvested-work ledger.

Movement, membrane changes, emission and remodeling share the paid `work_rate`
budget in proportion to controller demand. An exhausted budget scales all these
actions down together. A zero `work_rate` disables active actions; it does not
make them free. Division pays its own `birth_energy` overhead, requires at least
`birth_material` storage and twice the overhead in energy, and is suppressed at
capacity. `birth_material` is a storage threshold; it is not destroyed on birth.

An attached medium provides local phase/interface cues and accumulated signals.
Permeability reduces nutrient access and active movement using the minimum of the
start and proposed endpoint cells. Displacement is capped at 0.45 cells per axis
per advance so agents cannot jump several cells across a narrow barrier. This is
porous drag, not a hard collision solver. Passive fluid drift remains possible.
Deposited membrane changes the next medium permeability update, erosion removes
that resistance, and paid signals can open it. Signed remodeling increments drive
the medium's chemical-potential forcing. Advance population before medium, then
phase, for this feedback to operate in the same stage. Observe/refine the medium
after external membrane edits when immediate permeability is needed.

Use `ecology_statistics(population)` for population/turnover, lineage count,
material and energy ledgers, action work, deposited/eroded totals, and plasticity
norm. `last_actions` contains raw controller targets rather than executed work;
it is useful for policy observations alongside the paid-action totals.
`record_process_fields` excludes the medium reference because it is a separate
process sharing world/phase arrays. Checkpoints must preserve or reconnect those
aliases. The constructor retains the old `Ecology` parameter names and adds
`medium`, `plasticity_rate`, `plasticity_decay`, `signal_rate`,
`remodeling_rate`, `neighbor_radius`, `catabolic_fraction`, and `energy_yield`.

Sustained population turnover needs replenished nutrient and a balance between
uptake, local permeability, work costs, lifespan and division overhead. Capacity
is a hard safety bound, not a carrying-capacity prediction. Compare matched-seed
runs with plasticity, mutation, signals and remodeling disabled to test which
feedbacks actually matter; a colorful field or nonzero mutation alone is not
such evidence.
