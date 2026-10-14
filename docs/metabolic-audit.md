# Coupled metabolic audit

The audit in `test/metabolic_composition.jl` checks the actual combined processes,
including finite external source/sink budgets, diffusion, funded film maintenance,
film decay, births and deaths. It uses Float64 sums of Float32 reservoirs with an
absolute tolerance of 5e-5 on the small fixture. These are numerical model checks,
not calibrated physical units or evidence of stable cooperation.

The material inventory is world nutrient and mineral deposit, organism structural
storage, cumulative inert waste, chemical field mass, and shared boundary film.
Resource source additions and sink removals correct its initial baseline. A
resource shared by the organism wrapper and scheduler is one reservoir.

The combined energy inventory is chemical potential remaining in the resources,
organism stored work energy, controller heat, and metabolic chemical heat. Source
and sink potential correct the initial baseline. `BoundaryMaintenance.spent_work`
is an audit of maintenance contributions already charged to controller heat; it
must not be added as a second energy reservoir. Inert material returned by death
or film decay does not replenish chemical fuel.

The inheritance fixture assigns contrasting substrate allocations and metabolic
traits, permutes organism order, forces division, and then removes a parent.
With zero trait mutation, each child's traits match its parent; lifetime recurrent
state resets. A separate downhill conversion chain accounts for useful energy and
heat across intermediate byproduct consumption and rejects conversion back uphill.

Recording is tested at the actual `record_state` export boundary. Checkpoints
must preserve shared resource, medium, phase, boundary and world-deposit aliases;
restored continuation must match the uninterrupted trajectory. Process-hook
references alone are not immutable captures.

These checks establish bounded bookkeeping and deterministic composition. They
do not show cooperation resists invasion, that boundaries form persistent
compartments, or that scale changes preserve behavior. Those claims require
separate competing-trait interventions, spatial diagnostics, and refinement.

The scheduler validates ownership before changing state: a boundary wrapper owns
its internally advanced medium, and a metabolic wrapper owns its internally
advanced controller. Scheduling either owned stage separately is rejected, as is
repeating the same mutable process instance. Borrowed resource and boundary
references remain valid separately scheduled stages. Repeated immutable stateless
operators remain permitted.
