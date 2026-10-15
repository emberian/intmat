# Designing a material with room to learn

Reading shelf and experimental design, 2026-10-15. The target is a repertoire of
responsive, persistent, composable organizations—not a single winning morphology.
The proposals below are transfers to Interior, not results established by these
papers for our ecology. [Machine-readable sources](learning-materials-sources.json)
record access, versions, local PDF hashes and licenses; [BibTeX](learning-materials.bib)
is available for reuse.

## Three different kinds of learning

The native garden already combines inherited controller weights, within-life
hidden/adaptation state, and bounded activity-dependent plastic marks. The
implemented mark update is Oja-like reinforcement with decay; only hidden-output
rows acquire marks. Action and catalytic readouts are inherited. Thus a cell can
change its response with experience, but the present rule does not explicitly
assign credit for an ecological outcome. Selection across descendants is a
separate process. Better future reward does not follow merely from changing marks.
See [kernel](../../src/device_cellular_kernels.jl),
[parameters](../../src/cellular_initialization.jl), and the
[existing feedback analysis](../cellular-feedback-analysis.md).

A third learning process can live outside the ecology: our experiment designer
builds a map of what laws and parameterizations can do. Keep those three claims
separate: **development in an individual; heritable change across generations;
our discovery of capable worlds**. A good outer search can discover an oscillator
without the simulated organism having learned to oscillate. Conversely, failure
of a particular outer search says little about what a law could support.

## Nine useful primary works

| Work | What it contributes | Transfer to Interior |
|---|---|---|
| [Arora et al., MAttr (2026)](https://arxiv.org/abs/2609.25518v1) | Learns nested masks across component budgets, optimizing behavior under interventions. Section 3 defines soft interchange interventions; section 4 evaluates hard masks. | Order mechanisms by the response they causally sustain. Start with cheap explicit block ablations; our small controller does not require an attribution neural network. Learned rankings depend on chosen behavior, counterfactual source, component basis and horizon. |
| [Dambre et al., Information processing capacity (2012)](https://doi.org/10.1038/srep00514) | Measures reconstructability of orthogonal functions of input history with linear readouts, including recurrent and reaction–diffusion examples. Discusses tradeoffs between nonlinear transformation and memory. | Probe delayed nutrient signals and products of signals using hidden state, adaptation, marks and local chemistry separately. Held-out readout accuracy is a capability assay; it is not evidence that the cells use that information. |
| [Natschläger, Bertschinger & Legenstein (2004)](https://proceedings.neurips.cc/paper/2671-at-the-edge-of-chaos-real-time-computations-and-self-organized-criticality-in-recurrent-neural-networks.pdf) | Relates driven threshold-network dynamics and time-series computation; separates input-induced differences from persistent initial-state differences. | Run paired copies with identical input but perturbed state, and paired inputs from identical state. Measure both forgetting and discrimination. The useful regime need not maximize either divergence or memory alone. |
| [Cohen et al., Edge of Stability (2021)](https://arxiv.org/abs/2103.00065) | Full-batch gradient descent often operates near a loss-Hessian eigenvalue of `2/learning_rate`. | Distinguish optimization stability from ecological/recurrent criticality. Our Hebbian-like marks and evolved weights are not that training process; this threshold is not a recipe for choosing the simulator timestep. |
| [Stricker et al., synthetic gene oscillator (2008)](https://doi.org/10.1038/nature07389) | Constructs a persistent tunable oscillator using a feedback architecture informed by modeling and experimental calibration. | Design and test positive drive plus delayed inhibition as a controllable motif; then couple motifs through habitat and bonds. Use periods, phase response and recovery to characterize it rather than asking whether pictures look alive. |
| [Nisioti & Moulin-Frier, niche construction (2023)](https://arxiv.org/abs/2305.09369v1) | Studies interacting niche construction, adaptability and environmental diversity in an explicit evolutionary model. | Treat catalyst and retained-food fields as candidate ecological inheritance. Compare living cells in their own material history with matched cells in exchanged histories. Their abstract environmental/tolerance-curve model is not our physical chemistry. |
| [Cully & Demiris, Quality-Diversity (2017)](https://arxiv.org/abs/1708.09251) | Organizes search around collections of diverse high-performing solutions rather than one global optimum. | Archive lawful parameterizations across memory, motion, transport remodeling, recovery and topology. Avoid survival alone as the ranking function; maintain several descriptor views and retain surprising failures. |
| [Etcheverry, Oudeyer & Reinke, HOLMES (2020)](https://arxiv.org/abs/2005.06369v1) | Grows a hierarchy of representation modules for exploration, freezing and splitting modules to retain old distinctions while discovering finer ones. Tested on self-organized patterns. | Our discovery archive can gain descriptors as new behavior appears. Freeze descriptor versions and keep raw trajectories so an early image embedding does not decide permanently what counts as novelty. |
| [Mordvintsev et al., Growing Neural Cellular Automata (2020)](https://distill.pub/2020/growing-ca/) | Shows growth, persistence and repair under distinct training protocols, including a pool of prior states and damage during training. | A useful separate engineered-material substrate: optimize persistent/recoverable local rules, then compose with accountable chemistry and optics. Target-trained morphogenesis and autonomous ecological evolution are different achievements. |

MAttr's identity and method were verified in its primary full text. The quoted
persona-unlearning paper [also exists](https://arxiv.org/abs/2609.39882v3), but the
claimed 214-example connection was not verified and is not used here. Persona
editing is not our proposed objective. Scry found the MAttr authors' public
[announcement](https://x.com/aryaman2020/status/2102800933659000846) and
[method explanation](https://x.com/aryaman2020/status/2102800937769427281);
these led back to the paper rather than substituting for it.

## A clock map for the actual native recipe

These are estimates from [constructed-garden-64.toml](../../experiments/native/constructed-garden-64.toml)
and defaults, not measured relaxation times of the whole garden. In particular,
transport geometry, uptake, saturation, paid motion and changing contacts alter them.

| Process | Nominal scale in model time | What to measure |
|---|---:|---|
| Hidden response | `.8` | local impulse gain and decay/overshoot |
| Adaptation | `8` | phase lag and mode switching |
| Plastic marks | `1/(.02+.002*h²) ≈ 45–50` | persistence after input removal; erase/reinstate response |
| Free desorption | `1/.12 ≈ 8.33` | bound/free response with local capacity held fixed |
| Catalyst lifetime | `1/.008 = 125` | chemical history after constructor departure |
| Bare diffusion over exchange radius | `.18²/.006 = 5.4` | local spreading time; convention here is `L²/D`, not a 3D RMS first-passage time |
| Exchange-radius crossing at nominal speed | `.18/.06 = 3` | actual residence distribution; paid speed can be lower |
| Feedstock supply horizon | `1500/1 = 1500` each | total supplied inventory and location, separately from accessible food |

The five-unit observation cadence of the completed garden undersamples its
fastest controller clock. The response-time and adaptation-time assays need
small high-cadence windows (for example `.05–.1` spacing), while whole-volume
photographs can remain sparse. Record rates and fluxes, not only states.

Explore dimensionless ratios such as `tau_adapt/tau_hidden`,
`tau_marks/tau_adapt`, `tau_catalyst/tau_residence`, and `L²/(D*tau_adapt)`.
Vary a ratio while holding the relevant integrated supply or investment fixed
where feasible. A joint rescaling of all rates mainly changes the clock; changing
one loop can create a different system. Local eigenvalues, impulse curves and
observed trajectories should agree before naming an oscillatory mechanism.

## Experiments that can change our design

**1. A driven controller atlas.** Freeze a bank of actual founder and recovered
controllers; supply identical recorded or synthetic sensory histories. Cross
adaptation strength, mark rate/decay and recurrent gain. For each controller
measure saturation, impulse recovery, fading of initial conditions, delayed-input
readout and nonlinear-input readout. Train readouts on one signal sequence, test
on another, and include a shuffled-time control and a direct-input baseline.
First use prescribed inputs to separate recurrence from ecological forcing;
then test selected regimes in the full coupled material. Include bistability,
excitable transients and stable differentiation as successes alongside oscillation.

**2. Paired causal windows.** Branch new experimental copies from verified states;
never alter historical records. Compare intact dynamics with adaptation clamped,
marks frozen, marks erased, neighbor hidden input exchanged, or catalytic
regulation held fixed. Erasing information is an intervention, not a conservative
physical process: report what was changed. Match same-ID state where possible,
use multiple windows and seeds, and report immediate responses separately from
later survival/reproduction divergence. Also test reinstatement and combinations:
a mechanism can compensate for another, so one-at-a-time null effects do not show
irrelevance. A nested retained-component curve is a practical MAttr-inspired
measurement; calling it MAttr requires implementing and validating that method.

**3. Ecological memory transplantation.** Take cell states and environmental fields
from two matched histories and cross them in separate test worlds. Preserve the
chosen inventory explicitly; do not silently add food during exchange. Separate
geometry, retained food and catalyst history. Quantify the benefit/cost to the
constructor, descendants and nearby nonconstructors through actual uptake/work
fluxes. This tests where usable memory resides and whether it remains accessible,
not merely whether a long-lived field exists.

**4. A repertoire, not a winner.** Use a small initial space with clock ratios and
construction strength as axes. Store every executed recipe, seed, source and
trajectory summary. Archive multiple behaviors: recoverable injury, spatially
localized resource use, phase relationships, persistent differentiated states,
channel formation and propagating fronts. Promote diverse candidates to larger
and longer GPU tests; retain random exploration and a few unpromising controls so
early descriptors do not close the search. Measure complete simulation throughput
and experiments-per-hour separately. Render actual selected trajectories as clear
volumetric inclusions, layered refractive material and emission from measured
activity; the visual archive is another way to notice structure, not its proof.

## Architecture implied by the experiments

Keep physics/accounting, controller state update, heritable parameters,
interventions, observations and outer search as explicit interfaces. The same
Julia law should run a cheap prescribed-input assay and the coupled GPU world.
Batch independent small assays; keep device state resident; copy compact summaries
at the fast observation cadence. A long trajectory is valuable once we know what
its state can express and have instruments capable of seeing it. No requirement
that every candidate be differentiable follows: discrete events can be explored
by finite interventions or derivative-free search, while smooth controller
windows can expose Jacobian-vector products. We can deliberately engineer a
rich substrate without prescribing the exact organization that must inhabit it.
