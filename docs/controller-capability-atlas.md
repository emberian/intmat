# A driven controller capability map

The [interactive atlas](https://emberian.github.io/intmat/laboratory/) compares
nine recurrent-gain/adaptation regimes across eight founder controllers. It uses
the native `dc_controller!` and `dc_regulate!` Julia kernels, with prescribed
inputs and clamped physiology. No ecological, mechanical or lifecycle steps run.
This is a developmental capability assay, not an arena benchmark or demonstration
of useful learning. [Measured data](../site/laboratory/report.json).

## Experiment

Two independently seeded uniform input sequences drive complementary free-food
sensory channels (`.5±.4u`) and a tonic bias. Neighbor and other physiological
inputs are zero. Each sample lasts.4 model time; hidden response time is.8 and
adaptation time8. Eight founders from seed211 are paired with weight-identical
twins whose initial hidden state differs by.001. Hidden plastic marks retain the
existing production rule (`rate=.002`, decay.02); action weights are inherited.

Each sequence has200 warm samples and1200 retained samples. A linear ridge
readout uses hidden and adaptation state (16 features) to predict inputs delayed
1–12 samples or the product of the previous two inputs. Centering and scaling use
training data only; ridge is fixed at.001 with no test tuning. A direct-current-
input quadratic baseline and a circularly time-shifted held-out-state control
retain negative R² scores. No test-set selection or clipping contributes to the
displayed results. The outer readout is our measurement, not a learned cell policy.

| Gain / adaptation | Median lag1 R² | Median lag8 R² | Median lag12 R² |
| --- | ---: | ---: | ---: |
| .6 / 0 | .221 | .051 | .039 |
| 1 / .3 | .231 | .050 | .040 |
| .6 / .9 | .229 | .095 | .077 |
| 1.4 / .9 | .027 | .011 | .018 |

For this input family, stronger slow adaptation with lower recurrent gain retains
more readable distant input than the nominal regime. Strong adaptation combined
with stronger recurrence performs worse on these readouts. More feedback is not
monotonically better. Nonlinear product prediction has negative median R² in all
nine conditions. This particular encoding, amplitude, readout regularization and
founder bank may miss other capabilities; the result does not exclude nonlinear
computation, useful multistability or more suitable sensory histories.

Most same-input twin differences collapse near numerical precision; strong-
adaptation regimes retain small differences. These are finite perturbation
responses, not Lyapunov estimates. Population dynamics, local field memory and
network communication are absent from this assay.

## Discretization and checks

Repeating the full matrix at dt.025 instead of.05, with unchanged physical sample
duration, changes a regime's median score by at most.01094 across the13 targets.
The largest changes occur at adaptation.9. This supports the coarse qualitative
comparison over this window; it does not establish trajectory convergence.
[Score comparison](../site/laboratory/timestep-comparison.json).

Ten checks cover held-out linear recovery, training-constant features, constant
targets, deterministic inputs/trajectories, retained-sample alignment, perturbed
twins and cadence validation. Full coarse matrix execution98.28s excludes package
load and prior test compilation; this is CPU controller-only timing. It must not
be represented as ecological throughput.

Full input/feature/twin traces and executed-driver copies remain under
`outputs/controller-capabilities-20261016` and its `-half` companion. Compact
reports retain source hashes. The browser only selects recorded conditions;
the [live instrument design](exploratory-instrument.md) describes the next bridge.

```sh
julia --project=. bin/interior.jl study controller_capability_atlas outputs/new-atlas
julia --project=. studies/controller_capability_atlas.jl outputs/new-atlas-half --half-step
julia --project=. test/controller_capability_atlas.jl
```

Next: vary sensory encoding and receptor scales; couple selected regimes back
into native material; design an action-learning rule with an explicit local
outcome signal and eligibility horizon. An instrument should discover capable
mechanisms, not just rank whatever happened to be inherited.
