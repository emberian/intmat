# Where the garden's feedback can act

An analytical map of the native recipe, using its actual world-total source
rates and physical extent. These are operating-point estimates, not observed
relaxation times or conditions sufficient for rich behavior. The executable
[Julia study](../studies/garden_operating_map.jl) reads the recipe and shared
cellular defaults; [JSON](garden-operating-map.json) records assumptions,
source hashes and nine candidate contrasts.

| Clock | Model time |
| --- | ---: |
| Neural response | 0.8 |
| Nominal minimum travel over exchange radius | 3 |
| Bare diffusion over exchange radius, `L²/D` | 5.4 |
| Adaptation | 8 |
| Free desorption | 8.33 |
| Mark relaxation, frozen hidden activity | 45.45–50 |
| Catalyst decay | 125 |
| Initial raw concentration / mean input rate | 128 |
| Finite raw supply horizon | 1500 |

These scales overlap enough to make delayed feedback plausible. The question is
which paths have appreciable gain and remain accessible to a cell as it moves.
Residence time can be shorter than catalyst memory; retained traces can outlive
their authors. A long-lived trace is useful memory only if a later response
depends on it. The [matched interventions](native-memory-interventions.md)
measure that dependence; the [research shelf](research/learning-materials.md)
connects it to causal attribution and driven information-processing capacity.

## Catalysis and remodeling need not share an operating point

Each raw source supplies one mass unit per time across volume64, hence mean
concentration rate `1/64`, not1. With initial raw concentration2 and catalytic
rate6, matching that input requires estimated enzyme concentration0.002865 per
species. At this estimate, the total enzyme lowers permeability only to0.897.
Each binding capacity is0.01719 versus initial free-food concentration0.3.
Efficient food production can therefore occur before strong transport remodeling.
Actual local enzymes may accumulate above this estimate, and local raw depletion
changes it; the calculation does not establish weak feedback in the recorded run.

Three analytically chosen contrasts, with all other recipe values held fixed:

| Catalytic rate | Transport barrier | Permeability at supply match | Diffusion / adaptation |
| ---: | ---: | ---: | ---: |
| 6 | 20 | 0.897 | 0.752 |
| 0.6 | 20 | 0.466 | 1.448 |
| 0.6 | 100 | 0.149 | 4.542 |

Slower catalysis demands more paid construction for the same release; a stronger
barrier makes that construction affect transport sooner. Either can prevent
establishment or trap resources. Those are hypotheses to explore in small
matched worlds, alongside unaltered controls. They are not recommended production
settings on the evidence so far.

The map explicitly distinguishes world size from spatial resolution: refining
the grid alone changes none of these physical clocks. Domain size changes the
long-wavelength modes, and source totals must scale with volume to preserve local
input. Changing a timestep changes the discretization, not the intended physical
timescale. The sixteen checks cover those distinctions and null-source/transport
limits.

```sh
julia --project=. bin/interior.jl study garden_operating_map \
  experiments/native/constructed-garden-64.toml outputs/my-operating-map.json
julia --project=. test/garden_operating_map.jl
```

Use a fresh output path. Transport uses the characteristic convention `L²/D`,
not a three-dimensional first-passage time. The nine candidate points vary
catalysis and barrier independently; nonlinear uptake, evolving expression,
finite funding and spatial heterogeneity remain for actual simulation to resolve.
