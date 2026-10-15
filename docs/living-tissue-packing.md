# Space to return

`PackingTissue` is a separate constitutive option composed around the existing
turnover state and reserve buffers. It changes a material assumption: biomass B
and matrix M share a solid capacity C, so B+M≤C. Precursor N, reserve R and spent
material W are dissolved species and do not consume that solid occupancy.
This is a mixture convention, not incompressible mechanics or a stoichiometric
chemical model. The original turnover law and native 128³ art specimen remain
available under their exact source pins.

The option keeps the original seven tissue fields and finite reserve. Persistent
storage stays 64B/site; it adds no species, per-agent objects, contacts or face
arrays. Base neural/metabolic evolution and precursor diffusion are reused.
Base growth and passive B diffusion are disabled because the new phase stage
owns all growth and biomass transport. The same native exporter, observations,
budgets and recording interfaces serve both laws through host callbacks.

## Funded material changes

The phase stage first renews local precursor with bounded R→N transfer. Its
requested operating work is

```
dt * phase_work * (B+M) * (basal_activity + positive + negative)
```

Only available N pays N→W work; actual paid/requested fraction gates actions.
Basal activity makes quiescent thin tissue capable of regrowth. This charge is
an operating-cost convention and can be paid even when space blocks an attempted
action; it is not calibrated work per unit displacement. Original base metabolic
work is charged in its preceding stage, explicitly rather than folded into this
new convention.

Hopping moves existing biomass first, using old donor states. Optional paid basal
transport can move quiescent tissue; optional slow basal remodeling can dissolve
matrix after local living biomass and its neural support have disappeared. Both
default to zero, retaining the first packing recipe's behavior. These kinetic
rates share the explicitly charged basal operating convention; they are not
claims of calibrated energy per unit displaced or dissolved.

Funded growth then
converts N→B, limited by precursor, the local B template and remaining solid
space. Positive phases convert B→M; negative phases dissolve M→N and withdraw
B→N. Construction preserves solid occupancy; the negative transfers open space.
Growth and matrix remodeling have separate rate parameters. A newly arriving
parcel cannot fund another outgoing hop in the same stage. Later local growth
or conversion follows the declared sequential order.

N/B/M/R each retain declared potential 1; W has potential 0. All material transfers
conserve N+B+M+W+R in exact arithmetic. Only paid N→W changes that declared
potential. Latent signals are dimensionless controls; their integration and
diffusion are not an additional chemical mass or a thermodynamic accounting.
There is no learning, inheritance or selection mechanism.

## Competing inflows and capacity

Write ρᵢ=Bᵢ+Mᵢ, free fraction fᵢ=max(0,1−ρᵢ/C), and donor activation
αᵢ=paidᵢ·(basal_transport+negativeᵢ)/(1+basal_transport).
With bounded signals and nonnegative basal transport, 0≤αᵢ≤1. For neighboring sites,
Pᵢⱼ is the reused symmetric matrix permeability in [0,1], and
tᵢⱼ=tanh(χ(Nⱼ−Nᵢ)). The flux into i is

```
D Pᵢⱼ / (2h²) * [αⱼ Bⱼ fᵢ (1−tᵢⱼ) − αᵢ Bᵢ fⱼ (1+tᵢⱼ)]
```

It is antisymmetric under i↔j, hence conserves biomass. If all sites initially
obey 0≤B and B+M≤C, total outgoing biomass is bounded by
`2D·dt·sum(h_axis^-2)·Bᵢ`. All competing incoming faces together are bounded by
the same factor times `C−ρᵢ`, because each neighboring Bⱼ≤C. A factor≤1 therefore
preserves both donor positivity and recipient capacity; reasoning about one
incoming face alone would not suffice. The implementation's stricter stability
limit includes this summed bound and inhibitory diffusion. All faces read old
buffers, including their paid fractions and recipient free space.

The subsequent growth transfer is explicitly bounded by `max(0,C−ρ)` and N.
Construction preserves ρ; dissolution and withdrawal reduce it. The preceding
base stage has no B exchange or growth and only returns B to precursor, so it
also preserves the occupancy precondition. Native statistics report maximum
occupancy and excess; capture validation allows 16 Float32 eps times C for arithmetic
roundoff, without clipping transport, deleting excess material or masking a
negative inventory.

## Runnable checks and studies

`studies/living_tissue_packing.jl` implements the option; the small adapter
`studies/living_tissue_packing_run.jl` uses the generic turnover exporter. The
parameter recipe `experiments/art/living-tissue-packing64.toml` requests a nearly
uniform 64² sheet, extent 24², dt .005 and t0–80 with finite reserve concentration 12.
Its stronger, faster inhibitory communication is a separate recipe from the
original native art specimen. An interesting result must show actual opening,
filling or connectivity changes, not merely reduced high-frequency power.

`test/living_tissue_packing.jl` passes 50 focused checks: donor/free-space growth,
zero-fuel suppression, all six competing neighbor inflows, fully occupied
blocking, repeated exchange, no re-spending arrivals, anisotropic spacing,
material/potential balance, child scheduling, recording, exact restart and paid
basal transport/remodeling, including a large-baseline capacity challenge.
`test/living_tissue_exporters.jl` passes 26 composition checks for both original
turnover defaults and packing. Strict packing checkpoints include this new law's
hash and the reused turnover/base/shared hashes; they cannot silently restore
as the original law. No accelerator result is claimed for packing yet.
