# Feed the patches, then let them persist

2026-10-16. Separate bounded recipe and causal probe using the existing native
Julia law. This strengthens an actual material feedback without prescribing a
shape or assuming that a short quiet scene establishes incapability.

The [candidate recipe](../../experiments/native/living-patch-recruitment-16.toml)
changes only raw diffusion from zero to 0.006 and catalyst decay from 0.2 to 0.05.
All founders, initial inventories, finite source budgets, source-mask schedule,
food diffusion, sorption, motion, contacts and paid learning remain the same as
the moving/dividing scene. The [probe](../../studies/living_patch_recruitment.jl)
crosses those two changes in four independently initialized matched worlds.

## Mechanism in the implemented law

The native catalytic reaction consumes raw material at
`kcat × enzyme × raw/(K×voxel_volume + raw)` and delivers soluble precursor.
Precursor can fund battery, growth and new extracellular catalyst; actual
affordability decides which requests are fulfilled. This closes a possible
positive loop: paid catalyst → local food access → retained cellular resources
and descendants → more paid catalyst. The catalyst itself has no autonomous
reproduction law; this is a cell–material feedback.

Three mechanisms oppose or redirect it. Raw depletion reduces catalytic return;
finite world funding ultimately limits replenishment. Sorption removes food
from immediately importable form and releases it later. Motion can carry a
producer away from its patch. Native founder weights already favor soluble-food
gradients in all three motor directions; descendants can remain near a productive
birthplace, but no retention, adhesion or cooperation is forced. Catalyst changes
chemical permeability, not mechanical drag or anchoring.

With zero raw diffusion, a patch can exhaust the material beneath itself while
remote raw stock accumulates. This is especially visible early; the separate
living authority replay finds local raw abundance again late in the baseline,
while stores and funded anabolism remain low. Raw availability is therefore not
the sole late bottleneck. Keeping an exhausted catalyst alive longer is
not enough. Raw diffusion makes a patch recruit surrounding stock; longer
catalyst life then gives its location time to become a revisitable resource
and a legacy for nearby descendants. Local depletion, competing patches and
finite supply can still limit or dissolve it. This supports a hypothesis of
productive patches and trails; the current equations do not justify naming it
a Turing instability or predicting a selected wavelength.

## Operating scales and constraints

| Quantity | Baseline / candidate | Interpretation |
|---|---|---|
| Catalyst lifetime | 5 / 20 | Material outlasts a two-time credit horizon in both laws. |
| Radius / nominal speed | 0.18 / 0.06 = 3 time | A producer can leave before its catalyst decays; actual motion is paid and contact-limited. |
| Raw `R²/(6D)` | infinite / 0.9 time | Candidate raw stock can replenish a support-sized depletion zone. |
| Raw `sqrt(D×tau_catalyst)` | 0 / 0.346 | A characteristic diffusion length, not a measured capture radius or pattern spacing; RMS displacement is sqrt(6) times this. |
| Mean founder spacing scale | `(8/24)^(1/3) ≈ 0.693` | Candidate recruitment can begin to couple neighboring patches over their lifetime. |
| Food diffusion | 0.006, further reduced locally | Body packing and matrix permeability alter soluble access. |
| Desorption clock | 1/0.12 = 8.33 | Food can remain bound after the producer leaves; binding is not immediate nutrition. |
| Hidden/adaptation clocks | 0.8 / 8 | The native representation can respond on transport and residence timescales. |

At uniform average source density, a stationary radius-0.18 sphere intercepts
only `0.8/8 × (4π×0.18³/3) ≈ 0.00244` raw mass/time. Even ideal conversion and
respiration yield about 0.00635 usable energy/time, below founder basal work
0.0128 plus neural work 0.0016. This is an ongoing-source footprint calculation,
not a bound on a moving resident's total intake: existing food, transport and
past construction outside the footprint can all contribute. It explains why
adding raw transport can matter more than increasing a construction command.

Sorption equilibrium has `Kd=0.12/4=0.03`; at low free concentration the local
retardation factor can be large. Matrix permeability is
`1/(1+20×catalyst_concentration)` and multiplies body permeability. We deliberately
leave these laws unchanged to see whether replenishment makes retention useful,
as the earlier [paid-patch access assay](../paid-patch-access.md) suggested.
That earlier assay had no consuming cells; the present probe includes them.

## Eight-time causal prefix

One founder bank, common exploration seed, dt 0.02; all four received the same
6.40000037 supplied mass. Maximum conservation residual was below 1.6e−13.
All starts had 24 residents; no deaths occurred in this prefix.

| Raw diffusion | Decay | Living / births | Converted precursor | Funded catalyst | Retained catalyst | Catalyst near living |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.2 | 27 / 3 | 1.1887 | 0.4176 | 0.1047 | 46.5% |
| 0 | 0.05 | 27 / 3 | 1.2946 | 0.4182 | 0.2921 | 39.6% |
| 0.006 | 0.2 | 37 / 13 | 4.4991 | 0.5558 | 0.1783 | 73.3% |
| 0.006 | 0.05 | 37 / 13 | 5.2520 | 0.5587 | 0.4102 | 65.3% |

Raw recruitment increased realized conversion and funded construction in living
worlds. Slower decay alone mainly retained catalyst; its additional conversion
was much larger when raw material could move. The candidate also retained 1.1594
bound food versus baseline 0.2019, so some new production is stored rather than
immediately used. Catalyst-weighted raw saturation rose from baseline 0.150 to
0.262 in the candidate (0.491 with raw transport and fast decay). Extra catalyst
can therefore make local substrate scarcer per enzyme even while total output
rises. No single mean concentration describes that feedback.

At time 8 the candidate's maximum material receptor magnitudes were 0.218,
0.642 and 0.784 for raw balance, raw total and catalyst total. These do not show
hard saturation of all receptors; actual per-cell sensitivity and learned action
offsets still matter. Candidate funded/requested catalyst was only about 8.6%,
so finite precursor and competing native anabolic requests remain active
constraints. This recipe is not an unlimited-construction workaround.

“Near living” is the fraction of catalyst mass inside the union of geometric
exchange-radius spheres at the recorded positions; it is not producer identity,
actual uptake attribution or a claim that every nearby patch is profitable.
Effective catalyst volume is also recorded as `(sum(e))²/sum(e²)` divided by
voxel count; it describes localization without a thresholded picture.

## How learning becomes a spatial lever

The existing paid joint policy controls reserve, investment and A/B allocation,
not motion. By reconnecting constructed patches to remote stock, this recipe
changes the consequences available to that policy. A useful next comparison
holds the recipe fixed and contrasts paid correct credit with matched noise/no
actor credit; it should ask whether residents invest selectively at productive
patches and reduce investment after depletion. Do not attribute the four-arm
physical improvement to learning, which is shared across all arms.

If catalyst lifetimes of 20 create delayed benefits beyond the current horizon,
measure them with the existing decision-aligned lag probe. Do not lengthen the
credit queue from the decay constant alone: the producer may leave, and benefits
may go to descendants or neighbors. The paired local trend sensor is a candidate
for distinguishing a replenishing niche from an exhausted one at the same
instantaneous raw reading; its gain-controlled causal test remains separate.

## Completed 32-time comparison

The [compact report](../figures/living-patch-recruitment.json) contains every
two-time sample and executed study hashes. All four received 25.59999999 supplied
mass, exhausting the same declared budgets. The largest conservation residual
across these records is below 5.6e−13. All runs used one founder bank and common
exploration seed, so this is a matched causal comparison, not broad replication.

| Raw diffusion | Decay | Living / births | Converted precursor | Funded catalyst | Retained catalyst | Bound food | Usable cell inventory |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.2 | 35 / 11 | 3.9914 | 0.5867 | 0.0340 | 0.0597 | 20.7341 |
| 0 | 0.05 | 39 / 15 | 6.3055 | 0.6143 | 0.2048 | 0.7513 | 21.7082 |
| 0.006 | 0.2 | 48 / 24 | 18.4237 | 1.4348 | 0.1905 | 0.4412 | 29.7833 |
| 0.006 | 0.05 | 48 / 24 | 28.6224 | 1.5248 | 0.6984 | 2.2403 | 32.6344 |

The combined recipe retains about 20.6 times as much catalyst as baseline and
converts 7.17 times as much raw feedstock. It has 0.4145 catalyst mass near living
supports versus baseline 0.0247. Its catalyst effective volume fraction is 8.12%
versus 1.83%: a more extensive material distribution, not proof of coherent
morphology. All populations still have maximum generation one and no deaths;
neither descendant division nor ongoing demographic renewal has been shown.
Fast-decay/mobile-raw has greater fractional overlap with living supports (86.8%
versus candidate 59.4%), so longer memory also leaves more material behind.

This makes the two-parameter recipe a useful separate spatial operating regime.
Keep it for actual material captures and the selective-investment comparison;
do not replace the original living baseline. A fresh capture can use the existing
`NativeLivingConstruction.run` with this recipe and duration 32, requiring no
new engine or artistic geometry. Final actual fields/cells of all four probe
worlds are preserved locally in `outputs/living-patch-recruitment-32/`; the source
snapshots there preserve the exact executed study drivers.

The bounded extension recorded descendant generation, actual material,
conversion, funded requests and patch–resident overlap. Genuine spatial
organization would require persistent structured material coupled to living
activity, not merely more cells or a brighter cloud. A useful result can still
be a rich transient: retained trails, repeated patch use or turnover of local
construction. No old run is resumed or modified.
