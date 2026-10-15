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
owns all growth and biomass transport. Optional matrix redistribution uses that
same fused stage and its existing buffers. The same native exporter, observations,
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

With optional `matrix_transport>0`, gel also exchanges into free occupancy:

```
D_M / h² * [paidⱼ Mⱼ fᵢ − paidᵢ Mᵢ fⱼ]
```

This mobility is distinct from nutrient/signal permeability P. It represents
paid gel redistribution, rather than smoothing a rendered image or silently
changing immobile matrix in the original law. It conserves M separately during
transport and stops unfunded donors. A newly arriving parcel cannot leave in
that stage. Let Dmax=max(D,D_M): combined incoming B+M is bounded by
`2Dmax·dt·sum(h_axis^-2)·(C−ρᵢ)`, since neighboring B+M≤C and both donor
activations are≤1. Per-species outgoing bounds preserve nonnegative B and M.
Thus the shared capacity bound holds for simultaneous competing inflows of
both species, rather than separate allowances that could overfill a recipient.
The stability limit includes D_M. Its default 0 retains immobile matrix.

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

`test/living_tissue_packing.jl` passes 67 focused checks: donor/free-space growth,
zero-fuel suppression, all six competing neighbor inflows, fully occupied
blocking, repeated exchange, no re-spending arrivals, anisotropic spacing,
material/potential balance, child scheduling, recording, exact restart and paid
basal transport/remodeling, including a large-baseline capacity challenge.
Additional mixed-species checks cover all six incoming faces, separate B/M
inventories, repeated exchange, full-site blocking, matrix stability limits and
funded matrix donors without re-spending arrivals.
`test/living_tissue_exporters.jl` passes 30 composition checks for both original
turnover defaults and packing. Strict packing checkpoints include this new law's
hash and the reused turnover/base/shared hashes; they cannot silently restore
as the original law. No accelerator result is claimed for packing yet.

## An oscillatory finite-band candidate

The first four 64² packing probes are retained as failed or limited regimes.
Capacity alone did not repair the fast-inhibition recipe: its excitation retained
20–57% of fluctuation power beyond .7 of the axis Nyquist. Paid basal renewal
reopened some regions, while matrix mobility by itself filled thin spaces and
made a nearly stationary grid-scale maze. These are separate recipes and records,
not replacements for the original native 128³ specimen.

`living-tissue-packing64-oscillatory.toml` instead strengthens slow negative memory
feedback to 3.5, speeds its time constant to 8 and increases matrix dissolution
to 2.0, retaining finite reserve, paid growth and mixed-species transport. The
64², extent24², dt.005, t0–80 pilot develops moving partitions at finite spatial
scales. Excitation's radial power peak lies at wavelengths 2.18–3.43 from t16–80;
its fraction beyond .7 axis Nyquist grows from .00057 at t16 to .073 at t64 and
.058 at t80. This is a refinement candidate, not a convergence result.

The unchanged legacy absolute .10/.24 observer records zero openings here. A
separate, explicitly declared occupancy review uses 30%/70% of solid capacity
(.27/.63). A site must first exceed .63 before an opening below .27 is counted;
refilling above .63 requires an observed opening. At capture spacing .8, the
coarse pilot records 14,445 openings and 12,509 refillings across 4,088 sites,
with at most five openings at a site. These are sampled solid concentration
cycles, not mechanical cavities or organisms. Initial low-density tissue is
never counted as a prior opening. `studies/living_tissue_review.py` reproduces
this readback and physical spectra directly from the native fields.

The matching 128² recipe retains extent24² and the exact same continuous random
initial field, reduces dt fourfold to .00125, runs 64,000 steps to t80 and keeps
capture spacing .8. Its explicit host `step_budget=64000` extends the exporter's
default20,000 limit for this coordinated refinement, with an upper limit100,000;
there is no change to simulation rates. Resident arrays occupy1MiB and declared
capture budget is128MB. A separate controlled hbox CPU8,9 rate check decides
whether the full run can fit a600s projected body within720s,8GiB/no-swap/nice15.
The immutable overlay's source inventory and host controls are preserved in
[`living-tissue-packing-refinement-results.json`](living-tissue-packing-refinement-results.json).

The CPU refinement completed in169.97s of timed advances and noninitial captures;
the full output occupies59.74MB. Capacity excess stayed zero and relative material
drift was4.29e−6. Restricting initial fine cell averages back to64² gives maximum
N/B/x errors5.96e−8/2.98e−8/7.45e−9. At128², radial excitation wavelength is3.43
at t16/32,2.4 at t48 and2.667 at t64/80. Its late high-frequency fraction stays
below.001; the final solid high-frequency fraction is.0036, compared with.0711
on64². The moving partitions and sampled occupancy cycling survive refinement:
48,278 openings and42,191 refillings across15,174/15,014 sites, at most five
openings at a site. Both retain the explicit.27/.63 occupancy thresholds.

Final spent material differs by+2.46%, biomass by+5.20% and matrix by+3.60%
from the coarse pilot. Refining time and space together does not establish
formal spatial convergence; a small high-frequency tail alone is insufficient.
The retained finite scales and actual material cycling support a separate3D
pilot at the fine physical spacing. They do not justify reusing the coarse
grid as a resolved substitute or reclassifying low-density fluid as empty space.

### Frozen reaction derivation

The new recipe is motivated by a uniform reaction that is stable while a finite
spatial band has oscillatory growing modes. This is an analytical hypothesis
about frozen material/payment coefficients, not the complete material Jacobian,
a discrete-time theorem, a nonlinear pattern guarantee or a steady state of the
finite-reserve law. Spatially varying permeability and material dynamics remain
part of the actual nonlinear experiment.

Let s=B/(.12+B), v=s·paid·N/(.15+N), t=1−s and z be the tanh slope. Capacity .9
implies s≤15/17, hence t≥2/17; 0≤v≤s and 0≤z≤1. In the absence of spatial terms,

```
a=v(2.8z−1)−t, b=2.3vz, c=3.5vz
J=[a −b −c; 2 −2 0; 1/8 0 −1/8]
```

The characteristic polynomial is λ³+Aλ²+Bλ+C, where

```
A=17/8−a                        ≥89/136
B=2(b−a)+(1/8)(2−a+c)          ≥(9/8)v+(17/8)t+1/4
C=(1/4)(b−a+c)                 >0, and C≤v+t/4
AB−C                           ≥1203/18496 >0
```

Indeed a≤25/17, b−a=v(1−z/2)+t and b−a+c=v(1+3z)+t.
The displayed A/B/C bounds give
`AB−C≥−287v/1088+1241t/1088+89/544`; substituting v≤1−t and t≥2/17 yields
the final lower bound. Positive coefficients and AB>C exclude every root in
the closed right half-plane by the elementary cubic argument in
[`living-tissue-mode-tests.md`](living-tissue-mode-tests.md).

Using actual physical diffusion/stencil symbols, the exported frozen dispersion
at B=.15/.26/.5/.9 has negative uniform growth but a growing complex pair at
wavelengths3.13/3.35/3.54/3.68 on the128² axis scan, with oscillation periods
14.91/14.94/16.14/17.82. The B=.1 case has no growing band. These calculations
use N=.6,M=x=y=q=0 and full payment; derivatives of renewal, phase work, support,
material transport and changing payment are intentionally omitted.

### Numerical gate

`studies/living_tissue_device_gate.jl` shares only numerical/protocol checks,
not a simulation law. `living_tissue_packing_gpu.jl` supplies a dense randomized
16³ fixture with mixed B/M occupancy, unfunded donors, anisotropic spacing and
six full mixed-species neighbors competing for an empty recipient. It checks
CPU/device parity, nonnegative bounded material, capacity, recording, strict
checkpoint restore and exact same-device11-step continuation before a requested
bake can start. Its CPU plumbing passed with zero parity errors and exact
continuation. That plumbing is not an accelerator result.

The subsequent immutable `e2a2f16` AMD gate passed with maximum state error
5.96e−8, reserve error0 and exact checkpoint/11-step continuation. A128³ short
probe measured three synchronized256-step batches in.36128/.36480/.36109s;
actual native eight-field captures including transfer and diagnostics took
.42466/.07411/.08608s. Taking the slowest step/capture plus20% margin projects
129.82s for64ksteps and40 noninitial captures. This is a short rate projection,
not a sustained bake timing. The first host wrapper's missing capture parent
directory is retained as a failed probe; a fresh wrapper fixed the parent only.
Source hashes, strict gate, timing conventions and host controls are recorded in
[`living-tissue-packing-amd-results.json`](living-tissue-packing-amd-results.json).

The approved full native128³ run then completed t80 in99.574s of timed body
work (64,000 whole steps and40 noninitial native captures) and119.458s for the
complete controlled command. Initial capture, compilation/warmup and final
checkpoint/record/maps are excluded from body time but included in command
time. Persistent device field buffers total134,217,728bytes; this is not a
measured peak VRAM value. Output before the small host execution record is
2,921,497,079bytes, below3GB, with41 native states every two physical time units.
Physical spacing is.1875 in all axes and dt.00125 is below the3D stability bound.

All sampled values remain finite/nonnegative, relative material drift is3.94e−6,
and maximum solid-capacity excess is one Float32 unit5.96e−8, within the declared
roundoff allowance. Final excess is zero. Of the initial165,888 integrated
reserve,132,067.18 remains;33,156.50 material is spent. This is a finite budget,
not an imposed activity movie. Fixed-threshold volume review is separate from
this numerical completion evidence. The unchanged legacy.10/.24 observer
records no openings and is not reinterpreted as the declared.27/.63 readback.

The single-parameter64² `oscillatory-immobile` comparison sets only matrix
transport.04→0. Both runs retain finite scales and repeated solid-space cycling;
matrix mobility is not established as the origin of oscillations. It changes
persistent matrix and future communication: final integrated M27.86→25.16 and
mean permeability.33179→.43879. These are coarse exploratory comparisons,
retained in the packing result record rather than claims of converged causality.

The subsequent whole-volume readback is retained in
[`living-tissue-packing-volume-review.json`](living-tissue-packing-volume-review.json).
It reads exact native3D quantiles and spectra at t20/40/60/80 and exports
unresampled orthogonal Float32 planes, using the same fixed optics for review.
Excitation's radial wavelength is3.43 at t20/40/60 and2.667 at t80; high-frequency
power beyond.7 axis Nyquist is.00025/.00072/.00044/.00053. Solid-density tails
are.0021–.00324. This supports finite-scale3D structure, not formal3D convergence.
Fixed solid threshold.45 leaves lower-density fractions.000187/.42593/.43170/
.41653 at those times. These are concentration regions, not mechanical cavities.
At t40, solid q10/50/90 is.1875/.6087/.9000 and matrix is.00514/.03210/.15911;
at t60 they are.2189/.5343/.9000 and.00871/.03246/.11323. Sampled center-site
excitation upcrossings imply periods12.4–14.5, about six to seven retained states
per cycle. Optical threshold choices remain separate from these native measures.
Across all41 retained native states, the declared.27/.63 observer records
5,857,712 openings and5,130,463 refillings;87.54% of sites open and86.75% refill,
with at most five observed openings at a site. The two-unit cadence can miss
short events. These are sampled concentration cycles, not organism counts.
