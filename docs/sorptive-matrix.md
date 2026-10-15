# A habitat that holds chemistry

`SorptiveMatrix` makes constructed catalytic material consequential beyond food
release. The same paid enzyme inventory provides finite binding sites and lowers
extracellular solute permeability. Food attached to those sites is immobile and
unavailable to cellular import or the resource pool's outflow. It becomes
accessible again through desorption or loss of scaffold capacity. There is no
new uncharged material, preferential access by ancestry, or reward for a group.

This is a synthetic material law. It is not a claim that our enzyme material is a
calibrated protein gel. The analytical bet is that making a place changes the
space of ecological possibilities: construction can retain useful chemistry,
delay its availability, obstruct replenishment, and outlive the cells that built
it. The same effects can help a neighbor, feed an exploiter or starve a producer.

## Law and accounting

At each site, each soluble food species has its own binding-site class. With
bound mass B, total local food M, capacity C and voxel volume V, the update solves

    dB/dt = (k_on / V) (M - B) (C - B) - k_off B

exactly for fixed M and C during this split stage. A stable quadratic root and a
separate repeated-root case avoid cancellation in ordinary low-concentration
regimes. The full coupled model still uses first-order process splitting; an
exact local reaction does not make the whole integrator exact.

Capacity is `capacity_per_enzyme` times local enzyme mass. By default A enzyme
binds A food and B binds B; `cross_binding` continuously exchanges this affinity
between the two constructed materials. It does not assign biological roles to
cells. Separate site classes do not compete for a common capacity. Capacity loss
first displaces excess bound food into the soluble pool. If its Float32 field
cannot represent the entire transfer, remaining displaced mass becomes explicit
inert residue and its chemical potential becomes accounted heat. Ordinary slow
desorption can leave an unrepresentable amount bound until a later step.

Bound and free precursor have equal **bookkeeping** potential. This is reversible
kinetic trapping; binding does not create stored work or consume an unrecorded
energy subsidy. Concentration entropy and molecular binding free energy are not
modeled. As elsewhere in this substrate, actual representable field differences
pay for transfers. The arena counts retained material and potential exactly once
through borrowed process references, including when only the matrix wrapper is
passed to the inventory functions.

Matrix relative permeability is

    1 / (1 + transport_barrier * (enzyme_A + enzyme_B) / V).

The stage first reconstructs the existing cellular/background permeability and
then multiplies it by this factor. It never multiplies yesterday's matrix effect
again. The default ordered composition is `(feedstock, soluble, cells, habitat,
matrix)`; current matrix geometry affects the next soluble diffusion stage.
Construction observes existing scaffold geometry immediately when a matrix is
attached. This changes chemical transport, not cell drag, bond stiffness or an
unmodeled deformable wall.

The existing inherited secretion readout controls construction. Cells perceive
free chemistry and their own material/work state; the secretion readout also
senses enzyme and raw feedstock. There is no additional direct retained-food
sensor in this version. Feedback from retained food reaches development through
its subsequent release and uptake. Snapshots expose free food, constructed
catalyst and retained food separately, so a viewer need not confuse them.

## Running and observing

```sh
julia --project=. --threads=2 --gcthreads=1 studies/sorptive_garden.jl \
  experiments/arena/sorptive-garden.toml outputs/new-garden 10000 100
julia --project=. --threads=2 --gcthreads=1 studies/sorptive_garden.jl \
  experiments/arena/sorptive-refuges.toml outputs/new-refuges 30000 250
```

Each study runs two founding seeds with four matched settings: full matrix,
sorption only, transport resistance only, and catalytic function alone. Turning
off binding sets capacity to zero; turning off resistance sets its coefficient
to zero. Other mechanisms and founder randomness are preserved within each seed.
These are ecological/developmental explorations, not inherited-adaptation assays.
Use fresh output directories and pin the source for any long run.

The garden deliberately makes newly accessible food depend entirely on catalytic
construction after its soluble bootstrap. The refuges recipe adds a modest
continuing soluble background, so construction modifies an inhabitable
environment rather than serving as its sole ongoing food gate. Both retain
finite accounted supply budgets. The slower paid motility and longer matrix
lifetime in these recipes give construction time to affect a location; neither
settlement nor continued survival is imposed.

Snapshots add `retained_a`, `retained_b` (bound food mass per physical volume) and
`matrix_permeability` (dimensionless). Matrix metrics include bound potential,
cumulative adsorption/desorption, capacity-displaced food and inert residue.
`scripts/render_sorptive_diagnostic.py` renders those fields alongside catalyst,
free food, actual body disks and bonds. Chemistry uses fixed logarithmic exposure
over the supplied sequence; resistance uses the same linear 0–1 scale in every
frame. Recorded observations are never interpolated in time.

## First longer explorations

The eight garden trajectories under `outputs/sorptive-garden-20261015` completed
10,000 steps each (t=1000), retaining 101 observations and three checkpoints per
arm. Their immutable implementation is in the adjacent `-source` directory.
The full matrix retained peak food inventories of 25.03 and 17.85 mass units in
the two seeds. Peak populations across all arms were 179–247 cells, with no
capacity-blocked births. Connected growth and material retention occurred before
all eight populations eventually disappeared. Last sampled living states ranged
from t650 to t880; these are observation bounds, not exact extinction times.

The first full seed at t20 contains 136 living cells, 256 actual bonds and 6.18
units of constructed catalyst holding 25.03 units of food. It accumulated 154
births across the trajectory. These observations establish consequential
retention, release and transport changes, not evolved cooperation or an enduring
organism. The decline motivates a different ecological regime without rejecting
the architectural possibility. Maximum absolute residuals across the eight arms
were about 2.0e-5 material and 1.6e-5 chemical potential.

The eight separate refuges trajectories reached t3000 with 72–85 living cells
and 304–480 cumulative births. Continued soluble input sustained populations,
but it did not sustain much construction: seven endpoint enzyme inventories
were below 0.00013 and the remaining one was 0.00483. In seed211/full, cumulative
secretion stayed at approximately10.0312 from t100 through t1500, while cells
continued growing and reproducing. Source and observations are retained under
`outputs/sorptive-refuges-20261015` and its adjacent source directory. No historical
arena, garden trajectory or previous source pin is migrated or replaced.

That observation exposed a concrete scheduling asymmetry: the original cellular
growth reaction consumes complementary precursors before the separately
scheduled secretion effector receives them. A requested construction flux can
therefore receive exactly zero, even when growth is active. It is not an
impossibility result—mutations could alter growth demand—but it is a strong
architectural priority that deserves an explicit alternative.

## Shared growth and construction

`ConstructiveEcology` is the opt-in alternative, selected with
`[catalysis] shared_anabolism=true`. It owns the cells and their catalytic habitat;
scheduling either child again is rejected. The old `(cells, habitat)` schedule
remains the default and retains its original arithmetic.

After the cellular controller and uptake/junction stages, the inherited
construction readout proposes A/B secretion. Its work cost joins the respiratory
forecast once; the same cached requests are used after maintenance. Growth and
both construction requests then receive one common scale factor, bounded by
available complementary precursor and available work. No reaction gets an
unconditional first claim. Actual growth/deposition remains donor-bounded, and
unrepresentable construction stays intracellular. Metrics distinguish requested
amounts, planned allocations and actual secreted material.

The construction readouts all observe the same pre-construction external fields;
later cell IDs do not see material just deposited by earlier IDs in that sweep.
Intracellular observations are taken before respiration. Death and division
happen after anabolism; offspring inherit the correct parent readout, while
per-step secretion observations remain attached to the cells that actually
secreted. Catalytic release/decay follows the cellular step and sorption follows
that. This is a new physiological scheduling law, not a migration of old runs.

`experiments/arena/constructive-refuges.toml` combines the shared budget with the
same funded environment as the refuges recipe. Matched mechanism toggles and a
residual-allocation control can be run with `sorptive_study`'s `conditions` keyword
(`full`, `sorption_only`, `transport_only`, `catalytic_only`, `residual_full`).
These comparisons examine continuing material construction and its consequences;
they do not establish evolutionary improvement.

## Completed shared-budget comparison

The ten source-pinned trajectories under `outputs/constructive-refuges-20261015`
now each reach t3000, with 121 recorded states and three full checkpoints. Both
founding seeds run all four matrix settings plus the original growth-first
allocation control. All share the same physical funding and population guard
of 1,536. No arm reaches that guard. Compact results and source hashes are in
`constructive-refuges-results.json`; `comparison.png` in the output directory
shows every arm on common axes.

| Setting | Final cells, seed211 / seed307 | Construction during t1500–3000, seed211 / seed307 |
|---|---:|---:|
| Shared budget, full matrix | 299 / 306 | 92.886 / 81.391 |
| Shared budget, binding only | 295 / 324 | 75.018 / 86.108 |
| Shared budget, resistance only | 315 / 315 | 93.529 / 84.334 |
| Shared budget, catalytic function only | 299 / 316 | 62.769 / 84.809 |
| Growth first, full matrix | 79 / 85 | 0.00358 / <0.000001 |

The shared variants finish with 4.31–9.44 units of catalyst and 213.32–263.19
units of cumulative paid construction. Their late intervals also contain
247–333 births and 224–328 deaths: continuing construction accompanies turnover,
rather than a frozen founder arrangement. Seed211 shows several large rebuilding
and population episodes; seed307 has smaller fluctuations around continued
construction. The growth-first controls finish with approximately 4.35e-5 and
6.88e-10 catalyst units. Maximum observed material and chemical-energy residuals
across all ten arms are 7.80e-4 and 6.03e-5.

This comparison concerns the complete alternative physiology, including shared
anabolic allocation and construction in respiratory work forecasting. It does
not isolate either change, rank the ecological value of binding versus transport
resistance, or establish inherited adaptation. Continued soluble feeding remains
present. The new raw-food-only recipes remove that ongoing soluble support in
separate runs; neither those experiments nor their interpretation replace these
recorded trajectories. A large contact component alone is not a cooperative
organism, and repeated bond observations do not establish uninterrupted contact.

## Checks

Forty-three focused assertions cover binding bounds, equilibrium and semigroup
composition; owned mass/potential; release after capacity loss; initial and
repeated permeability composition; forbidden array aliases; complete arena
observations; and exact checkpoint continuation. Fifty-seven existing catalytic
checks also pass. An independent 120-case comparison of the analytic binding
update against small-step RK4 found maximum normalized difference about 2.1e-15.
A separate-process restore preserved the actual shared habitat/cell/pool aliases
and resumed exactly under a frozen source snapshot. These checks establish the
implemented contracts, not an ecological conclusion.

Fifty-five additional shared-allocation checks cover simultaneous material/work
scarcity, work forecasting, zero donors, simultaneous death/division with reordered identities, terminal-cell cleanup, duplicate-owner
rejection, complete observations and exact restart. Independent comparison against
the previous unrefactored catalytic loop preserves its observations exactly.

One separate developmental-law question remains open: the existing controller
reconstructs preactivation through `atanh(clamp(tanh(raw)))` before subtracting
adaptation, which caps strong inputs. A future explicit variant could use raw
preactivation directly. This is an analytical possibility, not a silent change
to the present developmental dynamics or a claim that they cannot be rich.
