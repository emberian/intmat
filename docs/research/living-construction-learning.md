# From instrumented residents to renewing material

The design below led to the implemented bounded scene reported at the end.
The four immobile residents isolate actuator authority and credit timing; they
are not the intended limit of Interior. The next scene uses the same Julia
native schedule in a modest 3D world with motion, growth and division restored,
finite raw-feed budgets, local material sensing and paid construction learning.
The purpose is a working spatial demonstration: cells can move through material
they and their predecessors changed, rather than merely optimize a sealed dish.

A first scene should use a 16³ or 32³ grid and tens of founders, with capacity
chosen by the existing conservative funded-population bound. Keep native motion,
contacts, inherited recurrent/genetic state, catalyst costs and chemistry. Begin
from a separate recipe and new state; do not migrate a historical arena. A short
instrumented preview establishes stable dt, actual turnover/births, source use
and wall cost before a bounded multi-generation capture. Use source addition
and removal ledgers, not a reset that quietly replenishes the world. Save actual
material/cell frames and compact histories for an artwork alongside the study.
A division-free or motionless preview is a diagnostic result, not a completed
renewing-world demonstration.

The descendant contract must be explicit. The current assay clears private
learned state when a slot acquires a new ID, which protects identity but does not
establish that learning should never be inherited. Compare two declared laws:

* **Private reset:** a newborn starts with zero learned action coefficients and
  calibrated local observation state; native inherited genome/recurrent state
  and material made by previous residents still persist.
* **Paid partial parental copy:** during the existing native inheritance hook,
  a newborn requests a fixed fraction of the parent's learned readout. Price the
  exact installed Float32 coefficient changes with the same turnover law and
  debit actual child battery into native work. If the complete copy is
  unaffordable, install the declared zero fallback. Do not duplicate battery,
  precursor, rewards, old action scores or pending credit targets.

Both laws clear the newborn's eligibility queue and return predictor history,
initialize observation filters from its first local observations, and record
parent/child IDs, copied norm, paid work and declined copies. The parent's live
private learning remains its own. Existing native inherited physiology and
neural state are not erased by this sidecar choice. A later separate experiment
could retain a value predictor or mutate inherited learned coefficients, but
those are additional laws and should not slip into this comparison.

Environmental legacy is another route to retention: paid catalyst, bound
substrate and remodeled transport survive a resident's disappearance according
to their actual chemistry. Track this legacy independently of genetic or private
learning inheritance. Useful retention would appear as persistent spatial
opportunities and continuing activity after individual turnover; it is not
established merely by nonzero coefficients, lineage age or a diversity score.

The first scene should carry the credit mechanism selected by the bounded lag
experiment, not combine several untested rules at once. Compare matched finite
funding and initial physical state, preserve failures as outcomes, and report
birth/death/movement, material persistence, actual funded construction, learning
work and complete-step throughput. A small set of independent bounded worlds is
preferable to a long unattended soak. This is the bridge to richer living
matter; neither the instrument nor its metrics are the final artwork.

## First implemented scene: measured motion, division and paid inheritance

The new `experiments/native/living-construction-16.toml` now runs through
`studies/native_living_construction.jl`. This is a 16³ periodic domain of extent
2, initially 24 residents, with native growth, motion, division, contacts and
mutation. The original 23-feature policy is retained; no temporal-filter feature
was added to this first scene. Two finite raw sources supply 0.4 each per time,
with budgets 12.8 each (32-time funding). The conservative funded live-population
bound is 340, below allocated capacity 512. Initial and supplied material remain
in the strict native ledger.

Three matched 16-time CPU runs compared paid half-readout inheritance, private
reset, and private reset with raw unadjusted birth credit. All reached 30 living
residents through six births, with no deaths in this window. Original residents'
mean periodic displacement was 0.435 in the extent-2 domain. Actual funded
secretion totaled approximately 0.476. These observations establish a moving,
dividing and remodeling scene; they do not establish selection, ecological
renewal through death, or an advantage for inherited learning.

The first three births preceded actor calibration. Three later births at nominal
times 11.4, 12.2 and 14.82 copied nonzero learned readouts with installed L1 norms
0.005355, 0.005748 and 0.002680. Their actual child-battery debits were
5.355e−6, 5.748e−6 and 2.680e−6, entered in native work. The reset law copied no
private learned coefficients. Birth counts and total funded construction were
nearly identical in this small comparison; no inheritance benefit is claimed.

The living law also makes birth-credit semantics explicit. Immediately after
native partition, before optional copy work, the hook measures the child's
actual battery-plus-usable-stores U. For each pending queue entry owned by the
actual parent ID, it subtracts that exported U from the entry's starting value.
This cancels the direct partition dip in the parent's eventual return and its
predictor target. It does not add physical resources, value future descendants,
or reimburse growth, copy work, later starvation or death. The once-per-child-ID
journal prevents repeated credit; unrelated or reused identities cannot receive
it. Raw-credit controls leave the target unchanged. In this comparison the raw
law reached a maximum learned norm of 0.0899 versus 0.00843 for adjusted credit,
showing that partition bookkeeping materially changes the learning signal;
that difference is not a fitness result.

Newborn predictor history and queues reset. Inherited action coefficients are
installed only after their complete payment is affordable; otherwise the child
gets the declared zero fallback. Its first local readings initialize observation
differences to zero. An off-cadence first action is deterministic from the actual
inherited readout, then the common decision clock resumes scored exploration.
The native inherited genome, hidden state and material legacy remain active.

The focused inheritance/credit gate passed 20 checks, including unchanged
physical ledgers under target adjustment, parent/child partition U consistency,
identity filtering, one-time handling, recycled-slot pending-credit invalidation and funded/unaffordable copies. Across the
three CPU scenes the maximum conservation residual was 1.57e−13. A separate
four-time CPU/AMD living gate seeded a small nonzero readout solely to exercise
paid copies at early births: three births, maximum checked difference 2.09e−7,
residual 3.91e−14. This parity fixture is distinct from the genuinely learned
coefficients in the 16-time scenes.

Actual frames carry native cell positions/radii/IDs/parents, eight material fields,
source/work/heat ledgers and birth-copy events. The public viewer uses fixed scales
and exact recorded frames. `docs/figures/native-living-construction-comparison.json`
retains all three measured cell/event/ledger histories and material totals;
`native-living-construction-gpu.json` records the native gate. Exact executed scene
and recipe copies are in `studies/outputs/native-living-construction/`. Full private
states and raw arrays remain outside Git; the public field replay is compressed.

A fresh 32-time copy/reset pair used the same recipe and finite budgets. Its
first 17 measured frames exactly equal the corresponding preserved 16-time
runs. Both reached 35 residents through 11 births, without deaths, and still had
maximum generation 1: all births were founder-to-child. The copy law completed
eight nonzero paid copies. Retained catalyst and bound material remained at
0.03399 and 0.05973, respectively, at time 32. Funded secretion increased only
from 0.47563 at time 16 to 0.58668 at time 32. This establishes continued paid
activity and material persistence but does not yet demonstrate learning-state
transmission through multiple descendant generations. The next question is
physical flux and actuator authority, not whether a longer unattended run can
produce a favorable scalar. See `native-living-construction-32-comparison.json`.
