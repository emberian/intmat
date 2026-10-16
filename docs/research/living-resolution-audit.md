# Physical refinement of the living material fixture

The16³,32³ and64³ checks preserve the same physical extent2³, founder positions and genomes, exchange radius.18, body density128, motor speeds, reaction rates, raw/food diffusion.006, material sensor scales, and finite12.8-per-channel source budgets. Only grid shape changes. This is a refinement of one physical law, not an enlarged habitat or additional funding.

`native_constructed_garden.jl` converts initial concentrations to voxel amounts by multiplying voxel volume. Its raw niche mask samples the same sinusoidal physical profile at voxel centers. Source rates are world totals: each per-site rate is divided by the mask sum, while budgets already denote total material. Opposed masks have equal totals and retain the existing4-time exchange schedule. Source rates therefore remain .4 per channel at every resolution.

The native habitat keeps this distinction through evolution:

- Material sensing and regulatory sensing sum support-weighted amounts divided by voxel volume, then apply fixed saturating receptor scales.
- Catalytic release uses enzyme amount times `raw/(half_saturation*voxel_volume+raw)`; decay acts on amount. Total paid secretion is distributed using normalized physical-support weights.
- Sorptive capacity is an amount proportional to enzyme amount. The sorption law receives voxel volume explicitly. Transport barrier depends on enzyme concentration.
- Food and raw diffusion use physical spacing in their finite-difference fluxes. Physical cell exchange support and mechanical parameters retain their lengths.

The bounded CPU audit initializes each resolution, records actual material receptors, advances.2 model time and checks both ledgers. Founder positions agree exactly. Initial material is59.36000052,59.36000038,59.36000032 respectively; the small difference is Float32 niche-mask/import rounding, not extra supply. Raw totals remain6.4 per channel and source budgets12.8. Maximum initial receptor difference is.005075 from16→32 and.000195 from32→64. This demonstrates consistent units and improving initial quadrature; it does not establish long-trajectory convergence in a nonlinear learning/dividing world.

The actual coupled timestep limits are .25,.09765625,.02441406. The existingdt=.02 passes all three, with less diffusion margin at64³. The materialized support storage grows from884,736 to2,097,152 to11,239,424 bytes at capacity512 because the same physical exchange ball spans more sites. Dense chemistry additionally scales with voxel count. These are real representation costs, not more biology. The short audit timings include first-call compilation and are not throughput evidence.

All three short trajectories close material/potential accounting within3.56e−14. The executed audit is `studies/living_resolution_audit.jl`; compact measurements are `docs/figures/native-living-resolution-audit.json`. Subsequent higher-resolution scenes should retain these physical parameters and budgets, label their grid, and avoid treating improved visual detail as new biological capability.

The subsequent32³ scene runs64 model time with captures every4, preserving the same25.6 total external material. It ends with48 residents, no deaths and maximum generation1, just as the16³ candidate;18 nonzero paid copies occur. Its actual graph has83 bonds at32 and105 at64. Retained catalyst is.659812→.346971 and bound precursor2.156775→1.195493 over that interval; cumulative funded secretion is1.488238→2.009731. Maximum closure residual is6.54e−13. These distinct measured trajectories support the qualitative finite-feed persistence observation, not a convergence rate or inheritance advantage. Compact cell/ledger/event histories are in `docs/figures/native-living-patch-grid32-comparison.json`; exact observer/recipe are in `studies/outputs/native-living-patch-grid32/`. The raw64.9MB field report remains outside Git for rendering.

The native graph also organizes beyond the population count. In32³ the largest connected component contains16 residents at16,26 at32 and46 of48 at64 (final component sizes46,1,1). In16³ the corresponding largest component is11,34,40 (final sizes40,6,1,1). These counts use only the recorded ID-validated bonds. They describe formation of a bonded assembly while material memory persists and decays; they do not attribute organization to learned policy or establish resolution-independent topology.

The 64³ scene subsequently completed 64 model time, captured every 8, using the
exact-gated empty-chemistry kernel
(`856670a8279e22c456012c8984036a3c6f0c981c95ca686be9d29653f0ecf39a`).
That optimization changes no state or physical law. The scene again has 48
residents, no deaths and maximum generation 1; 18 nonzero paid copies occur.
Native bonds grow from 89 at time 32 to 108 at time 64, with the largest component
38→46. Catalyst falls .642630→.325291 and bound precursor 2.093079→1.119267,
while cumulative funded secretion grows 1.458003→1.946225. Supply remains 25.6
and the maximum closure residual is 1.08e−12. These are separately measured finer
fields, not resampled 32³ images. Raw 266.4 MB capture data remain outside Git;
compact histories and executed sources are
`docs/figures/native-living-patch-grid64-comparison.json` and
`studies/outputs/native-living-patch-grid64/`. Its 389-second local run is not an
isolated performance benchmark.
