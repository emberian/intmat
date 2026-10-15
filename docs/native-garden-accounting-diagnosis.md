# Native garden accounting stop: read-only diagnosis

The exact `4b747919ec56` garden stopped at step 659,992 / time 659.992 with population 93 when the existing material guard failed. Material residual was **+0.2473634**, above tolerance **0.2064319**. Potential residual **−0.1310325** remained within tolerance **0.6874878**. The final checkpoint SHA256 is `171707c2a7a7e0a614e7e8d509f88a04264027591b33ca2869729d9c6a385716`. No state was resumed, no tolerance changed, and no GPU work was used for this diagnosis.

The append-only checkpoint history places the sign change between times 320 and 442:

| Step | Material residual | Potential residual | Population |
|---:|---:|---:|---:|
| 93,384 | −0.005404 | +0.001241 | 259 |
| 198,280 | −0.011722 | +0.003734 | 44 |
| 320,432 | −0.027060 | −0.030525 | 43 |
| 442,224 | +0.094594 | −0.022421 | 43 |
| 557,424 | +0.167250 | −0.065760 | 50 |
| 659,992 | +0.247363 | −0.131032 | 93 |

The final interval accumulated approximately **7.81×10⁻⁷ material/step**. The large positive interval already occurred while population stayed near 43, so births alone do not explain the onset.

## Concrete numerical contributor

A read-only CPU mirror evaluated one existing diffusion update on the saved arrays, without storing or advancing the checkpoint. It uses the frozen Float32 stencil arithmetic and timestep 0.001. At step 557,424, precursor endpoint changes summed to −1.2691×10⁻⁷, while zero-potential waste gained **5.9570×10⁻⁷**, for **+4.6879×10⁻⁷ net material**. At the final checkpoint, precursor changes summed to −1.8947×10⁻⁷ and waste gained **4.1537×10⁻⁷**, for **+2.2590×10⁻⁷ net material**.

In both snapshots, the sum of the computed diffusion flux increments was around 10⁻¹³. The discrepancy arises when each voxel independently rounds `value + dt*flux` back to Float32. The final grid contains **713.03 waste mass**, versus **135.73 soluble precursor mass**, so accumulated waste rounding can increase material without adding bookkeeping potential. This matches the observed signs and is a substantial contributor at the measured scale.

This is not complete attribution: CPU compilation can differ from AMD contraction/FMA behavior, and two frozen snapshots do not integrate the entire trajectory. The existing uptake gather and waste scatter also independently accumulate Float32 donor/recipient endpoints; these remain possible contributors. Habitat sorption and catalytic transfer explicitly account actual endpoints and residues, making them a less immediate suspect than the uncompensated diffusion field update.

The next repair should conserve transfers or retain sub-ULP field residuals, with species-wise material/potential accounting and native whole-step tests. Increasing the guard tolerance would obscure this evidence. Exact source-pinned history and the stopped checkpoint should remain immutable.

Evidence: [diagnostic JSON](../benchmarks/2026-10-15-native-garden-accounting-diagnosis.json). Raw append-only logs and component summaries are retained in `outputs/native-garden-accounting-20261015/`. CPU reads used one thread, CPUs 2–3, 8 GiB/no swap and 60-second command caps.
