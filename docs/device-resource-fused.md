# Experimental fused resource credits

This is a CPU-validated Julia study, not a production replacement or a measured AMD speedup. It combines source/sink exchange with the first deterministic paired reduction, and repairs the inherited rounded endpoint subtraction. Production source files are unchanged.

`studies/device_resource_fused.jl` derives its voxel helper from the source-pinned authoritative resource helper. The explicit changes are error-free Float32 endpoint differences (high plus low) and a comparison of that pair against the request before accepting an endpoint. For example, before 1 plus a requested 2²⁶ becomes the representable field 2²⁶: the actual credited increment is **2²⁶−1**, instead of the old rounded grant 2²⁶. A rare source rounding case also receives a smaller endpoint when the exact increment would exceed its request. This is a declared finite-precision repair, not pure scheduling equivalence.

Each `(block,species)` thread updates up to 128 strided sites, accumulating source and sink pairs directly into the existing first-level partials. Later reduction levels and the compensated budget finish are reused. A carrier intercepts the authoritative phase's grant launch and its now-redundant first reduction; it does not duplicate the complete cellular schedule. The existing bounded compensated-sum policy still applies: this does not promise arbitrary exact sums across unlimited dynamic range.

CPU evidence:

- **3,099 assertions passed**: 3,000 comparisons with 512-bit BigFloat; 54 endpoint, clamp and budget checks; five near-empty/NaN/poisoned-grant checks; 40 complete ecology checks.
- The full checks cover closed/fed sources with both Float32 and Float64 batteries, 20-step paid catalyst/sorption trajectories, study-pinned strict restore and continuation, and the failure latch. Ordinary/fused physical arrays agreed exactly in these fixtures.
- An independent read-only audit passed **626 additional BigFloat edge comparisons**, including subnormal boundaries, unit ULPs, 2²⁶ and floatmax.
- The complete benchmark driver passed a tiny 8³ CPU closed/fed plumbing run with zero field/restart/continuation differences. These timings are not performance evidence.

The experimental checkpoint wrapper pins both the study and base sources. Replaying only the underlying ordinary checkpoint would omit the experimental schedule identity; use `device_resource_fused_checkpoint` and `device_resource_fused_restore` together.

## Memory and native gate limits

The current carrier **still retains the ordinary V×6 grant allocation**. Its poisoned-grant test proves that allocation is no longer read or written, not that memory has been freed. Fusion avoids approximately 48 bytes/voxel of grant write-plus-read traffic. A future ordinary allocator could remove 24 bytes/voxel of dedicated grants; a compact adapter already reusing that arena for mandatory later V×3 transients can save at most another **12 bytes/voxel**, not 24. First-level partial storage remains necessary.

The serial 128-site thread loop may reduce available GPU parallelism. Native performance must be measured rather than inferred from traffic. AMD IEEE/subnormal behavior also needs a hardware gate before claiming error-free endpoint accounting across all tested ranges. The prepared `scripts/gate_device_resource_fused.jl` includes the large-increment endpoint probe, matched whole-loop closed/fed timing, strict restart/continuation and failure checks. No native probe is launched while the long garden has priority.
