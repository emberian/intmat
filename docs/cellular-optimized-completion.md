# Optimized cellular soak: completed

The source-pinned `361d0a9` continuation stopped normally at its unchanged
2026-10-15 21:54:05 UTC deadline: `stopped / wall_time_complete`, after
43,200.010 wall seconds. The Julia process is absent. No resume or reseed was
performed. All four bounded runs covered by the observatory are now handled.

| Terminal record | Value |
| --- | ---: |
| Simulation time / step | 8,497,388.626621 / 84,973,885 |
| Living cells / bonds | 0 / 0 |
| Cumulative births / deaths, including parent history | 42,106 / 42,170 |
| Complete metrics rows / chunks | 2,161 / 9 |
| Archived states | 72 |
| Actual artifact bytes / budget | 74,859,958 / 4,000,000,000 |
| Material / chemical-energy residual | 0.002730803 / 0.016143719 |

The continuation began with 207 cells at t=20,666.500308. It reached an
observed maximum of 231 cells and ancestry depth 208. Its first recorded empty
population was at t=522,081.507780, near the finite supply's nominal funding
boundary. This earlier collapse was already recorded; the remaining interval
advanced an empty trajectory. Cumulative births, deaths, work and heat at the
first empty observation match their terminal values. Empty-state throughput
does not measure the cost of a populated ecology. This history neither
demonstrates open-ended organization nor rules out richer behavior under other
laws or resource regimes.

The [compact record](cellular-optimized-completion.json) retains all-row extrema,
sampled longitudinal metrics, chunk hashes, archive inventory, terminal
inspection and the last live resource-control inspection. It complements the
[original two-soak review](cellular-soak-completion.md). Full raw evidence stays
at `/othersys/tmp/intmat-soaks-optimized-20261015/continuation` on hbox.

The final `checkpoint-1.jls` is 338,453 bytes. Its stable pointer, final metrics
and stopped status all identify step 84,973,885. Two independent read-only
inspections agree on SHA256
`3cf1e53f0fc79d6ff87705a513dc8f907d73ee380ea7f093340bc3d3963765f6`.
All 115 source-manifest files match; parent SHA256 remains
`d7f81c4cbbf21b05ce1918a142764eb0d22e1d6b68565c012c66cd759a6f01de`.
There are no malformed metrics rows, partial tails or failure files. This
terminal check verified bytes and metadata, without deserializing or advancing
the checkpoint.

The 21:30 live inspection verified CPUs 6–7, nice 15, an 8 GiB memory ceiling,
no swap and zero cgroup OOM events. These controls cannot be read from the
terminated process; their last live evidence is preserved explicitly. The
final status's pre-write artifact estimate differs from the actual directory
total by 30 bytes; both are well below budget.
