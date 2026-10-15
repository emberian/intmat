# Interior

Build composable, evolving matter in Julia: activity inside, remodeling, or
constituting material. Ecology is one substrate, not the whole project. Preserve
the original ambition of rich 3D interiors, persistent memory, changing topology,
waves, chemistry and optics. Make experiments and artworks, not just design notes.

Design for performance from the representation onward. Maximize useful throughput
across available hardware: concrete types in hot loops, preallocated state,
device-resident work, fused passes and appropriate numerical algorithms. Measure
complete simulation steps and sustained workloads; distinguish single-trajectory
speed from aggregate replica throughput. A controller-only benchmark is not an
arena speedup. Profile dense contacts as well as sparse populations.

Keep the simulation in Julia. GPU capability is an implementation requirement,
not a reason to duplicate algorithms in Python/Torch. Use the verified hbox Julia
environment in docs/hbox-runtime.md and existing host resource controls. Coordinate
GPU work and shared-host capacity so contention does not masquerade as progress.

Existing long runs have immutable source and independent state. Never silently
migrate, reseed or reset them. New implementations and new artistic laws may be
explored in separate runs. Short successes and failures do not bound the space of
possible behavior; distinguish a machine check from an ecological conclusion.

Commit and push working snapshots. Unsigned snapshots are authorized for the
current batch. Preserve unrelated files, especially dickson11.png and historical
outputs. Keep private transcripts, keys, large baked states and Blender scenes
outside Git; publish reviewed images and compact replay data through build_pages.py.
