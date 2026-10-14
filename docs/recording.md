# Record matter, then photograph it

`record_state(world, processes, directory; step)` captures the simulation between
complete steps. It transfers arrays to the host and writes a fresh directory with
`state.json` and native-precision floating binaries. Recording has no camera, transfer function,
render exposure, or integration timestep. Change capture cadence while keeping
the same `step!` schedule and the final world and persistent process state are
identical. The test exercises this with recurrent fields and moving carriers.

The manifest stores all eight world fields, simulation time and step, grid shape,
physical extents and spacing, cell-center coordinates, backend, world seed, Julia
version, and SHA-256 hashes of the implementation sources. Binaries use the
manifest's byte order and Julia column-major layout (x varies fastest for grid
arrays). Each binary has its shape and source element type in the manifest.
Schema version 2 preserves Float16, Float32 and Float64 arrays exactly, using
`.f16`, `.f32` and `.f64` filenames and corresponding dtype metadata. Other
floating types fail explicitly. Version 1 always converted arrays to Float32;
existing version 1 records retain that conversion.
Particle IDs remain exact integers in JSON.

Neural records include both persistent latent channels, scratch arrays, effective
recurrence weights, biases and feedback/rate parameters. Its constructor does
not retain the initialization seed or mode; the manifest says this explicitly
instead of guessing. Particle records include positions, hidden state, IDs,
storage, all rule matrices and scalar rates. Positions are periodic cell
coordinates: physical position is `position * spacing - extent / 2`. Columns
identify particles, with the same ordering across their state arrays. All other
built-in process parameters are stored in process execution order. Unsupported
processes fail explicitly instead of disappearing from the recording. New process
types implement `record_process_fields(process)` returning a NamedTuple of every
persistent array and effective parameter; the hook must be read-only. Scalars,
tuples and numeric arrays are supported.

The destination must be new. Files are prepared in a sibling temporary directory
and moved into place after the manifest is written; failures remove only that
new temporary export. Call recording when simulation work is quiescent: it is
not a concurrent snapshot or a restart checkpoint, and there is no restart
loader. Scratch arrays are included for inspection, but are overwritten by the
integrator and are not persistent rule memory.

The existing renderer samples its optical arrays in normalized `[-1,1]^3`
coordinates. This is a separate optical representation, not the state archive's
physical coordinate system. To preserve a world aspect ratio, an optical
consumer must map physical extents to its camera geometry explicitly. The raw
state remains independent of that choice; `optical_fields` supplies a subsequent
material-to-optics interpretation rather than an additional simulation step.
