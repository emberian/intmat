# From photographs to an instrument

Keep one Julia simulation implementation. Let a browser display its observations
and submit bounded experiments. The same interface can open a recorded experiment
on GitHub Pages or connect to a local Julia worker. A public remote worker is a
separate service, not part of Pages hosting.

## Already available

- `scripts/build_viewer.py`: standalone WebGL2 recorded volume exploration.
- `web/native_observatory.html`: actual 3D fields, cells, bonds, identities, slices,
  time selection and optical mappings.
- `web/arena.html` + `scripts/serve_arena.py`: local acknowledged pause, checkpoint
  and parameter commands for the older arena runtime.
- Native source-aware checkpoints and independent constructive-world batches.
- `site/laboratory/`: recorded controller capability matrix, per-founder memory
  curves and matched controls. This is an analytical view, not live simulation.

The gap is a common native control/observation adapter and a workbench that joins
these pieces. The older arena control loop does not yet control the native garden.

## Central interaction

**Select a state → branch it → change one mechanism → compare both futures.**

The parent remains immutable. A branch records source, checkpoint, parameters,
interventions and its finite step/time/resource budget. Sliders for display optics
respond locally. Sliders for physical laws create a branch or an acknowledged
step-boundary command; their effect is never represented by silently relighting
old data. Branching from an observation requires its corresponding complete
checkpoint—an averaged volume capture is insufficient.

Keep high-resolution state on its Julia device. Publish bounded, selected field
channels, cell records and scalar diagnostics at a separate display cadence.
Transfer block means or selected slices for responsiveness; request larger regions
or a high-quality photograph explicitly. Asynchronous display needs immutable or
double-buffered snapshots so GPU steps cannot race with an export. Backpressure
may skip display frames, never simulation steps or intervention acknowledgements.

The worker admits a finite vocabulary: describe parameters, observe, step, pause,
checkpoint, fork and apply a declared intervention. It must not evaluate arbitrary
Julia sent by a webpage. Local service binds loopback; a remote service needs its
own authentication, origin checks and resource scheduler. Public recorded viewing
works without access to compute. Existing artifact and host limits still apply.

## Technology choice

Reuse our WebGL viewer and Julia worker first: it preserves current rendering and
does not introduce another plotting stack into the hot path. The network/UI layer
may use ordinary web tools; the numerical laws remain Julia. Bonito/WGLMakie is a
credible alternative for integrated Julia-driven controls and scientific plots;
GLMakie is a credible native desktop front end. Benchmark transfer and rendering
latency before choosing a replacement for the existing volume renderer.

[GitHub Pages](https://docs.github.com/en/pages/getting-started-with-github-pages/what-is-github-pages)
hosts static HTML/CSS/JavaScript. [WGLMakie](https://docs.makie.org/stable/explanations/backends/wglmakie.html)
supports browser visualization and static export, with Julia-dependent callbacks
requiring a live connection. [GLMakie](https://docs.makie.org/stable/explanations/backends/glmakie.html)
provides native OpenGL windows. Julia-to-WebAssembly is a separate compiler/runtime
effort; this plan does not depend on porting AMDGPU or Metal kernels into a browser.

First live vertical slice: two fresh short native branches, one receptor edit,
fixed field scales, same-ID response plots, acknowledged completion and a saved
comparison. Then generalize the adapter to other material processes.
