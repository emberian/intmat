# Inner shore

This photograph enters the lower interior of *What the light persuades*. Instead of viewing the complete specimen against a distant background, the camera sits close to its recorded absorption structures. Surfaces cross the frame and extend beyond its edges; depth, transmission and shallow focus organize the image.

The source is the unchanged native 256³ optical-feedback state at step 480, t≈12. Its retained mineral contour at 0.32 forms translucent green walls. The absorbed-power contour at 0.25 forms pale sheets; this captured illumination field was held from the preceding refresh at t≈11.75. Both are extracted directly from saved Float32 arrays, without procedural displacement, fabricated texture, geometry warping or another simulation. The small ray accents come from the saved 550nm Julia GRIN observation. The existing experiment's energy balances, controls and fine outgoing-path sensitivity still apply; see [the source study](illuminated-memory.md).

Camera position, lens, focus, area lights and surface optical materials deliberately change. Their colors and Cycles refraction are observation choices, not claims about the mineral's calibrated optical constants or luminescence. The composition is therefore a new artwork from an existing state, not evidence of additional material evolution.

`scripts/render_inner_shore.py` reuses the established geometry/ray renderer and installs a separate macro observation hook. Each new render pins its observation script beside the image; final metadata includes script hashes, source dose hash, camera and contour settings. The publication manifest pins the original producer manifest, solid, dose and GRIN observation. Full PNG and Blender scene remain outside Git; the reviewed WebP is the publication image.

The final is 2000² at 128 Cycles samples on Metal. Q95 WebP is 183,970 bytes. No duplicated replay or film is included; the source artwork retains the interactive field exploration.
