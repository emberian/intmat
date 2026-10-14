# A shared film / amber and jade

The photograph uses an immutable copy of the metabolic arena at model time 1698. The full 96×96 phase, shared-film, interface and substrate fields control transmitted light. The 3000 captured inhabitants are a capped observation sample of 5599 agents; each point uses its recorded x/y position and founder lineage. Their small spherical glyphs are chosen for visibility, with no inferred organism shape or added connections.

Jade and amber absorption follow phase composition. Shared film adds copper tint and optical density; substrate availability clears local absorption. The scene has a diffused backlight, restrained studio lights, and a polished shallow optical cut. It adds no simulated-field emission or procedural texture. Thickness 0.08 and point radius 0.004 are display choices in normalized specimen coordinates. The underlying simulation is two-dimensional; no depth dynamics are implied.

The snapshot and hash remain in `outputs/art/metabolic-section-1698/captured-state.json` and `provenance.json`. The Blender scene packs its scalar-derived textures. `shared-film.optical.json` records transfer functions, source hash, counts, camera and optical constants. `scripts/render_metabolic_art.py` reproduces the photograph without accessing the live arena.

The accepted Blender 5.0.1 photograph emits `Maximum number of closures exceeded: 76 > 64` during scene preparation. A bounded simpler-glass probe retained the warning. The source mapping, agent coordinates and material assignments were checked, and the final image was visually inspected; combined shading may be limited, so this artwork is not a quantitative light-transport result. The approved photograph and matching shader remain frozen.
