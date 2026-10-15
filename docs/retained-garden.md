# The retained garden

A cutaway through a living chemical interior, photographed at simulation time 30. Cells pay to secrete catalysts; those catalysts bind food and change transport. Jade and wine surfaces distinguish the two retained chemistries, while faint glass spheres and fine bonds show the actual cells and their connections.

The source is the completed `density4-full64` run described in `docs/constructive-volumes.md`: a physical 12 × 12 × 12 periodic domain, native 64³ chemical fields, 320 initial cells, and 61 observations from time 0 to 300. The study source is pinned at `79fac8f054aacb6dca9b02eb79a1aff15ecc20976790537edee4e7f30877a6f4`; its profile hash is `0c4b27470204580bf7a7732683523ae60297295774d5298659261c7d973135f6`. All 427 raw field files were checked against their captured hashes before publication preparation. The later checkpoint continuation is a separate trajectory segment and is not mixed into this sequence.

At the chosen time 30 there are 2,085 cells, catalyst mass 393.60 and retained food mass 1,115.55. By time 300, the original run has 1,761 cells, catalyst mass 92.20 and retained food mass 187.89. The observed growth, depletion and persistence are one trajectory under this law, not a claim of stable homeostasis or evolved organismal anatomy. Soft contact permits overlapping spherical bodies; the optical renderer preserves their actual radii.

## Optical choices

`scripts/render_constructive_volume.py` displays the concentration contour `retained_a + retained_b = 0.45`, colored by the actual local ratio `A/(A+B)`. A-rich material is darker, rougher jade; B-rich material is clearer wine-colored glass. These materials make neighboring chemistries distinguishable. The contour is an observation surface, not a simulated membrane. Native fields are polygonized directly without increasing their spatial resolution.

The cutaway removes a front-right wedge (`x > -0.30`, `y < 0.20` in the normalized display box). Periodic cell and bond copies are clipped to the original domain. Bonds use minimum-image displacements and stable cell IDs. No sphere positions, radii or motion are invented. The faint cell skin, bond radius, colors, thin film, lighting and transparency are artistic choices. They exert no forces and do not feed back into the simulation.

The renderer freezes the source run manifest and checks the hashes of every consumed raw field. Each image records the geometry and field hashes, renderer hash, contour threshold and fixed camera. The collection wrapper preserves the original run manifest and links its immutable captures. The 2000-pixel, 128-sample photograph was visually reviewed and curated; its original image hash is `f5f72b8f0c4e5955215f7482660307acbde3f1e7789f1d90a300de8bc073a8ed`. The compact gallery photograph retains the original locally. The 61-frame actual-state film remains in preparation and has not yet been approved.

Two separate late-continuation observations at times 562 and 762 use the same threshold, camera, exposure and materials. They show sparse retained-food pockets among a substantial surviving cell population. This optical comparison does not imply that all late catalyst is surviving early residue: the independent turnover decomposition described with the constructive-volume results measures continued replacement. These observations remain separate from the original time 0–300 film.
