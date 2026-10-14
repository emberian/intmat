# Photographing the demixed phase

![Amber phase specimen](../outputs/phase/phase-glass.png)

The specimen is a geometric rendering of the Cahn–Hilliard export in `outputs/phase`. Its interwoven channels come from the saved `structure.f32` field, with the interface selected at a phase fraction of 0.5. OpenVDB extracts the continuous isosurface; the mesh contains 6,022 vertices and 6,038 polygons. There is no procedural detail or emission added to the specimen.

The optical interpretation is a cropped amber slab. The crop uses voxel indices `[2:30, 7:25, 2:30]` in the original 32³ grid, preserving most of its width and height while reducing the depth so transmitted light can reveal channels. One separable periodic smoothing pass with weights `[0.125, 0.75, 0.125]` removes the most obvious voxel faceting before extraction. Zero padding closes the selected phase at the crop boundary; these enclosing surfaces are display cuts, not physical walls predicted by the simulation. The extraction itself uses zero mesh adaptivity.

The simulated field establishes the geometry. Amber color, refractive index 1.46, transmission 0.94, roughness 0.08, absorption, three studio lights, and the camera are chosen for the photograph. They are neither measured material properties nor evidence of mechanical stability. The saved `.optical.json` records the choices and source; the `.blend` contains the inspectable mesh, materials, camera, and lighting.

Reproduce with Blender 5.0.1, including its `openvdb` and `numpy` modules:

```sh
/Applications/Blender.app/Contents/MacOS/Blender -b \
  --python scripts/render_phase.py -- \
  outputs/phase outputs/phase/phase-glass.png --samples 128 --size 1000
```

The interface and optical stage live in a dedicated script so they can be reused after a different Cahn–Hilliard run without replacing the computed morphology with a stock sphere or generated texture.
