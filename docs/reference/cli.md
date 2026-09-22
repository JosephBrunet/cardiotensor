# CLI Commands

Cardiotensor provides powerful features for analyzing 3D cardiac imaging data, from processing imaging volumes and calculating orientation to managing large-scale projects with support for multiple configurations and platforms.

## Orientation Computation

Compute myocyte orientation from a 3D volume using a configuration file.

- `cardio-tensor`
  Computes structure tensor, helical/intrusion angle, FA, and eigenvectors.

- `cardio-tensor-slurm`
  Submits chunked `cardio-tensor` runs as SLURM array jobs.
  See [SLURM Launcher](../advanced/slurm-launcher.md).

See the [example](../getting-started/examples.md) to get started.

---

## Transmural Analysis

Plot angle profiles across the heart wall using an interactive GUI.

- `cardio-analysis`
  Define transmural lines, adjust sampling, and export results.

See the [example](../getting-started/examples.md#visualizing-transmural-profiles) for details.

---

## 3D Visualization

Visualize results in 3D using vector fields and streamlines.

!!! note

    To use these modules, you must first compute and save the eigenvector field by setting `WRITE_VECTORS=True` in the configuration file **before running the structure tensor calculation**.



- `cardio-visualize-vector`
  Render vector fields using Fury. Optionally export to ParaView (VTK format).

- `cardio-generate-streamlines`
  Generate streamlines from the vector field. Outputs `.trk` files compatible with Amira.

- `cardio-visualize-streamlines`
  Visualize streamlines in 3D using Fury. Can also export to ParaView.

The FURY backend uses FURY 2 (PyGfx/WebGPU). In the streamline viewer:

- `O` enables the clipping plane; `H` shows its transform handles; `I` flips it and `R` resets it.
- `C` enables the crop box; `G` shows its transform handles; `X` resets it.
- Drag gizmo arrows to translate, spheres to rotate, and cubes to scale along an axis. Click the center to change the gizmo coordinate system.
- `L` opens material controls, `+`/`-` changes thickness, and `1`–`4` selects camera views. `5` switches projection.
- `Ctrl+S` saves the session; `P` saves a PNG, `Shift+P` saves at twice the pixel resolution, and `V` records an orbit.
- `B` changes the background and `S` toggles a ruler in dataset coordinate units.

Existing VTK-era sessions are read and converted when saved. New sessions use
format version 2 and are intended for the FURY 2 viewer. Overlays and gizmos have
changed appearance. FURY smoothing now uses a natural cubic spline; set
`--spline-subdiv 1` to preserve the original points. The PyVista backend retains
its VTK smoothing.

See the [example](../getting-started/examples.md#generating-and-visualizing-streamlines) for usage.

---

## Next Steps

Read the [example](../getting-started/examples.md) for an introduction to each feature.
