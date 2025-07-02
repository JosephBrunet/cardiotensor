# Commands

Cardiotensor provides powerful features for analyzing 3D cardiac imaging data, from processing imaging volumes and calculating orientation to managing large-scale projects with support for multiple configurations and platforms.

## Orientation Computation

Compute myocyte orientation from a 3D volume using a configuration file.

- `cardio-tensor`
  Computes structure tensor, helix/transverse angle, FA, and eigenvectors.

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

- `cardio-visualize-vector`
  Render vector fields using Fury. Optionally export to ParaView (VTK format).

- `cardio-generate-streamlines`
  Generate streamlines from the vector field. Outputs `.npz` files compatible with Amira.

- `cardio-visualize-streamlines`
  Visualize streamlines in 3D using Fury. Can also export to ParaView.

See the [example](../getting-started/examples.md#generating-and-visualizing-streamlines) for usage.

---

## Next Steps

Read the [example](../getting-started/examples.md) for an introduction to each feature.
