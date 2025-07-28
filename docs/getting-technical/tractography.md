# Tractography

This page describes the streamline tractography module in `cardiotensor`, used to trace paths through 3D vector fields of myocardial orientation.

Streamlines represent the continuous direction of cardiomyocyte organization by following the principal eigenvector of the structure tensor field. This process is adapted from diffusion MRI tractography and tailored for cardiac microstructure.

---

## Overview

The tractography module uses:

* FA thresholding to restrict tracing to well-aligned regions
* 3D vector field integration to follow fiber orientation
* Curvature and boundary criteria to terminate streamlines
* Optional binning for speed and scalability

Streamlines are computed using Runge-Kutta 4th order integration of the vector field.

---

## Algorithm

### Seeding

Streamline seeds are placed randomly within a binary mask, typically derived from the FA map:

* Seed points are voxel coordinates where FA exceeds a threshold (e.g., 0.4).
* A fixed number of seeds (default: 20,000) are sampled randomly.
* The number and location of seeds can be configured using CLI.

### Integration

Each streamline is computed using Runge-Kutta 4 (RK4) integration through the vector field.

* Input vector field must be of shape `(3, Z, Y, X)`, where the first dimension is (x, y, z) components.
* Step size is defined in voxel units (default: 0.5).
* Streamlines are traced bidirectionally from the seed (forward and backward).

### Termination Criteria

Tracing stops when:

* FA falls below the user-defined threshold (default: 0.1)
* Turning angle between steps exceeds a threshold (default: 60°)
* The streamline exits image bounds

### Filtering

Streamlines shorter than `min_length_pts` (default: 10 points) are discarded. You can also post-filter based on curvature or anatomical ROIs.

---

## Output

The streamline results are saved in `.npz` format as:

* `streamlines`: a list of arrays, each of shape (N, 3), where N is the number of points in that streamline
* `ha_values`: sampled HA (helix angle) values along each point

These can be loaded in Python or exported to `.vtk` for 3D visualization in ParaView.

---

## Example Command

```bash
cardio-generate config.conf --seeds 20000 --bin 2 --step 0.5 --fa-threshold 0.15 --angle 60 --min_len 10
```

---

## Notes

* Streamlines are computed from the 3rd eigenvector of the structure tensor, corresponding to the myocyte axis.
* FA is computed once from the structure tensor and optionally downsampled for speed.
* Helix angle (HA) is sampled at each point along the streamline and saved for further analysis.

---

## Advanced

* Trilinear interpolation is used to ensure sub-voxel accuracy during integration.
* RK4 integration ensures smoother trajectories compared to simpler Euler methods.
* The code is optimized for parallel execution and can scale to large volumes (e.g., 8000³ voxels).
