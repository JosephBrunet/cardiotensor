# Tractography

Tractography in `cardiotensor` is the process of generating streamlines that follow local myocyte bundle orientation within a 3D vector field. It enables visual and quantitative analysis of cardiac fiber organization.

---

## Principle

Streamlines are computed by tracing the path along the **principal orientation vector** (3rd eigenvector $\vec{v}_1$) from the structure tensor. Integration follows the local direction of this vector across voxels.

This method is analogous to diffusion MRI tractography, but uses structure tensor-based orientation fields derived from high-resolution imaging.

---

## Algorithm Overview

### **Seeding**:

   * Seed points are placed uniformly in a user-defined volume or mask.
   * Seeding can be restricted based on FA or binary masks.

### **Integration**:

   * Streamlines are generated using numerical integration (e.g., Euler method).
   * Step size and number of steps control streamline length and smoothness.

### **Termination**:

A streamline stops when:

   * FA is below a threshold (default 0.1)
   * Angle between successive steps exceeds a curvature threshold
   * It leaves the image bounds

### **Filtering**:

   * Short or low-quality streamlines can be removed based on length or curvature.

---

## Parameters

* `--seeds`: number of seed points
* `--bin`: Bin the vector field and angle values
* `--step`: integration step size (in voxel units)
* `--min-fa`: minimum FA to start/continue
* `--max-angle`: curvature constraint (in degrees)
* `--mask`: optional binary mask to restrict seeding

Example CLI usage:

```bash
cardio-generate conf.toml --seeds 10000 --bin 2
```

---

## Output

Streamlines are saved in `.npz` format as:

* `streamlines`: list of 3D coordinate arrays
* `HA`: helix angle per vertex (optional)

They can be visualized using tools like:

* **Fury** (Python-based 3D visualization)
* **ParaView** (via VTK export)

---

## Applications

* Visualizing transmural fiber architecture
* Identifying laminar structures and myocardial disarray
* Comparing healthy vs pathological fiber orientation
* Enabling tract-based statistics or mapping

---