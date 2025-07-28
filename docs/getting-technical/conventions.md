# Conventions

## Sign conventions

Cardiotensor processes full-field orientation data in 3D volumes. Consistent with NumPy matrix indexing, the position (0, 0, 0) refers to the corner at the top-left of the first slice. 

This convention is maintained across all processing pipelines, including structure tensor calculation and streamline generation.

## Orientation fields and eigenvectors

Eigenvector maps are saved in the shape `(3, Z, Y, X)`:

- 3 = x, y, z vector components
- Each voxel contains the orientation of the 3rd eigenvector (principal myocyte axis)

These vector fields are used for computing helix/intrusion angles, streamlines, and for visualization.


## Units of measurement

- Lengths are expressed in **pixels** (or voxels for 3D).
- Angles (helix, intrusion) are in **degrees**.
- Fractional anisotropy is **dimensionless**, ranging from 0 (isotropic) to 1 (highly anisotropic).


## Image format

All image volumes and masks must be provided as image stack files in formats such as TIFF JP2, PNG, or as single RAW file (MHD).

Internally, images are processed as NumPy arrays:
```python
volume.shape == (Z, Y, X)
```

The same format is used throughout the orientation pipeline.


---
