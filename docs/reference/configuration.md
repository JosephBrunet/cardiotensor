
# Configuration File

The configuration file in **cardiotensor** allows you to specify the parameters required for processing input datasets, calculating orientation tensors, computing angles, and saving results. Below is a detailed explanation of each section and parameter in the configuration file.

---

## Example Configuration File

This is an example of a configuration file as found in the `examples/` directory.

```ini
[DATASET]
# Path to the folder containing the input images (accepted formats: .tif, .jp2, .mhd)
IMAGES_PATH = ./data/635.2um_LADAF-2021-17_heart_overview_

# Voxel size of the input images in micrometers (µm).
VOXEL_SIZE = 635.2

# Path to the folder containing the segmentation mask (accepted formats: .tif or .jp2)
# If no mask is available, leave this field empty.
MASK_PATH = ./data/mask


[STRUCTURE TENSOR CALCULATION]
# Gradient scale (sigma) used for smoothing before gradient calculation.
# This Gaussian filter reduces noise while preserving important edges.
SIGMA = 0.4

# Integration scale (rho) used for smoothing the products of gradients.
# A larger value results in smoother, more coherent orientation fields by integrating over a larger neighborhood.
RHO = 0.6

# Multiple of RHO at which the gradients are truncated.
# A larger value requires a larger margin around the image, set in VERTICAL_PADDING.
TRUNCATE = 4

# Padding to avoid border artifacts
# By default this includes the support of both Gaussian filtering stages.
# VERTICAL_PADDING = 10

# Number of slices to load into memory at a time during processing.
# This affects memory usage and processing speed. Adjust based on system capacity.
N_CHUNK = 20

# Enable GPU computation during the structure tensor calculation (True/False)
USE_GPU = True

# Store temporary eigenvalue/eigenvector volumes in memory-mapped files.
# Enable this when tensor outputs do not comfortably fit in RAM. Temporary
# files default to OUTPUT_PATH/.cardiotensor_scratch, never the system /tmp.
# Unsafe N_CHUNK values are reduced automatically from the available memory.
LOW_MEMORY = False
# Optional explicit scratch root on your data-processing filesystem.
# LOW_MEMORY_DIR = /path/to/data-proposal/cardiotensor_scratch

# Whether to save the orientation vectors (True/False)
# Use for 3D vector/fiber visualisation
WRITE_VECTORS = True

# Vector storage format: npy or zarr
# Use zarr for compressed, chunked storage of very large volumes.
VECTOR_FORMAT = zarr

# Specify the processing direction:
#   - True: Process slices from the beginning (0) to the end.
#   - False: Process slices from the end to the beginning.
REVERSE = False


[ANGLE CALCULATION]
# Whether to save the helical and intrusion angles and fractional anisotropy (True/False)
WRITE_ANGLES = True

# Choose which angle pair to compute:
#   ha_ia  → Helical Angle / Intrusion Angle (standard myocardial architecture)
#   az_el  → Azimuth / Elevation (generic vector orientation in 3D)
ANGLE_MODE = ha_ia

# Coordinates of points along the left ventricle axis.
# The first point should be coordinates of the mitral valve point in the volume ([X, Y, Z])
# The last point should be coordinates of the apex point in the volume ([X, Y, Z])
# Intermediate points will be interpolated to create a curved centre line.
AXIS_POINTS = [104,110,116], [41,87,210], [68,95,162]


[TEST]
# Enable test mode:
#   - True: Process and plot only a single slice for testing.
#   - False: Perform the full processing on the entire volume.
TEST = True

# Specify the slice number to process when test mode is enabled.
N_SLICE_TEST = 155


[OUTPUT]
# Path to the folder where the results will be saved
OUTPUT_PATH =./output

# Output file format for the results (e.g., jp2 or tif).
# Default format is jp2
OUTPUT_FORMAT = tif

# Type of pixel values in the output file:
#   - "8bit" for grayscale 8-bit images
#   - "rgb" for 3-channel color images
OUTPUT_TYPE = 8bit

# Optional colormaps for RGB output and test plots.
# Use any Matplotlib colormap name, or helix_angle for the CardioTensor map.
# COLORMAP applies to both angles unless COLORMAP_ANGLE1 or COLORMAP_ANGLE2 is set.
# For ANGLE_MODE = az_el, defaults are helix_angle for AZ and viridis for EL.
# COLORMAP =
# COLORMAP_ANGLE1 =
# COLORMAP_ANGLE2 =
```

!!! note

    Modify the configuration file as needed to fit your dataset.
    Relative dataset, mask, and output paths are resolved from the directory
    containing the configuration file.

---

## Explanation of Parameters

### `[DATASET]`
- **`IMAGES_PATH`**: Path to the input dataset containing 3D images.
- **`VOXEL_SIZE`**: Voxel size of the dataset in micrometers.
- **`MASK_PATH`**: Path to the binary segmentation mask. Leave blank if no mask is available.

!!! note

    The mask volume can be downscaled compared to the heart volume. The binning factor will be estimated automatically and the mask will be upscaled accordingly.

---

### `[STRUCTURE TENSOR CALCULATION]`
- **`SIGMA`**: Noise scale before gradient computation. Helps reduce noise while preserving structures.
- **`RHO`**: Integration scale for smoothing tensor components. Larger values yield smoother orientation fields.
- **`TRUNCATE`**: Multiple of RHO for defining the gradient filter kernel size.
- **`VERTICAL_PADDING (optional)`**: Extra padding (in voxels) to avoid edge
  artifacts. If omitted, Cardiotensor includes the support radius of both the
  `SIGMA` and `RHO` Gaussian filters.
- **`N_CHUNK`**: Number of slices to process simultaneously.
- **`USE_GPU`**: Enable GPU computation (requires CuPy).
- **`LOW_MEMORY`**: Store temporary tensor outputs on disk using memory maps.
  This reduces RAM use but can be slower unless the temporary directory is on
  fast node-local storage. Cardiotensor treats `N_CHUNK` as an upper bound and
  reduces it when the padded chunk would exceed a conservative cgroup-aware
  memory estimate.
- **`LOW_MEMORY_DIR (optional)`**: Scratch root for tensor memory maps. If
  omitted, CardioTensor uses `OUTPUT_PATH/.cardiotensor_scratch`; it never
  falls back to the system `/tmp`. The required and available scratch
  capacities are printed before calculation. Temporary task directories are
  removed after normal completion. Because OOM and `SIGKILL` cannot run
  cleanup code, inspect this scratch root after failed jobs and remove orphaned
  `cardiotensor_tensor_*` directories.
- **`WRITE_VECTORS`**: Save the orientation vector field.
- **`VECTOR_FORMAT`**: Vector storage backend:
    - `npy`: one `float32` file per Z slice in `eigen_vec/`.
    - `zarr` (default): compressed, chunked and sharded storage in `eigen_vec.zarr/`.
      This is recommended for very large volumes and supports safe slice-level
      restart tracking.

!!! warning

    Orientation vectors are saved in `float32` format. Zarr uses Blosc Zstd
    compression with bit-shuffle; zero-filled masked regions generally compress well.

- **`REVERSE`**: Process volume from end to start if set to `True`.

!!! note

    The structure tensor calculation is performed using the [`structure-tensor`](https://github.com/Skielex/structure-tensor) Python package with optional CUDA support.

---

### `[ANGLE CALCULATION]`
- **`WRITE_ANGLES`**: Save helical and intrusion angles and fractional anisotropy values.
- **`ANGLE_MODE`**: Choose `ha_ia` for cardiac HA/IA or `az_el` for generic azimuth/elevation maps.
- **`PROJECTED_ANGLES`**: In `ha_ia` mode, keep `False` for unprojected 3D HA/IA. Set `True` only to write legacy projected `HA_projected`/`IA_projected` comparison maps.
- **`AXIS_POINTS`**: List of 3D points `[X, Y, Z]` along the ventricle axis. Typically, the first and last points correspond to the mitral valve and apex. Intermediate points help create a curved axis via interpolation.

---

### `[TEST]`
- **`TEST`**: Enable test mode to process a single slice.
- **`N_SLICE_TEST`**: Index of the slice to test.

!!! note

    Use test mode to verify processing on a small subset before running on the full volume.

---

### `[OUTPUT]`
- **`OUTPUT_PATH`**: Directory where results are saved.
- **`OUTPUT_FORMAT`**: Output format (`jp2` or `tif`).
- **`OUTPUT_TYPE`**: Image bit depth and channel type:
    - `8bit`: grayscale 8-bit image.
    - `rgb`: 3-channel color image.
- **`COLORMAP`**: Optional shared colormap for both angle outputs.
- **`COLORMAP_ANGLE1`**: Optional colormap for the first angle output (`HA` or `AZ`).
- **`COLORMAP_ANGLE2`**: Optional colormap for the second angle output (`IA` or `EL`).

---
