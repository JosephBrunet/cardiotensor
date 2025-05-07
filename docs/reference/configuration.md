
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
SIGMA = 1
# Integration scale (rho) used for smoothing the products of gradients.
# A larger value results in smoother, more coherent orientation fields by integrating over a larger neighborhood.
RHO = 2
# Number of slices to load into memory at a time during processing.
# This affects memory usage and processing speed. Adjust based on system capacity.
N_CHUNK = 20
# Enable GPU computation during the structure tensor calculation (True/False)
USE_GPU = True
# Whether to save the orientation vectors (as .npy) (True/False)
# Use for 3D vector/fiber visualisation
WRITE_VECTORS = False
# Specify the processing direction:
#   - True: Process slices from the beginning (0) to the end.
#   - False: Process slices from the end to the beginning.
REVERSE = False


[ANGLE CALCULATION]
# Whether to save the helical and intrusion angles and fractional anisotropy (True/False)
WRITE_ANGLES = True
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
```

!!! note

    Modify the configuration file as needed to fit your dataset.

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
- **`N_CHUNK`**: Number of slices to process simultaneously.
- **`USE_GPU`**: Enable GPU computation (requires CuPy).
- **`WRITE_VECTORS`**: Save orientation 3rd vector as `.npy` files.

!!! warning

    Orientation vectors are saved in `float32` format and may consume significant disk space.

- **`REVERSE`**: Process volume from end to start if set to `True`.

!!! note

    The structure tensor calculation is performed using the [`structure-tensor`](https://github.com/Skielex/structure-tensor) Python package with optional CUDA support.

---

### `[ANGLE CALCULATION]`
- **`WRITE_ANGLES`**: Save helix and intrusion angles and fractional anisotropy values.
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

---
