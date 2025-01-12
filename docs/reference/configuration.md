# Configuration File

The configuration file in **cardiotensor** allows you to specify the parameters required for processing input datasets, calculating orientation tensors, and saving results. Below is a detailed explanation of each section and parameter in the configuration file.

---
## Example Configuration File

This is an example of a configuration file as present in the `examples/` directory.

```ini
[DATASET]
# Path to the folder containing the input images (accepted formats: .tif, .jp2, .mhd)
IMAGES_PATH = ./data/635.2um_LADAF-2021-17_heart_overview_
# Voxel size of the input images in micrometers (µm).
VOXEL_SIZE = 635.2
# Path to the folder containing the segmentation mask (accepted formats: .tif or .jp2)
# If no mask is available, leave this field empty.
MASK_PATH = ./data/mask
# Whether to flip the image volume along the z-axis (True/False).
FLIP = True
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
# Whether to save the orientation vectors (as .npy) (True/False)
# Use for 3D vector/fiber visualisation
VECTORS = False

[STRUCTURE TENSOR CALCULATION]
# Noise scale window (sigma) used for structure tensor calculation.
# Adjust based on the level of noise in the dataset.
SIGMA = 1
# Integration scale window (rho) for structure tensor calculation.
# Larger values result in smoother orientation fields.
RHO = 2
# Number of slices to load into memory at a time during processing.
# This affects memory usage and processing speed. Adjust based on system capacity.
N_CHUNK = 20
# Enable GPU computation during the structure tensor calculation (True/False)
USE_GPU = False

[LV AXIS COORDINATES]
# Coordinates of the mitral valve point in the volume (X, Y, Z).
# This is required for defining the left ventricle axis.
POINT_MITRAL_VALVE = 104,110,116
# Coordinates of the apex point in the volume (X, Y, Z).
# This is required for defining the left ventricle axis.
POINT_APEX = 41,87,210

[RUN]
# Specify the processing direction:
#   - True: Process slices from the beginning (0) to the end.
#   - False: Process slices from the end to the beginning.
REVERSE = False

[TEST]
# Enable test mode:
#   - True: Process and plot only a single slice for testing.
#   - False: Perform the full processing on the entire volume.
TEST = True
# Specify the slice number to process when test mode is enabled.
N_SLICE_TEST = 155

```

!!! note

    Modify the configuration file as needed to fit your dataset.

---

## Explanation of Parameters

### `[DATASET]`
- **`IMAGES_PATH`**: Path to the input dataset containing 3D images.
- **`VOXEL_SIZE`**: Voxel size of the dataset in micrometers (e.g., 635.2 µm).
- **`MASK_PATH`**: Path to the binary segmentation mask. Leave blank if no mask is available.

!!! note

    The mask volume can be downscaled compare to the heart volume. The binning factor will be searched and the mask will be automatically upscaled to match the heart volume.

- **`FLIP`**: Boolean flag to flip the dataset along the z-axis.

### `[OUTPUT]`
- **`OUTPUT_PATH`**: Directory where processed results will be saved.
- **`OUTPUT_FORMAT`**: Format of output files (`jp2` or `tif`).
- **`OUTPUT_TYPE`**: Pixel format for the output:
    - `8bit`: Grayscale 8-bit.
    - `rgb`: 3-channel color.
- **`VECTORS`**: Boolean flag to save orientation vectors as `.npy` files for 3D visualization.

!!! warning

    The `.npy` vector files are saved in `float32` by default. The resulting disk size can be massive

### `[STRUCTURE TENSOR CALCULATION]`
- **`SIGMA`**: Noise scale for the structure tensor. Higher values smooth noisy data.
- **`RHO`**: Integration scale for the structure tensor. Larger values produce smoother orientation fields.

!!! note

    Adjust `SIGMA`, `RHO` based on dataset characteristics and voxel size.

- **`N_CHUNK`**: Number of slices to process at once. Adjust based on system memory.

!!! note

    Adjust `N_CHUNK` based on your system to avoid OOM error. If you obtain an OOM error try to reduce `N_CHUNK`.

- **`USE_GPU`**: Enable GPU computation during the structure tensor calculation (True/False)

!!! note

    The structure tensor calculation is performed using the [`structure-tensor`](https://github.com/Skielex/structure-tensor) python package. This package allows GPU computation using CUDA and [CuPy](https://github.com/cupy/cupy). 

### `[LV AXIS COORDINATES]`
- **`POINT_MITRAL_VALVE`**: (X, Y, Z) coordinates of the mitral valve.
- **`POINT_APEX`**: (X, Y, Z) coordinates of the apex. Both parameters are essential for defining the left ventricle axis.

### `[RUN]`
- **`REVERSE`**: Boolean to process slices in reverse order (end to beginning).

### `[TEST]`
- **`TEST`**: Boolean flag to enable test mode for processing a single slice.
- **`N_SLICE_TEST`**: Slice number to process in test mode.

!!! note

    Use test mode to quickly verify the pipeline before running full-volume processing.

---
