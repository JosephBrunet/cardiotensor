### index.md

# CardioTensor

**CardioTensor** is a Python package designed to compute and analyze 3D cardiomyocyte orientations in high-resolution heart images. This documentation provides a detailed guide on installation, usage, and contributing to the project.

## Features

- Automated orientation analysis
- Configurable workflows with `.conf` files
- Export support for visualization tools like Amira and VTK
- GPU support for faster processing

---

### installation.md

# Installation

## Prerequisites

- Python 3.11 or newer
- Required libraries (see `pyproject.toml`)

## Installation Steps

1. Clone the repository:

   ```bash
   git clone https://github.com/JosephBrunet/cardiotensor.git
   cd cardiotensor
   ```

2. Install the package:

   ```bash
   pip install .
   ```

3. Verify installation:

   ```bash
   cardio-tensor --help
   ```

---

### usage/index.md

# Usage

This section explains the primary functionalities of the CardioTensor package.

## Quick Start

1. Prepare a configuration file (e.g., `parameters_example.conf`).
2. Run the `cardio-tensor` script:
   ```bash
   cardio-tensor parameters_example.conf
   ```
3. Visualize results using the GUI or export them using `cardio-vtk` or `cardio-amira`.

---

### usage/scripts.md

# Command-Line Scripts

## `cardio-tensor`

Computes cardiomyocyte orientations.

```bash
cardio-tensor parameters_example.conf --gpu
```

### Options

- `--start_index` / `--end_index`: Subset of the volume to process
- `--gpu`: Use GPU for faster processing

## `cardio-analysis`

Launches the GUI for analyzing slices.

```bash
cardio-analysis parameters_example.conf 150
```

### Options

- `--N_line`: Number of profile lines
- `--angle_range`: Angle range in degrees
- `--image_mode`: Output mode (HA, IA, FA)

## `cardio-vtk`

Exports results in VTK format.

```bash
cardio-vtk parameters_example.conf --start_index 0 --end_index 100
```

## `cardio-amira`

Exports results in Amira-compatible format.

```bash
cardio-amira parameters_example.conf --num_steps 1000000
```

---

### usage/configuration.md

# Configuration File

Configuration files control the behavior of the package. Below is an example:

```ini
[IMAGES_PATH]
VOLUME_PATH = ./data/volume/
MASK_PATH = ./data/mask/

[PROCESSING]
FLIP = False
OUTPUT_PATH = ./output/
OUTPUT_TYPE = 8bit
SIGMA = 1.0
RHO = 2.0
N_CHUNK = 100
POINT_MITRAL_VALVE = (50, 50, 50)
POINT_APEX = (150, 150, 150)
TEST = False
```

---

### advanced.md

# Advanced Features

## GPU Support

Ensure `cupy-cuda` is installed to enable GPU support.

## Export Options

- **VTK**: Use `cardio-vtk` for exporting results to VTK.
- **Amira**: Use `cardio-amira` with customizable parameters (e.g., fiber length, angle thresholds).

---

### development.md

# Development Guide

## Testing

Run tests using `pytest`:

```bash
pytest
```

## Contributing

1. Fork the repository.
2. Create a feature branch.
3. Submit a pull request.

Follow PEP8 guidelines and ensure tests pass before submission.

---

### reference/index.md

# API Reference

Refer to the individual modules for detailed information about the functions and classes provided.

---

### reference/orientation\_computation.md

# `orientation_computation_functions.py`

### Functions

- `interpolate_points`: Generates interpolated points between two 3D points.
- `calculate_structure_tensor`: Computes the structure tensor for a 3D volume.
- `compute_fraction_anisotropy`: Computes Fractional Anisotropy from eigenvalues.

### Example

```python
from cardiotensor.orientation.orientation_computation_functions import (
    calculate_structure_tensor
)

tensor, eigenvalues, eigenvectors = calculate_structure_tensor(volume, sigma=1.0, rho=2.0)
```

---

This structure divides the documentation logically, providing clarity and accessibility for different user groups. Copy each file's content into a respective `.md` file for your MkDocs documentation.

