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
