# Python API

Cardiotensor can be used as a Python library in addition to the CLI.
The main workflow helpers are importable directly from `cardiotensor`, and
subpackage imports are available for grouped functionality.

## Installation

```console
$ pip install cardiotensor
```

---

## Quick Start

```python
from cardiotensor import compute_orientation

compute_orientation(
    volume_path="path/to/images/",
    mask_path="path/to/mask.tif",       # optional
    output_dir="./output",
    sigma=1.0,
    rho=3.0,
    axis_points=[(0, 0, 0), (256, 256, 256), (512, 512, 512)],
)
```

---

## Structure Tensor

```python
import numpy as np
from cardiotensor import calculate_structure_tensor, compute_fraction_anisotropy

volume = np.random.rand(32, 64, 64).astype(np.float32)

# Compute eigenvalues and eigenvectors
eigenvalues, eigenvectors = calculate_structure_tensor(
    volume,
    sigma=1.0,
    rho=3.0,
    use_gpu=False,
)
# eigenvectors shape: (3, Z, Y, X)
# eigenvalues  shape: (3, Z, Y, X)

# Fractional Anisotropy map
z_index = volume.shape[0] // 2
fa = compute_fraction_anisotropy(eigenvalues[:, z_index, :, :])
print(fa.shape)   # (Y, X), values in [0, 1]
```

---

## Helix and Transverse Angles

```python
import numpy as np
import matplotlib.pyplot as plt
from cardiotensor import compute_helix_and_transverse_angles
from cardiotensor.orientation import rotate_vectors_to_new_axis

# Use one slice of the eigenvector field computed above
vector_field_slice = eigenvectors[:, z_index, :, :]

# Example long-axis direction used to rotate the slice before angle computation
new_axis_vec = np.array([0.0, 1.0, 0.0], dtype=np.float32)
rotated_vectors = rotate_vectors_to_new_axis(vector_field_slice, new_axis_vec)

center_point = (
    vector_field_slice.shape[2] // 2,
    vector_field_slice.shape[1] // 2,
    z_index,
)

helix_angle, intrusion_angle = compute_helix_and_transverse_angles(
    rotated_vectors,
    center_point,
)
# Both arrays shape: (Y, X), values in degrees
# The second output corresponds to the transverse/intrusion companion angle map

fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)

im0 = axes[0].imshow(helix_angle, cmap="twilight", vmin=-90, vmax=90)
axes[0].set_title("Helix angle")
axes[0].axis("off")
fig.colorbar(im0, ax=axes[0], label="degrees")

im1 = axes[1].imshow(intrusion_angle, cmap="coolwarm", vmin=-90, vmax=90)
axes[1].set_title("Intrusion angle")
axes[1].axis("off")
fig.colorbar(im1, ax=axes[1], label="degrees")

plt.show()
```

---

## Available Imports

```python
from cardiotensor import (
    # Core
    DataReader,
    compute_orientation,
    read_conf_file,
    convert_to_8bit,
    # Orientation
    calculate_structure_tensor,
    compute_azimuth_and_elevation,
    compute_fraction_anisotropy,
    compute_helix_and_transverse_angles,
    # Tractography
    generate_streamlines_from_vector_field,
    generate_streamlines_from_params,
    # I/O
    load_npz_streamlines,
    load_trk_streamlines,
    write_spatialgraph_am,
    export_vector_field_to_vtk,
    # Analysis
    calculate_intensities,
    find_end_points,
    plot_intensity,
    save_intensity,
)
```

You can also use grouped imports:

```python
from cardiotensor.analysis import calculate_intensities, find_end_points
from cardiotensor.orientation import (
    calculate_structure_tensor,
    compute_orientation,
    rotate_vectors_to_new_axis,
)
from cardiotensor.tractography import generate_streamlines_from_params
from cardiotensor.utils import DataReader, load_trk_streamlines, read_conf_file
from cardiotensor.visualization import visualize_streamlines, visualize_vector_field
```


!!! note

    For the full API reference with parameter descriptions, see the [API reference](../reference/api.md) page.
