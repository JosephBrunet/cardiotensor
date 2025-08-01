from pathlib import Path

import numpy as np


def export_vector_field_to_vtk(
    vector_field: np.ndarray,
    color_volume: np.ndarray,
    voxel_size: float,
    stride: int = 32,
    save_path: Path | None = None,
) -> Path:
    if vector_field.ndim != 4 or vector_field.shape[-1] != 3:
        raise ValueError("vector_field must have shape (Z, Y, X, 3)")

    print("Preparing data for VTK export...")

    # Flip negative Z vectors
    vector_field[vector_field[..., 2] < 0] *= -1

    # Remove NaNs
    vector_field[~np.isfinite(vector_field)] = 0
    color_volume[~np.isfinite(color_volume)] = 0

    # Apply stride-based downsampling
    print(f"Downsampling vector field with stride={stride}...")

    vector_field = vector_field[::stride, ::stride, ::stride, :]
    color_volume = color_volume[::stride, ::stride, ::stride]

    # Create mask
    mask_volume = (np.linalg.norm(vector_field, axis=-1) > 0).astype(np.uint8)

    # Prepare output
    cellData = {
        "eigenVectors": vector_field,
        "color_angles": color_volume,
        "mask": mask_volume,
    }

    # Save to VTK
    if save_path is None:
        save_path = Path("paraview.vtk")
    print(f"Saving to VTK: {save_path}")
    writeStructuredVTK(
        aspectRatio=[voxel_size * stride] * 3,  # adjust spacing due to stride
        origin=[0.0, 0.0, 0.0],
        cellData=cellData,
        fileName=str(save_path),
    )
    print("Export complete.")
    return save_path


def writeStructuredVTK(
    aspectRatio: list[float] = [1.0, 1.0, 1.0],
    origin: list[float] = [0.0, 0.0, 0.0],
    cellData: dict[str, np.ndarray] = {},
    pointData: dict[str, np.ndarray] = {},
    fileName: str = "output.vtk",
) -> None:
    dimensions = []
    if not (cellData or pointData):
        print(f"No data provided to writeStructuredVTK() → {fileName}")
        return

    if cellData:
        shape = list(cellData.values())[0].shape[:3]
        for k, v in cellData.items():
            if v.shape[:3] != tuple(shape):
                raise ValueError(f"Inconsistent shape for field {k}")
        dimensions = [d + 1 for d in shape]

    with open(fileName, "w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write("VTK file from cardiotensor\n")
        f.write("ASCII\n\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write("DIMENSIONS {} {} {}\n".format(*reversed(dimensions)))
        f.write("ASPECT_RATIO {} {} {}\n".format(*reversed(aspectRatio)))
        f.write("ORIGIN {} {} {}\n\n".format(*reversed(origin)))

        if cellData:
            dims = shape
            f.write(f"CELL_DATA {dims[0] * dims[1] * dims[2]}\n\n")
            for name, data in cellData.items():
                _writeFieldInVtk({name: data}, f)


def _writeFieldInVtk(data: dict[str, np.ndarray], f, flat: bool = False) -> None:
    for key, field in data.items():
        if flat:
            raise NotImplementedError("Flat writing not supported here.")

        if field.ndim == 3:
            f.write(f"SCALARS {key} float\n")
            f.write("LOOKUP_TABLE default\n")
            for val in field.flatten():
                f.write(f"{val:.6f}\n")
        elif field.ndim == 4 and field.shape[3] == 3:
            f.write(f"VECTORS {key} float\n")
            for vec in field.reshape(-1, 3):
                f.write("{:.6f} {:.6f} {:.6f}\n".format(*reversed(vec)))
        else:
            raise ValueError(f"Unsupported shape for VTK field {key}: {field.shape}")
