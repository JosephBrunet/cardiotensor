"""
vector_field.py
---------------
High-level functions for loading, processing, and visualizing 3D vector fields
using FURY or exporting to VTK.
"""

import sys
import math
from pathlib import Path
import numpy as np

from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.downsampling import downsample_vector_volume, downsample_volume
from cardiotensor.visualization.fury_plotting_vectors import plot_vector_field_with_fury
from cardiotensor.utils.vector_vtk_export import export_vector_field_to_vtk


def visualize_vector_field(
    vector_field_path: str | Path,
    color_volume_path: str | Path | None = None,
    mask_path: str | Path | None = None,
    stride: int = 10,
    bin_factor: int = 1,
    size_arrow: float = 1.0,
    start: int | None = None,
    end: int | None = None,
    save_path: str | Path | None = None,
    voxel_size: float = 1.0,
    is_vtk: bool = False,
):
    """
    High-level visualization of a 3D vector field with FURY or optional VTK export.

    Parameters
    ----------
    vector_field_path : str or Path
        Path to the 3D vector field (directory or file).
    stride : int, optional
        Step size for downsampling during visualization. Default is 10.
    bin_factor : int, optional
        Spatial downsampling factor for the vector field. Default is 1.
    size_arrow : float, optional
        Scaling factor for the arrows in visualization. Default is 1.0.
    start, end : int, optional
        Z slice range to visualize. Default is full volume.
    save_path : str or Path, optional
        If provided, saves the visualization image.
    mask_path : str or Path, optional
        Optional mask volume to filter the vector field.
    color_path : str or Path, optional
        Optional HA volume for coloring the vector field.
    voxel_size : float, optional
        Voxel size for VTK export. Default is 1.0.
    is_vtk : bool, optional
        If True, exports the vector field to VTK.
    """

    vector_field_path = Path(vector_field_path)
    if not vector_field_path.exists():
        print(f"❌ Vector field path does not exist: {vector_field_path}")
        sys.exit(1)

    # Load input volume just for shape
    data_reader_vol = DataReader(vector_field_path)
    _, Z_full, Y_full, X_full = data_reader_vol.shape

    # If end slice is None, use the last slice
    start_idx = start or 0
    end_idx = end if end is not None else Z_full


    # Downsample if needed
    if bin_factor > 1:
        downsample_vector_volume(vector_field_path, bin_factor, vector_field_path.parent)
        vec_load_dir = vector_field_path.parent / f"bin{bin_factor}" / vector_field_path.name
        start_binned = start_idx // bin_factor
        end_binned = math.ceil(end_idx / bin_factor) if end_idx else None
    else:
        vec_load_dir = vector_field_path
        start_binned = start_idx
        end_binned = end_idx

    print(f"📥 Loading vector field from {vec_load_dir} ...")
    vec_reader = DataReader(vec_load_dir)
    vector_field = vec_reader.load_volume(start_index=start_binned, end_index=end_binned)

    # Ensure Z-component orientation
    print("🔄 Aligning vector orientations...")
    neg_mask = vector_field[0] > 0
    vector_field[:, neg_mask] *= -1
    del neg_mask

    # Convert from (3, Z, Y, X) → (Z, Y, X, 3)
    if vector_field.shape[0] == 3:
        vector_field = np.moveaxis(vector_field, 0, -1)

    # Optional mask
    if mask_path:
        print(f"🩹 Applying mask from {mask_path} ...")
        mask_reader = DataReader(mask_path)
        mask_volume = mask_reader.load_volume(
            start_index=start_binned,
            end_index=end_binned,
            unbinned_shape=vec_reader.shape[1:],  # (Z, Y, X)
        )
        mask = (mask_volume > 0).astype(np.uint8)
        vector_field[mask == 0, :] = np.nan

    # Optional for coloring  
    color_volume = None
    if color_volume_path:
        print(f"🎨 Loading color volume from {color_volume_path} ...")
        color_volume_path = Path(color_volume_path)
        if bin_factor > 1:
            downsample_volume(
                input_path=color_volume_path,
                bin_factor=bin_factor,
                output_dir=color_volume_path.parent,
                subfolder=color_volume_path.name,
                out_ext="tif",
                min_value=0,
                max_value=255,
            )
            color_load_dir = color_volume_path.parent / f"bin{bin_factor}" / color_volume_path.name
        else:
            color_load_dir = color_volume_path

        print(f"🎨 Loading color volume from {color_load_dir} ...")
        color_reader = DataReader(color_load_dir)
        color_volume = color_reader.load_volume(start_index=start_binned, end_index=end_binned)

    # Plot using FURY
    plot_vector_field_with_fury(
        vector_field,
        size_arrow=size_arrow,
        stride=stride,
        color_volume=color_volume,
        save_path=save_path,
    )

    # Optional VTK export
    if is_vtk:
        vtk_path = vector_field_path.parent / "paraview.vtk"
        print(f"💾 Exporting vector field to VTK at {vtk_path} ...")
        export_vector_field_to_vtk(
            vector_field=vector_field,
            color_volume=color_volume,
            voxel_size=voxel_size * bin_factor,
            stride=stride,
            save_path=vtk_path,
        )
