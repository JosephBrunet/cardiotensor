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
from cardiotensor.visualization.fury_plotting import plot_vector_field_with_fury
from cardiotensor.utils.vector_vtk_export import export_vector_field_to_vtk


def visualize_vector_field(
    vector_field_path: str | Path,
    stride: int = 10,
    bin_factor: int = 1,
    size_arrow: float = 1.0,
    start: int | None = None,
    end: int | None = None,
    save_path: str | Path | None = None,
    mask_path: str | Path | None = None,
    ha_path: str | Path | None = None,
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
    ha_path : str or Path, optional
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

    # Optional HA for coloring
    ha_volume = None
    if ha_path:
        print(f"🎨 Loading HA volume from {ha_path} ...")
        ha_path = Path(ha_path)
        if bin_factor > 1:
            downsample_volume(
                input_path=ha_path,
                bin_factor=bin_factor,
                output_dir=ha_path.parent,
                subfolder=ha_path.name,
                out_ext="tif",
                min_value=0,
                max_value=255,
            )
            ha_load_dir = ha_path.parent / f"bin{bin_factor}" / ha_path.name
        else:
            ha_load_dir = ha_path

        ha_reader = DataReader(ha_load_dir)
        ha_volume = ha_reader.load_volume(start_index=start_binned, end_index=end_binned)

    # Plot using FURY
    plot_vector_field_with_fury(
        vector_field,
        size_arrow=size_arrow,
        stride=stride,
        ha_volume=ha_volume,
        save_path=save_path,
    )

    # Optional VTK export
    if is_vtk:
        vtk_path = vector_field_path.parent / "paraview.vtk"
        print(f"💾 Exporting vector field to VTK at {vtk_path} ...")
        export_vector_field_to_vtk(
            vector_field=vector_field,
            HA_volume=ha_volume,
            voxel_size=voxel_size * bin_factor,
            stride=stride,
            save_path=vtk_path,
        )

















    # # Slice range
    # start_idx = start if start is not None else 0
    # end_idx = end

    # # Load input volume just for shape
    # data_reader_vol = DataReader(VOLUME_PATH)
    # Z_full, Y_full, X_full = data_reader_vol.shape
    # if end_idx is None:
    #     end_idx = Z_full

    # # Setup eigenvector volume
    # eigen_dir = OUTPUT_DIR / "eigen_vec"
    # if not eigen_dir.exists():
    #     print(f"⚠️ No eigenvector directory at {eigen_dir}")
    #     sys.exit(1)

    # bin_factor = args.bin
    # if bin_factor > 1:
    #     downsample_vector_volume(eigen_dir, bin_factor, OUTPUT_DIR)
    #     vec_load_dir = OUTPUT_DIR / f"bin{bin_factor}" / "eigen_vec"
    #     start_binned = start_idx // bin_factor
    #     end_binned = math.ceil(end_idx / bin_factor)
    # else:
    #     vec_load_dir = eigen_dir
    #     start_binned = start_idx
    #     end_binned = end_idx

    # print("Loading vector field...")
    # vec_reader = DataReader(vec_load_dir)
    # vector_field = vec_reader.load_volume(
    #     start_index=start_binned, end_index=end_binned
    # )
        
    # print("Ensuring Z-components are positive...")
    # neg_mask = vector_field[0] > 0  # Identify where Z component is negative
    # vector_field[:, neg_mask] *= -1  # Flip the entire vector at that location
    # del neg_mask
    
    

    # # If your vector_field is in shape (3, Z, Y, X), convert it:
    # if vector_field.shape[0] == 3:
    #     vector_field = np.moveaxis(vector_field, 0, -1)

    # if MASK_PATH:
    #     print("Applying mask from config...")
    #     mask_reader = DataReader(MASK_PATH)

    #     # Load the corresponding mask volume, resampled to match vector field shape
    #     mask_volume = mask_reader.load_volume(
    #         start_index=start_binned,
    #         end_index=end_binned,
    #         unbinned_shape=vec_reader.shape[1:],  # (Z, Y, X)
    #     )

    #     mask = (mask_volume > 0).astype(np.uint8)

    #     vector_field[mask == 0, :] = np.nan

    # # Load HA for sampling
    # HA_dir = OUTPUT_DIR / "HA"
    # if not HA_dir.exists():
    #     print(f"⚠️ No HA directory found at {HA_dir}")
    #     sys.exit(1)

    # if bin_factor > 1:
    #     downsample_volume(
    #         input_path=HA_dir,
    #         bin_factor=bin_factor,
    #         output_dir=OUTPUT_DIR,
    #         subfolder="HA",
    #         out_ext="tif",
    #         min_value=0,
    #         max_value=255,
    #     )
    #     ha_load_dir = OUTPUT_DIR / f"bin{bin_factor}" / "HA"
    # else:
    #     ha_load_dir = HA_dir

    # ha_reader = DataReader(ha_load_dir)
    # HA_volume = ha_reader.load_volume(start_index=start_binned, end_index=end_binned)

    # plot_vector_field_with_fury(
    #     vector_field,
    #     size_arrow=args.size_arrow,
    #     stride=args.stride,
    #     ha_volume=HA_volume,
    #     save_path=args.save,
    # )

    # if args.vtk:
    #     print("💾 Exporting vector field to VTK...")
    #     from cardiotensor.utils.vector_vtk_export import export_vector_field_to_vtk

    #     vtk_path = OUTPUT_DIR / "paraview.vtk"
    #     export_vector_field_to_vtk(
    #         vector_field=vector_field,
    #         HA_volume=HA_volume,
    #         voxel_size=VOXEL_SIZE * bin_factor,
    #         stride=args.stride,
    #         save_path=vtk_path,
    #     )
