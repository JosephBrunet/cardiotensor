import argparse
import math
import sys
from pathlib import Path

import numpy as np

from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.downsampling import downsample_vector_volume, downsample_volume
from cardiotensor.utils.plot_vector_field import plot_vector_field_with_fury
from cardiotensor.utils.utils import read_conf_file



def script():
    parser = argparse.ArgumentParser(
        description="Plot 3D vector field using FURY from configuration file."
    )
    parser.add_argument("conf_file", type=Path, help="Path to configuration file")
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Stride to downsample vectors for display (default: 10)",
    )
    parser.add_argument(
        "--bin",
        type=int,
        default=1,
        help="Binning factor used during preprocessing (default: 1)",
    )
    parser.add_argument(
        "--size-arrow",
        type=float,
        default=1,
        help="Factor used for arrow size (default: 1)",
    )
    parser.add_argument("--start", type=int, default=None, help="Start slice index")
    parser.add_argument("--end", type=int, default=None, help="End slice index")
    parser.add_argument(
        "--save", type=Path, help="Optional path to save rendered image"
    )
    parser.add_argument(
        "--vtk", action="store_true", help="Export the vector field to VTK for ParaView"
    )

    args = parser.parse_args()

    # Read parameters
    params = read_conf_file(args.conf_file)
    VOLUME_PATH = params.get("IMAGES_PATH", "")
    OUTPUT_DIR = Path(params.get("OUTPUT_PATH", "./output"))
    VOXEL_SIZE = params.get("VOXEL_SIZE", 1)
    MASK_PATH = params.get("MASK_PATH", "")

    OUTPUT_DIR = Path(params.get("OUTPUT_PATH", "./output"))

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Slice range
    start_idx = args.start if args.start is not None else 0
    end_idx = args.end

    # Load input volume just for shape
    data_reader_vol = DataReader(VOLUME_PATH)
    Z_full, Y_full, X_full = data_reader_vol.shape
    if end_idx is None:
        end_idx = Z_full

    # Setup eigenvector volume
    eigen_dir = OUTPUT_DIR / "eigen_vec"
    if not eigen_dir.exists():
        print(f"⚠️ No eigenvector directory at {eigen_dir}")
        sys.exit(1)

    bin_factor = args.bin
    if bin_factor > 1:
        downsample_vector_volume(eigen_dir, bin_factor, OUTPUT_DIR)
        vec_load_dir = OUTPUT_DIR / f"bin{bin_factor}" / "eigen_vec"
        start_binned = start_idx // bin_factor
        end_binned = math.ceil(end_idx / bin_factor)
    else:
        vec_load_dir = eigen_dir
        start_binned = start_idx
        end_binned = end_idx

    print("Loading vector field...")
    vec_reader = DataReader(vec_load_dir)
    vector_field = vec_reader.load_volume(
        start_index=start_binned, end_index=end_binned
    )
        
    print("Ensuring Z-components are positive...")
    neg_mask = vector_field[0] > 0  # Identify where Z component is negative
    vector_field[:, neg_mask] *= -1  # Flip the entire vector at that location
    del neg_mask
    
    

    # If your vector_field is in shape (3, Z, Y, X), convert it:
    if vector_field.shape[0] == 3:
        vector_field = np.moveaxis(vector_field, 0, -1)

    if MASK_PATH:
        print("Applying mask from config...")
        mask_reader = DataReader(MASK_PATH)

        # Load the corresponding mask volume, resampled to match vector field shape
        mask_volume = mask_reader.load_volume(
            start_index=start_binned,
            end_index=end_binned,
            unbinned_shape=vec_reader.shape[1:],  # (Z, Y, X)
        )

        mask = (mask_volume > 0).astype(np.uint8)

        vector_field[mask == 0, :] = np.nan

    # Load HA for sampling
    HA_dir = OUTPUT_DIR / "HA"
    if not HA_dir.exists():
        print(f"⚠️ No HA directory found at {HA_dir}")
        sys.exit(1)

    if bin_factor > 1:
        downsample_volume(
            input_path=HA_dir,
            bin_factor=bin_factor,
            output_dir=OUTPUT_DIR,
            subfolder="HA",
            out_ext="tif",
            min_value=0,
            max_value=255,
        )
        ha_load_dir = OUTPUT_DIR / f"bin{bin_factor}" / "HA"
    else:
        ha_load_dir = HA_dir

    ha_reader = DataReader(ha_load_dir)
    HA_volume = ha_reader.load_volume(start_index=start_binned, end_index=end_binned)

    plot_vector_field_with_fury(
        vector_field,
        size_arrow=args.size_arrow,
        stride=args.stride,
        ha_volume=HA_volume,
        save_path=args.save,
    )

    if args.vtk:
        print("💾 Exporting vector field to VTK...")
        from cardiotensor.utils.vector_vtk_export import export_vector_field_to_vtk

        vtk_path = OUTPUT_DIR / "paraview.vtk"
        export_vector_field_to_vtk(
            vector_field=vector_field,
            HA_volume=HA_volume,
            voxel_size=VOXEL_SIZE * bin_factor,
            stride=args.stride,
            save_path=vtk_path,
        )
