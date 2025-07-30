#!/usr/bin/env python3
"""
generate_streamlines.py

Trace streamlines from a 3D vector field and save results (streamlines + HA) into .npz.

Usage:
    cardio-generate <path_to_conf.conf> [--start <int>] [--end <int>] [--bin <int>] [--seeds <int>]
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np

from cardiotensor.tractography.generate_streamlines import (
    generate_streamlines_from_vector_field,
)
from cardiotensor.utils.DataReader import DataReader
from cardiotensor.utils.downsampling import downsample_vector_volume, downsample_volume
from cardiotensor.utils.utils import read_conf_file



def script():
    parser = argparse.ArgumentParser(
        description="Trace streamlines from a 3D vector field and save to .npz"
    )
    parser.add_argument("conf_file", type=str, help="Path to .conf config file.")
    parser.add_argument("--start-z", type=int, default=0, help="Start index in Z.")
    parser.add_argument("--end-z", type=int, default=None, help="End index in Z.")
    parser.add_argument("--start-y", type=int, default=0, help="Start index in Y.")
    parser.add_argument("--end-y", type=int, default=None, help="End index in Y.")
    parser.add_argument("--start-x", type=int, default=0, help="Start index in X.")
    parser.add_argument("--end-x", type=int, default=None, help="End index in X.")
    parser.add_argument("--bin", type=int, default=1, help="Downsampling factor.")
    parser.add_argument("--seeds", type=int, default=20000, help="Number of seeds.")
    parser.add_argument(
        "--fa-threshold",
        type=float,
        default=0.1,
        help="Minimum FA to continue tracing.",
    )
    parser.add_argument("--step", type=float, default=1, help="Step length in voxels.")
    parser.add_argument(
        "--max_steps", type=int, default=None, help="Max steps (default no max)."
    )
    parser.add_argument("--angle", type=float, default=60.0, help="Angle threshold.")
    parser.add_argument(
        "--min_len", type=int, default=10, help="Min streamline length."
    )
    # parser.add_argument("--save-trk", action="store_true", help="Also save streamlines in TrackVis .trk format")

    args = parser.parse_args()

    bin_factor = args.bin
    num_seeds = args.seeds

    # Load configuration
    conf_path = Path(args.conf_file)
    try:
        params = read_conf_file(conf_path)
    except Exception as e:
        print(f"⚠️ Error reading config: {e}")
        sys.exit(1)

    VOLUME_PATH = params.get("IMAGES_PATH", "")
    OUTPUT_DIR = Path(params.get("OUTPUT_PATH", "./output"))
    MASK_PATH = params.get("MASK_PATH", "")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load input volume info
    vec_reader = DataReader(OUTPUT_DIR / "eigen_vec")
    full_shape = vec_reader.shape  # (3, Z, Y, X)

    start_z = args.start_z
    end_z = args.end_z or full_shape[1]
    start_y = args.start_y
    end_y = args.end_y or full_shape[2]
    start_x = args.start_x
    end_x = args.end_x or full_shape[3]

    # Load input volume just for shape
    data_reader_vol = DataReader(VOLUME_PATH)
    Z_full, Y_full, X_full = data_reader_vol.shape
    if end_z is None:
        end_z = Z_full

    # Setup eigenvector volume
    eigen_dir = OUTPUT_DIR / "eigen_vec"
    if not eigen_dir.exists():
        print(f"⚠️ No eigenvector directory at {eigen_dir}")
        sys.exit(1)


    #---------------------------------
    # BINNING
    
    # Manage indexes if bin_factor > 1
    if bin_factor > 1:
        downsample_vector_volume(eigen_dir, bin_factor, OUTPUT_DIR)
        vec_load_dir = OUTPUT_DIR / f"bin{bin_factor}" / "eigen_vec"
        start_z_binned = start_z // bin_factor
        end_z_binned = math.ceil(end_z / bin_factor)
        
        start_y_binned = start_y // bin_factor
        end_y_binned = math.ceil(end_y / bin_factor)

        start_x_binned = start_x // bin_factor
        end_x_binned = math.ceil(end_x / bin_factor)

    else:
        vec_load_dir = eigen_dir
        start_z_binned = start_z
        end_z_binned = end_z
        start_y_binned = start_y
        end_y_binned = end_y
        start_x_binned = start_x
        end_x_binned = end_x
        

    FA_dir = OUTPUT_DIR / "FA"
    if bin_factor > 1:
        downsample_volume(
            input_path=FA_dir,
            bin_factor=bin_factor,
            output_dir=OUTPUT_DIR,
            subfolder="FA",
            out_ext="tif",
            min_value=0,
            max_value=255,
        )
        fa_load_dir = OUTPUT_DIR / f"bin{bin_factor}" / "FA"
    else:
        fa_load_dir = FA_dir
        
    HA_dir = OUTPUT_DIR / "HA"
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
        
    #---------------------------------


    print("Loading vector field...")
    vec_reader = DataReader(vec_load_dir)
        
    vector_field = vec_reader.load_volume(
        start_index=start_z_binned, end_index=end_z_binned
    )[:, :, start_y_binned:end_y_binned, start_x_binned:end_x_binned]

    if vector_field.ndim == 4 and vector_field.shape[0] == 3:
        vector_field = vector_field
    elif vector_field.ndim == 4 and vector_field.shape[-1] == 3:
        print("⚠️ Vector field is in (Z, Y, X, 3) format, moving channels to first dimension...")
        vector_field = np.moveaxis(vector_field, -1, 0)
    else:
        print(f"⚠️ Unexpected vector field shape: {vector_field.shape}")
        sys.exit(1)

    print("Ensuring Z-components are positive...")
    neg_mask = vector_field[0] < 0
    vector_field[:, neg_mask] *= -1
    del neg_mask

    if MASK_PATH:
        print("Applying mask from config...")
        
        print(f"MASK_PATH: {MASK_PATH}\nstart_binned: {start_z_binned}, end_binned: {end_z_binned}", vec_reader.shape[1:])
        
        mask_reader = DataReader(MASK_PATH)

        # Load the corresponding mask volume, resampled to match vector field shape
        mask = mask_reader.load_volume(
            start_index=start_z_binned,
            end_index=end_z_binned,
            unbinned_shape=vec_reader.shape[1:],  # (Z, Y, X)
        )
        mask = mask[:, start_y_binned:end_y_binned, start_x_binned:end_x_binned]

        mask = (mask > 0).astype(np.uint8)
        vector_field[:, mask == 0] = np.nan

    # Seed mask from FA
    print("Generating seed mask using FA > 0.4...")

    if not fa_load_dir.exists():
        print(f"⚠️ No FA directory found at {fa_load_dir}")
        sys.exit(1)

    print("Loading FA volume...")
    fa_reader = DataReader(fa_load_dir)
    # FA vlume is stored as 8 bits (0-255), so we need to scale it
    fa_volume = fa_reader.load_volume(start_index=start_z_binned, end_index=end_z_binned)    
    fa_volume = fa_volume[:, start_y_binned:end_y_binned, start_x_binned:end_x_binned]


    seed_mask = (fa_volume > (0.25 * 255))

    print(f"Selecting {num_seeds} random seed points from mask...")
    valid_indices = np.argwhere(seed_mask > 0)
    if valid_indices.shape[0] < num_seeds:
        print("⚠️ Not enough valid seed points.")
        sys.exit(1)

    chosen_indices = valid_indices[
        np.random.choice(valid_indices.shape[0], num_seeds, replace=False)
    ]
    del valid_indices

    # Trace streamlines
    streamlines = generate_streamlines_from_vector_field(
        vector_field=vector_field,
        seed_points=chosen_indices,
        fa_volume=fa_volume,
        fa_threshold=args.fa_threshold,
        step_length=args.step,
        max_steps=args.max_steps,
        angle_threshold=args.angle,
        min_length_pts=args.min_len,
    )

    # Load HA for sampling
    if not HA_dir.exists():
        print(f"⚠️ No HA directory found at {HA_dir}")
        sys.exit(1)

    ha_reader = DataReader(ha_load_dir)
    HA_volume = ha_reader.load_volume(start_index=start_z_binned, end_index=end_z_binned)
    HA_volume = HA_volume[:, start_y_binned:end_y_binned, start_x_binned:end_x_binned]

    def sample_ha_along(pts_list):
        ha_vals = []
        Z, Y, X = HA_volume.shape
        for z, y, x in pts_list:
            zi, yi, xi = map(int, [round(z), round(y), round(x)])
            zi = max(0, min(zi, Z - 1))
            yi = max(0, min(yi, Y - 1))
            xi = max(0, min(xi, X - 1))
            ha_vals.append(float(HA_volume[zi, yi, xi]))
        return ha_vals

    print("Sampling HA along streamlines...")
    all_ha = []
    for sl in streamlines:
        all_ha.extend(sample_ha_along(sl))

    # Save output
    out_path = OUTPUT_DIR / "streamlines.npz"
    np.savez_compressed(
        out_path,
        streamlines=np.array(streamlines, dtype=object),
        ha_values=np.array(all_ha, dtype=np.float32),
    )

    print(
        f"✅ Saved {len(streamlines)} streamlines and {len(all_ha)} HA values to {out_path}"
    )


if __name__ == "__main__":
    script()
