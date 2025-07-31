#!/usr/bin/env python3
"""
generate_streamlines.py
-----------------------
CLI script to generate streamlines from a vector field.
Reads a .conf file for paths and parameters, then calls the library function.
"""

import argparse
import sys
from pathlib import Path
from cardiotensor.utils.utils import read_conf_file
from cardiotensor.tractography.generate_streamlines import generate_streamlines_from_params


def script():
    parser = argparse.ArgumentParser(
        description="Generate streamlines from a 3D vector field and save to .npz"
    )
    parser.add_argument("conf_file", type=Path, help="Path to configuration .conf file.")
    parser.add_argument("--start-z", type=int, default=0, help="Start slice index in Z")
    parser.add_argument("--end-z", type=int, default=None, help="End slice index in Z")
    parser.add_argument("--start-y", type=int, default=0, help="Start slice index in Y")
    parser.add_argument("--end-y", type=int, default=None, help="End slice index in Y")
    parser.add_argument("--start-x", type=int, default=0, help="Start slice index in X")
    parser.add_argument("--end-x", type=int, default=None, help="End slice index in X")
    parser.add_argument("--bin", type=int, default=1, help="Downsampling factor")
    parser.add_argument("--seeds", type=int, default=20000, help="Number of seeds")
    parser.add_argument("--fa-threshold", type=float, default=0.1, help="FA threshold")
    parser.add_argument("--step", type=float, default=0.5, help="Step length in voxels")
    parser.add_argument("--max-steps", type=int, default=None, help="Max steps per streamline")
    parser.add_argument("--angle", type=float, default=60.0, help="Max turning angle (deg)")
    parser.add_argument("--min-len", type=int, default=10, help="Minimum streamline length (points)")

    args = parser.parse_args()

    # Load configuration
    try:
        params = read_conf_file(args.conf_file)
    except Exception as e:
        print(f"❌ Failed to read config file: {e}")
        sys.exit(1)

    output_dir = Path(params.get("OUTPUT_PATH", "./output"))
    vector_field_dir = output_dir / "eigen_vec"
    fa_dir = output_dir / "FA"
    ha_dir = output_dir / "HA"
    mask_path = params.get("MASK_PATH", None)

    # Pass paths and params to the library function
    generate_streamlines_from_params(
        vector_field_dir=vector_field_dir,
        output_dir=output_dir,
        fa_dir=fa_dir,
        ha_dir=ha_dir,
        mask_path=mask_path,
        start_xyz=(args.start_z, args.start_y, args.start_x),
        end_xyz=(args.end_z, args.end_y, args.end_x),
        bin_factor=args.bin,
        num_seeds=args.seeds,
        fa_threshold=args.fa_threshold,
        step_length=args.step,
        max_steps=args.max_steps,
        angle_threshold=args.angle,
        min_length_pts=args.min_len,
    )


if __name__ == "__main__":
    script()
